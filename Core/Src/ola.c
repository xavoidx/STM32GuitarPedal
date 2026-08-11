#include "ola.h"
#include <float.h>

// Normalized correlation score. ||target|| is constant across lags and cancels,
// so only the candidate energy has to be divided out. The floor guards a silent
// candidate (0/0) and the float cancellation you get differencing two large
// prefix sums when a loud attack shares the cache with a near-silent window.
static inline float ncc(float cross, float energy, float energy_floor) {
    return cross / sqrtf(energy > energy_floor ? energy : energy_floor);
}

// NOTE ON HYSTERESIS -- deliberately absent, do not "restore" it.
//
// There used to be a bias weighting each correlation by how close its lag was to
// the previous frame's chosen lag, on the theory that it damped lag jitter. It
// does the opposite. The nominal analysis position advances by OLA_SA while the
// output advances by OLA_SS, so waveform continuity requires
//     lag_k = lag_{k-1} + (OLA_SS - OLA_SA) - m*P
// for pitch period P and some integer m. The correct lag WALKS by 256 mod P every
// frame; it does not sit still. Biasing toward prev_lag therefore pulls away from
// the phase-correct choice, and a bias big enough to matter is big enough to
// override the true correlation peak. Measured cost when it was in: ~23 dB of
// extra artifact energy and several cents of pitch error on steady tones.
//
// A small bias toward lag == 0 would be defensible -- all period-equivalent lags
// splice equally well, so the one nearest the nominal position minimizes envelope
// mismatch -- but that is a different quantity from prev_lag, and steady-tone
// tests cannot show whether it helps. It would need transient material to judge.

void OLA_Init(OLA_State *s) {
    memset(s, 0, sizeof(*s));               // FIX: sizeof(*s), not sizeof(s) (pointer)
    for(size_t i = 0; i < OLA_W; ++i) {
        // FIX: standard Hann is 0.5 - 0.5cos (zero at edges, 1 at center).
        // The '+' form is an inverted window and destroys the signal.
        //
        // Divisor is OLA_W, not OLA_W - 1: the PERIODIC Hann satisfies
        // h[i] + h[i + W/2] == 1 exactly, so overlap-add at a hop of OLA_SS is
        // exactly flat. The symmetric form leaves ~0.13% COLA ripple.
        s->hann[i] = 0.5f - 0.5f * cosf((2.0f * M_PI * i) / OLA_W);
    }
    // Prime the read pointer to trail the write frontier by one synthesis hop.
    // See the note on the OLA_SS assert in ola.h for why this is safe.
    s->out_read = OLA_OUT_SIZE - OLA_SS;
}

__attribute__((optimize("fast-math")))
void OLA_Process(OLA_State* s, const float* input_buffer, float* output_buffer, size_t buffer_size) {

    // Push the new block into the input ring buffer, DC-blocked on the way in.
    //
    // Why this matters here specifically: the correlation below is a RAW dot
    // product, so two signals each carrying a DC offset mu contribute a term
    // ~ OLA_L * mu_cand * mu_tail to every lag. While the note is loud that term
    // is buried under the AC energy. As the note decays toward the ADC noise
    // floor it starts to dominate, the correlation peak flattens, and the chosen
    // lag goes essentially random frame-to-frame -> warble that gets worse the
    // quieter the note. y[n] = x[n] - x[n-1] + R*y[n-1] costs 2 mul-adds/sample.
    for(size_t i = 0; i < buffer_size; ++i) {
        const float xn = input_buffer[i];
        const float yn = xn - s->dc_x1 + OLA_DC_R * s->dc_y1;
        s->dc_x1 = xn;
        s->dc_y1 = yn;

        s->in_buffer[s->in_write] = yn;
        s->in_write++;
        if(s->in_write >= OLA_IN_SIZE) s->in_write = 0;
    }
    s->in_count += buffer_size;  
  
    // Fire one analysis frame for every SA samples accumulated.
    while(s->in_count >= OLA_SA) {

        // #4: hoist the read-only input ring into a restrict-qualified local.
        //     in_buffer is only WRITTEN in the push loop above, never inside
        //     this search, so restrict is valid and lets the compiler assume
        //     it never aliases the out_buffer writes.
        const float * restrict in = s->in_buffer;
        const size_t out_write = s->out_write;

        // Logical newest sample (exclusive) available to THIS frame. Equals
        // in_write whenever the caller passes blocks of exactly OLA_SA, which is
        // what the pedal does. It differs only if one call carries several
        // frames' worth of input, and then it is what keeps consecutive frames
        // exactly OLA_SA apart instead of letting them share one anchor -- which
        // would collapse the analysis hop and with it the 2x stretch ratio.
        // All the causality bounds below are stated against this, not in_write.
        const size_t frame_end = s->in_write - (s->in_count - OLA_SA);

        // ---- Correlation target: the natural continuation of the segment the
        // PREVIOUS frame used, read straight from the input ring.
        //
        // The previous frame wrote in[prev_start + i] * hann[i] at
        // out_write_prev + i, so out_buffer at out_write holds exactly
        // in[prev_start + OLA_SS + j] * hann[OLA_SS + j]. Correlating against the
        // output tail (what this used to do) is therefore the same data with a
        // Hann taper falling 1.0 -> 0.02 across the window, which collapsed the
        // effective correlation length to ~200 samples -- under half a period of
        // low E -- and hard-capped OLA_L at OLA_SS. Reading the raw input gives a
        // flat window and no cap.
        //
        // OLA_L now exceeds the OLA_SS overlap, so the tail of this correlation
        // looks past the seam. That is deliberate: it is what disambiguates the
        // pitch period when several period-equivalent lags sit in the search range.
        static float target[OLA_L];
        float t_energy = 0.0f;
        for(size_t j = 0; j < OLA_L; ++j) {
            const float t = in[(s->prev_start + OLA_SS + j) % OLA_IN_SIZE];
            target[j] = t;
            t_energy += t * t;
        }

        /**
         * Cache the searchable span of input, flattened so lag arithmetic is a
         * plain array offset.
         *
         * Frame start for lag L is  frame_end - OLA_ANCHOR + L,  L in [-LAG, +LAG],
         * so the candidate that begins at cache index (L + OLA_LAG) is exactly
         * that frame. The span runs
         *   [frame_end - (OLA_W + 2*OLA_LAG), frame_end - (OLA_W - OLA_L))
         * whose newest sample is frame_end - 384 -- always in the past.
         *
         * NOTE: only the FIRST OLA_L samples of a candidate are ever correlated.
         * The rest of the frame is copied unverified, so the causality bound has
         * to be enforced by construction rather than discovered by the search.
         */
        static float input_cache[2 * OLA_LAG + OLA_L];
        for(int i = 0; i < 2 * OLA_LAG + OLA_L; ++i) {
            input_cache[i] = in[(frame_end + i - OLA_W - 2 * OLA_LAG + 2 * OLA_IN_SIZE) % OLA_IN_SIZE];
        }

        // ---- Candidate energies for the normalized correlation.
        //
        // A raw dot product scales with ||candidate||, so among the several
        // period-equivalent lags in a +-OLA_LAG window it systematically picks
        // the loudest -- the earliest one, on a decaying note. That crossfades a
        // louder segment onto a decaying tail every frame, i.e. an amplitude
        // sawtooth at the 187.5 Hz frame rate. ||target|| is still constant
        // across lags and still cancels, so only the candidate energy is needed,
        // and a prefix sum makes it O(1) per lag.
        static float psq[2 * OLA_LAG + OLA_L + 1];
        psq[0] = 0.0f;
        for(int i = 0; i < 2 * OLA_LAG + OLA_L; ++i) {
            psq[i + 1] = psq[i] + input_cache[i] * input_cache[i];
        }
        const float e_floor = 1e-12f + 1e-7f * psq[2 * OLA_LAG + OLA_L];

        // ---- Build decimated copies for the COARSE search (phase-0, by 4).
        //
        // Raw decimation, no anti-alias filter. That was worth checking rather
        // than assuming: adding a [1,2,1]/4 lowpass first changed the output by
        // literally nothing on every test signal. The reason is that the coarse
        // stage only has to land within +-3 of the true peak for the full-rate
        // fine stage to finish the job, and it does. (Widening the fine window to
        // +-6 likewise changed nothing, for the same reason.) The signal is also
        // already band-limited to ~10.6 kHz by the analog Sallen-Key and the
        // 192->48 kHz decimation FIR, so there is little up there to alias.
        //
        // Decimating by 8 keeps each evaluation short enough that the +-1092 sweep
        // costs LESS than the old +-450 sweep at /2 did. Lag step is 8, so the
        // fine window below has to be at least +-4 to close the gap.
        static float target_dec[OLA_L / 8];
        for(size_t k = 0; k < OLA_L / 8; ++k) target_dec[k] = target[8 * k];

        static float cache_dec[OLA_LAG / 4 + OLA_L / 8];
        static float psq_dec[OLA_LAG / 4 + OLA_L / 8 + 1];
        psq_dec[0] = 0.0f;
        for(size_t m = 0; m < OLA_LAG / 4 + OLA_L / 8; ++m) {
            cache_dec[m] = input_cache[8 * m];
            psq_dec[m + 1] = psq_dec[m] + cache_dec[m] * cache_dec[m];
        }
        const float e_floor_dec = 1e-12f + 1e-7f * psq_dec[OLA_LAG / 4 + OLA_L / 8];

        // ---- COARSE: decimated dot products. dec offset p -> lag 8p - OLA_LAG,
        // so stepping p by 1 sweeps the range at a lag step of 8.
        // Score is pure normalized correlation -- see the hysteresis note up top.
        float best_coarse = -FLT_MAX;
        int coarse_offset = -OLA_LAG;
        for(int p = 0; p <= OLA_LAG / 4; ++p) {
            float c;
            arm_dot_prod_f32(cache_dec + p, target_dec, OLA_L / 8, &c);
            const float score = ncc(c, psq_dec[p + OLA_L / 8] - psq_dec[p], e_floor_dec);
            if(score > best_coarse) { best_coarse = score; coarse_offset = 8 * p - OLA_LAG; }
        }

        // ---- FINE: full rate. Coarse scores come off decimated data and are not
        // comparable to full-rate ones, so every candidate in the +-3 window --
        // including the coarse winner itself -- is rescored here at full rate.
        // Window is +-5, one more than the +-4 the lag step of 8 strictly needs.
        float best_score = -FLT_MAX;
        int best_lag = coarse_offset;
        for(int i = coarse_offset - 5; i <= coarse_offset + 5; ++i) {
            if(i < -OLA_LAG || i > OLA_LAG) continue;   // stay inside input_cache
            float c;
            const int idx = i + OLA_LAG;
            arm_dot_prod_f32(input_cache + idx, target, OLA_L, &c);
            const float score = ncc(c, psq[idx + OLA_L] - psq[idx], e_floor);
            if(score > best_score) { best_score = score; best_lag = i; }
        }

        // ---- Confidence gate.
        //
        // best_score is cross/||candidate||, so it still scales with ||target||
        // and cannot be compared against a fixed number. Dividing it out here --
        // once, after the search, since it is identical for every lag and so
        // cannot change the argmax -- makes it a cosine similarity in [-1, 1].
        //
        // Once a note decays into the pickup and 12-bit ADC noise there is no real
        // peak left, and an ungated search picks a fresh arbitrary lag every frame.
        // That is a modulating artifact, the worst kind. lag == 0 is the nominal
        // no-time-warp position: its seam error is fixed rather than wandering.
        if(t_energy < OLA_SILENCE * (float)OLA_L) {
            best_lag = 0;                       // nothing there to lock onto
        } else if(best_score / sqrtf(t_energy) < OLA_MIN_CORR) {
            best_lag = 0;                       // peak is not trustworthy
        }

        // overlap-add the best-matching frame. Writes out_buffer (not in_buffer),
        // so the restrict qualifier on `in` is not violated.
        //
        // Frame start is frame_end - OLA_ANCHOR + lag. At the worst case
        // lag = +OLA_LAG the last sample read is frame_end - 1, i.e. the newest
        // sample this frame is allowed to see. (The old anchor of in_write - OLA_W
        // let positive lags read past in_write into the un-overwritten part of the
        // ring, returning x[n - OLA_IN_SIZE + k] -- audio from ~36 ms earlier,
        // spliced in at up to 0.63 Hann weight on roughly half of all frames.)
        const size_t frame_start =
            (frame_end + best_lag - OLA_ANCHOR + 2 * OLA_IN_SIZE) % OLA_IN_SIZE;
        for(size_t i = 0; i < OLA_W; ++i) {
            size_t dst = (out_write + i) % OLA_OUT_SIZE;
            s->out_buffer[dst] += in[(frame_start + i) % OLA_IN_SIZE] * s->hann[i];
        }
        s->prev_start = frame_start;   // next frame's correlation target anchor
        s->out_write = (out_write + OLA_SS) % OLA_OUT_SIZE;
        s->in_count -= OLA_SA;
    }

    // Decimate by 2 (read every other sample) => octave up.
    for (size_t i = 0; i < buffer_size; ++i) {
        output_buffer[i] = s->out_buffer[s->out_read];
        s->out_buffer[s->out_read] = 0.0f;                          // clear read sample
        s->out_buffer[(s->out_read + 1) % OLA_OUT_SIZE] = 0.0f;     // clear skipped sample
        s->out_read = (s->out_read + 2) % OLA_OUT_SIZE;
    }
}
