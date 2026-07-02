#include "ola.h"
static inline float square(float x) { return x * x; }

// Lag-continuity bias: weight a correlation by how close its lag is to the
// previous frame's chosen lag. Returns 1.0 at prev_lag, falling linearly to
// (1 - LAG_BIAS) at the search-range edge. Damps frame-to-frame lag jitter
// (period vs 2*period jumps) that shows up as warble. LAG_BIAS is the knob:
// larger = stickier (resists jitter more, slower to track real pitch changes);
// set it to 0 to disable.
#define LAG_BIAS 0.2f
static inline float lag_weight(int lag, int prev_lag) {
    return 1.0f - LAG_BIAS * fabsf((float)(lag - prev_lag)) / (float)OLA_LAG;
}

void OLA_Init(OLA_State *s) {
    memset(s, 0, sizeof(*s));               // FIX: sizeof(*s), not sizeof(s) (pointer)
    for(size_t i = 0; i < OLA_W; ++i) {
        // FIX: standard Hann is 0.5 - 0.5cos (zero at edges, 1 at center).
        // The '+' form is an inverted window and destroys the signal.
        s->hann[i] = 0.5f - 0.5f * cosf((2.0f * M_PI * i) / (OLA_W - 1));
    }
    // Prime the read pointer to trail the write frontier by one full frame (W),
    // so every output sample has received all overlapping frames before read.
    s->out_read = OLA_OUT_SIZE - OLA_W;
}

__attribute__((optimize("fast-math")))
void OLA_Process(OLA_State* s, float* input_buffer, float* output_buffer, size_t buffer_size) {

    // Push the new block into the input ring buffer.
    for(size_t i = 0; i < buffer_size; ++i) {
        s->in_buffer[s->in_write] = input_buffer[i];
        s->in_write++;
        if(s->in_write >= OLA_IN_SIZE) s->in_write = 0; 
    }
    s->in_count += buffer_size;  
  
    // Fire one analysis frame for every SA samples accumulated.
    while(s->in_count >= OLA_SA) {

        // #4: hoist the read-only input ring into a restrict-qualified local.
        //     in_buffer is only WRITTEN in the push loop above, never inside
        //     this search, so restrict is valid and lets the compiler assume
        //     it never aliases the out_buffer writes. in_write is fixed for the
        //     whole frame; out_write only changes at the end, so cache both.
        const float * restrict in = s->in_buffer;
        const size_t in_write  = s->in_write;
        const size_t out_write = s->out_write;

        // #4: cache the output tail ONCE. It is identical for every lag, so this
        //     kills the per-lag modulo on out_buffer and reads it contiguously
        //     in the hot loop. static => not on the ISR stack (OLA_Process runs
        //     in the DMA callback; non-reentrant, so a shared scratch is fine).
        static float tail[OLA_L];
        for(size_t j = 0; j < OLA_L; ++j) {
            tail[j] = s->out_buffer[(out_write + j) % OLA_OUT_SIZE];
        }

        // Raw cross-correlation search. ||tail|| is constant across lags and
        // cancels; for sustained notes the candidate energy barely varies, so
        // the normalization is dropped (standard real-time WSOLA). Coarse-to-
        // fine: sweep at step 4 to span a full pitch period cheaply, then
        // refine ±3 around the coarse winner. `cross > best_cross` handles
        // anti-correlation naturally (picks the least-bad seam), so no guards.
        /**
         * Cache input; in write - OLA_LAG - OLA_W -> OLA in write + OLA_LAG - OLA_W + OLA_L
         */
        static float input_cache[2 * OLA_LAG + OLA_L];
        for(int i = 0; i < 2 * OLA_LAG + OLA_L; ++i) {
            input_cache[i] = in[(in_write + i - OLA_W - OLA_LAG + OLA_IN_SIZE) % OLA_IN_SIZE];
        }

        // ---- Build decimated copies for the COARSE search (phase-0, by 2).
        // Raw decimation (no AA filter): the periodicity that sets the peak is
        // below the new Nyquist; the full-rate fine stage fixes the residue.
        // (If a high-gain/distorted tone ever misbehaves, replace the copy with
        //  a 2-tap box average: tail_dec[k] = 0.5f*(tail[2k] + tail[2k+1]).)
        static float tail_dec[OLA_L / 2];
        for(size_t k = 0; k < OLA_L / 2; ++k) tail_dec[k] = tail[2 * k];

        static float cache_dec[OLA_LAG + OLA_L / 2];
        for(size_t m = 0; m < OLA_LAG + OLA_L / 2; ++m) cache_dec[m] = input_cache[2 * m];

        // Hysteresis: bias every correlation toward the previous frame's lag
        // (see lag_weight). Applied in BOTH stages because the big jitter jumps
        // happen in the coarse sweep. best_* now hold biased scores, not raw.
        const int prev_lag = s->prev_lag;

        // ---- COARSE: decimated dot products. dec offset p -> lag 2p - OLA_LAG.
        // Step p by 2 (= lag step 4) so each eval is half length (~2x cheaper).
        float best_coarse;
        arm_dot_prod_f32(cache_dec, tail_dec, OLA_L / 2, &best_coarse);
        best_coarse *= lag_weight(-OLA_LAG, prev_lag);   // seed: lag -OLA_LAG
        int coarse_offset = -OLA_LAG;
        for(int p = 2; p <= OLA_LAG; p += 2) {
            float c;
            arm_dot_prod_f32(cache_dec + p, tail_dec, OLA_L / 2, &c);
            int lag = 2 * p - OLA_LAG;
            c *= lag_weight(lag, prev_lag);
            if(c > best_coarse) { best_coarse = c; coarse_offset = lag; }
        }

        // ---- FINE: full rate, RE-SEEDED. Coarse scores are decimated and not
        // comparable to full-rate, so recompute (and re-bias) the incumbent at
        // full resolution here, then refine +-3 around it.
        float best_cross;
        arm_dot_prod_f32(input_cache + coarse_offset + OLA_LAG, tail, OLA_L, &best_cross);
        best_cross *= lag_weight(coarse_offset, prev_lag);
        int max_fine_correlation_offset = coarse_offset;
        for(int i = coarse_offset - 3; i <= coarse_offset + 3; ++i) {
            if(i < -OLA_LAG || i > OLA_LAG) continue;   // stay inside input_cache
            if(i == coarse_offset) continue;            // already the incumbent
            float cross;
            arm_dot_prod_f32(input_cache + i + OLA_LAG, tail, OLA_L, &cross);
            cross *= lag_weight(i, prev_lag);
            if(cross > best_cross) { best_cross = cross; max_fine_correlation_offset = i; }
        }
        s->prev_lag = max_fine_correlation_offset;   // remember for next frame
        // overlap-add the best-matching frame. Writes out_buffer (not in_buffer),
        // so the restrict qualifier on `in` is not violated.
        for(size_t i = 0; i < OLA_W; ++i) {
            size_t src = (in_write + max_fine_correlation_offset - OLA_W + i + OLA_IN_SIZE) % OLA_IN_SIZE;
            size_t dst = (out_write + i) % OLA_OUT_SIZE;
            s->out_buffer[dst] += in[src] * s->hann[i];
        }
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
