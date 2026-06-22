#include "ola.h"
static inline float square(float x) { return x * x; }

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

        // #3: ||tail|| is gone. The search maximizes cross/(||cand||*||tail||),
        //     and ||tail|| is constant across lags, so it cancels from every
        //     comparison. Track the winner as the (cross, ||cand||^2) pair and
        //     defer the sqrt/divide entirely. (No 1e-12 guard needed anymore.)
        float sum_cand_squared = 0.0f;   // running ||cand||^2 of the current lag
        float best_cross   = 0.0f;       // cross-correlation of the best lag
        float best_cand_sq = 0.0f;       // ||cand||^2 of the best lag

        // Seed from the first lag (-OLA_LAG).
        for(size_t j = 0; j < OLA_L; ++j) {
            float c = in[(in_write + j - OLA_W - OLA_LAG + OLA_IN_SIZE) % OLA_IN_SIZE];
            sum_cand_squared += c * c;
            best_cross       += c * tail[j];   // tail from cache, no modulo
        }
        best_cand_sq = sum_cand_squared;
        int max_course_correlation_offset = -OLA_LAG;

        for(int i = -OLA_LAG + 1; i <= OLA_LAG; i += 4) {

            // slide ||cand||^2: drop old leading edge, add new trailing edge
            sum_cand_squared -= square(in[(in_write + i - 1         - OLA_W + OLA_IN_SIZE) % OLA_IN_SIZE]);
            sum_cand_squared += square(in[(in_write + i + OLA_L - 1 - OLA_W + OLA_IN_SIZE) % OLA_IN_SIZE]);

            // full cross-correlation against the cached (contiguous) tail
            float cross = 0.0f;
            for(size_t j = 0; j < OLA_L; ++j) {
                cross += in[(in_write + i + j - OLA_W + OLA_IN_SIZE) % OLA_IN_SIZE] * tail[j];
            }

            // #3: compare cross/sqrt(cand_sq) WITHOUT sqrt or divide.
            //     for positives:  a/sqrt(p) > b/sqrt(q)  <=>  a^2 * q > b^2 * p
            //     a non-positive cross is anti-correlation: never a good seam,
            //     so it can't beat a positive best.
            int better;
            if(cross <= 0.0f)           better = 0;
            else if(best_cross <= 0.0f) better = 1;   // any positive beats non-positive
            else better = (cross * cross * best_cand_sq) >
                          (best_cross * best_cross * sum_cand_squared);

            if(better) {
                best_cross   = cross;
                best_cand_sq = sum_cand_squared;
                max_course_correlation_offset = i;
            }
        }

        //fine search
        int max_fine_correlation_offset = 0;
        for(int i = max_course_correlation_offset - 4; i <= max_course_correlation_offset + 4 ; ++i) {

            sum_cand_squared -= square(in[(in_write + i - 1         - OLA_W + OLA_IN_SIZE) % OLA_IN_SIZE]);
            sum_cand_squared += square(in[(in_write + i + OLA_L - 1 - OLA_W + OLA_IN_SIZE) % OLA_IN_SIZE]);

            float cross = 0.0f;
            for(size_t j = 0; j < OLA_L; ++j) {
                cross += in[(in_write + i + j - OLA_W + OLA_IN_SIZE) % OLA_IN_SIZE] * tail[j];
            }
            int better;
            if(cross <= 0.0f)           better = 0;
            else if(best_cross <= 0.0f) better = 1;   // any positive beats non-positive
            else better = (cross * cross * best_cand_sq) >
                          (best_cross * best_cross * sum_cand_squared);

            if(better) {
                best_cross   = cross;
                best_cand_sq = sum_cand_squared;
                max_fine_correlation_offset = i;
            }
        }
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
