#include "ola.h"

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

void OLA_Process(OLA_State* s, float* input_buffer, float* output_buffer, size_t buffer_size) {

    // Push the new block into the input ring buffer.
    for(size_t i = 0; i < buffer_size; ++i) {
        s->in_buffer[s->in_write] = input_buffer[i];
        s->in_write++;
        if(s->in_write >= OLA_IN_SIZE) s->in_write = 0;   // FIX: >= not > (off-by-one)
    }

    s->in_count += buffer_size;             // FIX: accumulate across calls, don't reset

    // Fire one analysis frame for every SA samples accumulated.
    while(s->in_count >= OLA_SA) {
        for(size_t i = 0; i < OLA_W; ++i) {
            size_t src = (s->in_write - OLA_W + i + OLA_IN_SIZE) % OLA_IN_SIZE;
            size_t dst = (s->out_write + i) % OLA_OUT_SIZE;
            s->out_buffer[dst] += s->in_buffer[src] * s->hann[i];
        }
        s->out_write = (s->out_write + OLA_SS) % OLA_OUT_SIZE;
        s->in_count -= OLA_SA;
    }

    // FIX: emit buffer_size samples every call, OUTSIDE the frame loop.
    // Decimate by 2 (read every other sample) => octave up.
    for (size_t i = 0; i < buffer_size; ++i) {
        output_buffer[i] = s->out_buffer[s->out_read];
        s->out_buffer[s->out_read] = 0.0f;                          // clear read sample
        s->out_buffer[(s->out_read + 1) % OLA_OUT_SIZE] = 0.0f;     // clear skipped sample
        s->out_read = (s->out_read + 2) % OLA_OUT_SIZE;
    }
}
