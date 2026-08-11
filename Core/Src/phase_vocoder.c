#include "phase_vocoder.h"

arm_rfft_fast_instance_f32 rfft;

float princarg(float in) {
    float mod = fmodf(in + M_PI, 2 * M_PI);
    if(mod < 0) mod += 2 * M_PI; 
    return mod - M_PI;
}

void pv_init(Phase_Vocoder* pv) {
    memset(pv, 0, sizeof(*pv));
    for(size_t i = 0; i < PV_FRAME_SIZE; ++i) {
        pv->hann[i] = 0.5f - 0.5f * cosf((2 * M_PI * i) / PV_FRAME_SIZE);
    }
    arm_rfft_fast_init_f32(&rfft, PV_FRAME_SIZE);
}

__attribute__((optimize("fast-math")))
void pv_process(Phase_Vocoder* pv, const float* in, float* out, size_t buffer_size) {

    const size_t H_a = buffer_size; // Analysis hop 256
    const size_t H_s = buffer_size * 2; // Synthesis hop 512 for octave up
    
    for(size_t i = 0; i < buffer_size; ++i) {
        pv->in_buffer[pv->in_wr++] = in[i];
        if(pv->in_wr == PV_IN_SIZE) pv->in_wr = 0; 
    }
    pv->in_count += buffer_size;
    if(pv->in_count >= H_a) {
        for(size_t i = 0; i < PV_FRAME_SIZE; ++i) {
            pv->fft_in[i] = pv->hann[i]
             * pv->in_buffer[(pv->in_wr - PV_FRAME_SIZE + i) & (PV_IN_SIZE - 1)];
        }
        arm_rfft_fast_f32(&rfft, pv->fft_in, pv->fft_out, RFFT);
        for(size_t i = 1; i < PV_FRAME_SIZE / 2; i++) {
            float phase, phase_delta, phase_true;
            arm_atan2_f32(pv->fft_out[2 * i + 1], pv->fft_out[2 * i], &phase);
            phase_delta = phase - pv->phase_prev[i];
            pv->phase_prev[i] = phase;
            
            phase_true = princarg(phase_delta - (H_a * i * 2 * M_PI / PV_FRAME_SIZE));
            pv->phase_acc[i] += (H_s * 1.0 / H_a) * (phase_true + H_a * i * 2 * M_PI / PV_FRAME_SIZE);
            pv->phase_acc[i] = princarg(pv->phase_acc[i]);
            
            float magnitude;
            arm_cmplx_mag_f32(pv->fft_out + 2 * i, &magnitude, 1);
            pv->fft_in[0] = pv->fft_out[0];
            pv->fft_in[1] = pv->fft_out[1];
            pv->fft_in[2 * i] = magnitude * arm_cos_f32(pv->phase_acc[i]);
            pv->fft_in[2 * i + 1] = magnitude * arm_sin_f32(pv->phase_acc[i]);
        }
        arm_rfft_fast_f32(&rfft, pv->fft_in, pv->fft_out, IRFFT);

        for(size_t i = 0; i < PV_FRAME_SIZE; ++i) {
            pv->out_buffer[(pv->out_wr + i) & (PV_OUT_SIZE - 1)] += pv->hann[i] * pv->fft_out[i];
        }
        pv->in_count -= H_a;
        pv->out_wr += H_s;

    }
    for(size_t i = 0; i < buffer_size; ++i) {
        out[i] = pv->out_buffer[(pv->out_rd - H_s) & (PV_OUT_SIZE - 1)];
        pv->out_buffer[(pv->out_rd - H_s) & (PV_OUT_SIZE - 1)] = 0.f;
        pv->out_buffer[(pv->out_rd + 1 - H_s) & (PV_OUT_SIZE - 1)] = 0.f;
        pv->out_rd += 2;
    }
}