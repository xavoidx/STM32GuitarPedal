#include "phase_vocoder.h"
#pragma GCC diagnostic warning "-Wdouble-promotion"

arm_rfft_fast_instance_f32 rfft;
static const float pi_f = (float) M_PI;
static const float phase_exp[8] = {
    0.f, 
    pi_f / 4,
    pi_f / 2,
    3 * pi_f / 4,
    pi_f, 
    5 * pi_f / 4,
    3 * pi_f / 2,
    7 * pi_f / 4
};
static const float e_floor = 1e-12;
float princarg(float in) {
    int div = (int) (in / (2 * pi_f) + 10.5f); 
    div -= 10;
    return in - (2 * pi_f * div);
}
void pv_init(Phase_Vocoder* pv) {
    memset(pv, 0, sizeof(*pv));
    for(size_t i = 0; i < PV_FRAME_SIZE; ++i) {
        pv->sine[i] = 0.5f - 0.5f * cosf((2 * pi_f * i) / PV_FRAME_SIZE);
        if(i < PV_FRAME_SIZE / 2) {
            pv->phase_acc[2 * i] = 1.0f;
            pv->phase_acc[2 * i + 1] = 0.f;
            pv->phase_prev[2 * i] = 1.0f;
            pv->phase_prev[2 * i + 1] = 0.f;
        }
    }
    arm_rfft_fast_init_f32(&rfft, PV_FRAME_SIZE);
    pv->is_first = 1;

}

__attribute__((optimize("fast-math")))
void pv_process(Phase_Vocoder* pv, const float* in, float* out, size_t buffer_size) {

    const size_t H_a = buffer_size; // Analysis hop 512
    const size_t H_s = H_a * 2; // Synthesis hop 1024 for octave up
    
    for(size_t i = 0; i < buffer_size; ++i) {
        pv->in_buffer[pv->in_wr++] = in[i];
        if(pv->in_wr == PV_IN_SIZE) pv->in_wr = 0; 
    }
    pv->in_count += buffer_size;
    if(pv->in_count >= H_a) {
        for(size_t i = 0; i < PV_FRAME_SIZE; ++i) {
            pv->fft_in[i] = pv->sine[i]
             * pv->in_buffer[(pv->in_wr - PV_FRAME_SIZE + i) & (PV_IN_SIZE - 1)];
        }
        arm_rfft_fast_f32(&rfft, pv->fft_in, pv->fft_out, RFFT); 

        for(size_t i = 1; i < PV_FRAME_SIZE / 2; i++) {
            float phase_diff[2] = {0.f};
            float phase_sqr[2] = {0.f};
            float phase_acc_curr[2] = {0.f};
            float prev_conj[2] = {
                pv->phase_prev[2 * i],
                -1 * pv->phase_prev[2 * i + 1]
            };
            if(pv->is_first) {
                pv->phase_acc[2 * i] = pv->fft_out[2 * i];
                pv->phase_acc[2 * i + 1] = pv->fft_out[2 * i + 1];
            }
            /*Multiply current bin by conjugate of previous*/
            phase_diff[0] = pv->fft_out[2 * i] * prev_conj[0] - pv->fft_out[2 * i + 1] * prev_conj[1];
            phase_diff[1] = pv->fft_out[2 * i] * prev_conj[1] + pv->fft_out[2 * i + 1] * prev_conj[0];

            pv->phase_prev[2 * i] = pv->fft_out[2 * i];
            pv->phase_prev[2 * i + 1] = pv->fft_out[2 * i + 1];

            /*Square that result as alpha = 2, multiplying phase by 2 is same as squaring phasor*/
            phase_sqr[0] = phase_diff[0] * phase_diff[0] - phase_diff[1] * phase_diff[1];
            phase_sqr[1] = 2 * phase_diff[0] * phase_diff[1];

            /*Accumulate total phase*/
            phase_acc_curr[0] = pv->phase_acc[2 * i];
            phase_acc_curr[1] = pv->phase_acc[2 * i + 1];
            pv->phase_acc[2 * i] = phase_sqr[0] * phase_acc_curr[0] - phase_sqr[1] * phase_acc_curr[1];
            pv->phase_acc[2 * i + 1] = phase_sqr[0] * phase_acc_curr[1] + phase_sqr[1] * phase_acc_curr[0];

            /*Normalize accumulator*/
            float acc_magnitude = sqrtf(pv->phase_acc[2 * i] * pv->phase_acc[2 * i] + pv->phase_acc[2 * i + 1] * pv->phase_acc[2 * i + 1]);
            float mag_divisor = 1.0f / (acc_magnitude + e_floor);
            pv->phase_acc[2 * i] *= mag_divisor;
            pv->phase_acc[2 * i + 1] *= mag_divisor; 
            /*Shift bin output by this correct phase*/
            float bin_magnitude = sqrtf(pv->fft_out[2 * i] * pv->fft_out[2 * i] + pv->fft_out[2 * i + 1] * pv->fft_out[2 * i + 1]);
            pv->fft_in[2 * i] = bin_magnitude * pv->phase_acc[2 * i];
            pv->fft_in[2 * i + 1] = bin_magnitude * pv->phase_acc[2 * i + 1];
        }
        pv->fft_in[0] = pv->fft_out[0];
        pv->fft_in[1] = pv->fft_out[1];
        arm_rfft_fast_f32(&rfft, pv->fft_in, pv->fft_out, IRFFT);
    
        for(size_t i = 0; i < PV_FRAME_SIZE; ++i) {
            pv->out_buffer[(pv->out_wr + i) & (PV_OUT_SIZE - 1)] += pv->sine[i] * pv->fft_out[i];
        }
        pv->out_wr += H_s;
        pv->in_count -= H_a;
        pv->is_first = 0;

    } 
    for(size_t i = 0; i < buffer_size; ++i) {
        out[i] = pv->out_buffer[(pv->out_rd - 2 * H_s) & (PV_OUT_SIZE - 1)];
        pv->out_buffer[(pv->out_rd - 2 * H_s) & (PV_OUT_SIZE - 1)] = 0.f;
        pv->out_buffer[(pv->out_rd + 1 - 2 * H_s) & (PV_OUT_SIZE - 1)] = 0.f;
        pv->out_rd += 2;
    }
}