#include "phase_vocoder.h"
#pragma GCC diagnostic warning "-Wdouble-promotion"

arm_rfft_fast_instance_f32 rfft;
static const float pi_f = (float) M_PI;
static const float e_floor = 1e-12;


void pv_init(Phase_Vocoder* pv) {
    memset(pv, 0, sizeof(*pv));
    for(size_t i = 0; i < PV_FRAME_SIZE; ++i) {
        pv->sine[i] = 0.5f - 0.5f * cosf((2 * pi_f * i) / PV_FRAME_SIZE);
        if(i < PV_FRAME_SIZE / 2) {
            pv->phase_acc[2 * i] = 1.0f;
            pv->phase_acc[2 * i + 1] = 0.f;
            pv->fft_prev[2 * i] = 1.0f;
            pv->fft_prev[2 * i + 1] = 0.f;
        }
    }
    arm_rfft_fast_init_f32(&rfft, PV_FRAME_SIZE);
}

static int fcompare(const void* arg1, const void* arg2) {
    if(*(float*)arg1 > *(float*)arg2) return 1;
    if(*(float*)arg1 < *(float*)arg2) return -1;
    return 0;
}

static inline float pv_get_median(Phase_Vocoder* pv) {
    float flux_scratch[PV_FLUX_HISTORY_SIZE];
    memcpy(flux_scratch, pv->flux_history, PV_FLUX_HISTORY_SIZE);
    qsort(flux_scratch, PV_FLUX_HISTORY_SIZE, sizeof(*(flux_scratch)), fcompare);
    return flux_scratch[PV_FLUX_HISTORY_SIZE / 2]; 
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

        float flux = 0.f; 
        for(size_t i = 1; i < PV_FRAME_SIZE / 2; i++) {
            float bin_magnitude = sqrtf(pv->fft_out[2 * i] * pv->fft_out[2 * i] + pv->fft_out[2 * i + 1] * pv->fft_out[2 * i + 1]);
            if(i > 4 /*Guard against energy spikes near DC*/) {
                flux += fmaxf(bin_magnitude - pv->prev_magnitude[i], 0.f);
            }
            pv->prev_magnitude[i] = bin_magnitude;   
        }
        /**
         * Alias the prev_magnitude array for readability; 
         * prev_magnitude holds the current magnitude at this stage.
         */
        float* curr_magnitude = pv->prev_magnitude; 

        if(flux > fmaxf(PV_TRANSIENT_MULT * pv_get_median(pv), PV_TRANSIENT_FLOOR)) {    
            /*Short-circuit input to output leaving phase unchanged*/
            memcpy(pv->fft_in, pv->fft_out, PV_FRAME_SIZE * sizeof(float));

            /*Reset phase accumulation to current phase*/
            for(size_t i = 1; i < PV_FRAME_SIZE / 2; i++) {
                float mag_divisor = curr_magnitude[i] + e_floor;
                pv->phase_acc[2 * i] = pv->fft_out[2 * i] * mag_divisor;
                pv->phase_acc[2 * i + 1] = pv->fft_out[2 * i + 1] * mag_divisor;

                pv->fft_prev[2 * i] = pv->fft_out[2 * i];
                pv->fft_prev[2 * i + 1] = pv->fft_out[2 * i + 1];
            }
        }
        else {
            for(size_t i = 1; i < PV_FRAME_SIZE / 2; i++) {    
                /*Complex floats*/
                float phase_diff[2] = {0.f};
                float phase_sqr[2] = {0.f};
                float prev_conj[2] = {
                    pv->fft_prev[2 * i],
                    -1 * pv->fft_prev[2 * i + 1]
                };

                /*Multiply current bin by conjugate of previous*/
                phase_diff[0] = pv->fft_out[2 * i] * prev_conj[0] - pv->fft_out[2 * i + 1] * prev_conj[1];
                phase_diff[1] = pv->fft_out[2 * i] * prev_conj[1] + pv->fft_out[2 * i + 1] * prev_conj[0];
    
                pv->fft_prev[2 * i] = pv->fft_out[2 * i];
                pv->fft_prev[2 * i + 1] = pv->fft_out[2 * i + 1];
    
                /*Square that result as alpha = 2, multiplying phase by 2 is same as squaring phasor*/
                phase_sqr[0] = phase_diff[0] * phase_diff[0] - phase_diff[1] * phase_diff[1];
                phase_sqr[1] = 2 * phase_diff[0] * phase_diff[1];
    
                /*Accumulate total phase*/
                float phase_acc_curr[2] = {0.f};
                phase_acc_curr[0] = pv->phase_acc[2 * i];
                phase_acc_curr[1] = pv->phase_acc[2 * i + 1];
                pv->phase_acc[2 * i] = phase_sqr[0] * phase_acc_curr[0] - phase_sqr[1] * phase_acc_curr[1];
                pv->phase_acc[2 * i + 1] = phase_sqr[0] * phase_acc_curr[1] + phase_sqr[1] * phase_acc_curr[0];
    
                /*Normalize accumulator*/
                float acc_magnitude = sqrtf(pv->phase_acc[2 * i] * pv->phase_acc[2 * i] + pv->phase_acc[2 * i + 1] * pv->phase_acc[2 * i + 1]);
                float mag_divisor = 1.0f / (acc_magnitude + e_floor); /*Pre-calculate divisor so there's only one UDIV*/
                pv->phase_acc[2 * i] *= mag_divisor;
                pv->phase_acc[2 * i + 1] *= mag_divisor; 

                /*Shift bin output by this correct phase*/
                pv->fft_in[2 * i] = curr_magnitude[i] * pv->phase_acc[2 * i];
                pv->fft_in[2 * i + 1] = curr_magnitude[i] * pv->phase_acc[2 * i + 1];
            }
        }        
        pv->fft_in[0] = pv->fft_out[0];
        pv->fft_in[1] = pv->fft_out[1];
        arm_rfft_fast_f32(&rfft, pv->fft_in, pv->fft_out, IRFFT);
        pv->flux_history[pv->flux_history_wr++ & (PV_FLUX_HISTORY_SIZE - 1)] = flux;
        
        /*Overlap-add @ intervals of H_s*/
        for(size_t i = 0; i < PV_FRAME_SIZE; ++i) {
            pv->out_buffer[(pv->out_wr + i) & (PV_OUT_SIZE - 1)] += pv->sine[i] * pv->fft_out[i];
        }
        pv->out_wr += H_s;

        pv->in_count -= H_a;
    } 
    /*Read samples at 2x rate*/
    for(size_t i = 0; i < buffer_size; ++i) {
        out[i] = pv->out_buffer[(pv->out_rd - 2 * H_s) & (PV_OUT_SIZE - 1)];
        pv->out_buffer[(pv->out_rd - 2 * H_s) & (PV_OUT_SIZE - 1)] = 0.f;
        pv->out_buffer[(pv->out_rd + 1 - 2 * H_s) & (PV_OUT_SIZE - 1)] = 0.f;
        pv->out_rd += 2;
    }
}