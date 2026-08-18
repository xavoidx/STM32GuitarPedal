#ifndef PHASE_VOCODER_H
#define PHASE_VOCODER_H

#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <complex.h>
#include <string.h>
#include "arm_math.h"

/**
 * NOTE: Ring buffers are implemented with bit masks and intentional
 * counter overflow; size_t indexes wrap at 2^32 and thus can safely
 * overflow with buffer sizes of 2^N. Therefore buffers must be of size
 * 2^N or else counters will break.  
 */
#define PV_FRAME_SIZE (1U << 11) /* (2048) Defines the number of samples of latency introducted*/
#define PV_IN_SIZE (PV_FRAME_SIZE * 2)
#define PV_OUT_SIZE (PV_FRAME_SIZE * 2)

#define PV_TRANSIENT_MULT 6
#define PV_TRANSIENT_FLOOR 50
#define PV_FLUX_HISTORY_SIZE (1U << 4)

#define RFFT 0x00
#define IRFFT 0x01

typedef struct {
    float sine[PV_FRAME_SIZE]; 
    float phase_acc[PV_FRAME_SIZE];
    
    /*Ping pong buffers for spectral operations*/
    float fft_in[PV_FRAME_SIZE]; 
    float fft_out[PV_FRAME_SIZE];
    float fft_prev[PV_FRAME_SIZE];

    float in_buffer[PV_IN_SIZE];
    float out_buffer[PV_OUT_SIZE];
    size_t in_wr;
    size_t in_count;
    size_t out_wr;
    size_t out_rd;
    
    size_t flux_history_wr;
    float flux_history[PV_FLUX_HISTORY_SIZE];
    float prev_magnitude[PV_FRAME_SIZE / 2]; 
} Phase_Vocoder;

void pv_init(Phase_Vocoder* pv);
void pv_process(Phase_Vocoder* pv, const float* in, float* out, size_t buffer_size);

#endif //end ifndef PHASE_VOCODER_H
