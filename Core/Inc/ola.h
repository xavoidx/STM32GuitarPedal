#ifndef OLA_H
#define OLA_H
#include <stddef.h>
#include <string.h>
#include <math.h>
#include "arm_math.h"

#define OLA_W        1024   // frame size (21ms @ 48kHz)
#define OLA_SA       256    // analysis hop  (input advances by this)
#define OLA_SS       512    // synthesis hop (output advances by this, = W/2 for COLA)
#define OLA_IN_SIZE  2048   // input ring buf  (2×W)
#define OLA_OUT_SIZE 4096   // output ring buf (4×W, pre-decimation)
#define OLA_L 		 ( 450 )    // correlation length (tail of output is Ss)
#define OLA_LAG      300    // Search correlation windows +-OLA_LAG from in_write

typedef struct _OLA_State {
    float in_buffer[OLA_IN_SIZE];
    float out_buffer[OLA_OUT_SIZE];
    float hann[OLA_W];
    size_t in_write;
    size_t in_count;
    size_t out_write;
    size_t out_read;
    int    prev_lag;        // last frame's chosen lag, for hysteresis
} OLA_State;

void OLA_Init(OLA_State *instance);
void OLA_Process(OLA_State* ola_state, float* input_buffer, float* output_buffer, size_t buffer_size);

#endif //end OLA_H ifndef
