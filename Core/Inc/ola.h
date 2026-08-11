#ifndef OLA_H
#define OLA_H
#include <stddef.h>
#include <string.h>
#include <math.h>
#include "arm_math.h"

#define OLA_W        1024   // frame size (21ms @ 48kHz)
#define OLA_SA       256    // analysis hop  (input advances by this)
#define OLA_SS       512    // synthesis hop (output advances by this, = W/2 for COLA)
#define OLA_IN_SIZE  4096   // input ring buf; must hold OLA_W + 2*OLA_LAG
#define OLA_OUT_SIZE 4096   // output ring buf (4×W, pre-decimation)
// Correlation length. Must cover at least one period of the lowest note or the
// search cannot tell one period from another: low E (82.4 Hz) is 582.5 samples
// at 48 kHz, so 640 gives 1.10 periods.
//
// This used to be capped at OLA_SS = 512 (0.88 periods) because the correlation
// target was read out of out_buffer, which only holds OLA_SS valid samples past
// out_write. Targeting the raw input instead removes that cap; the new ceiling
// is the causality bound asserted below.
#define OLA_L 		 ( 640 )
// Search window is +-OLA_LAG around the nominal frame start. The requirement is
// that it span the full repeat period of whatever is being played, because the
// splice can only be clean where the waveform actually repeats.
//
//   single note   -> the pitch period. Low E (82.4 Hz) = 582 samples, needs +-291.
//   two notes p:q -> the COMPOSITE period, 1/gcd(f1,f2), which is much longer.
//                    octave (1:2) 110+220 -> 436 samples, needs +-218.
//                    fifth  (2:3) 110+165 -> 873 samples, needs +-437.
//                    fourth (3:4)         -> 1309 samples, needs +-655.
//                    third  (4:5)         -> 1745 samples, needs +-873.
//
//                    min3rd (5:6)         -> 2182 samples, needs +-1091.
//                    min6th (5:8)         -> 2182 samples, needs +-1091.
//
// 1092 covers every common dyad up to and including minor thirds and sixths.
// Measured frame-rate sideband energy, +-450 -> +-1092:
//
//   octave  -58.3 -> -60.5 dB     fourth  + 0.3 -> -19.7 dB
//   fifth   -56.5 -> -47.7 dB     maj6th  - 2.0 -> -19.1 dB
//   maj3rd  - 2.2 -> -35.8 dB     min3rd  - 5.9 -> -27.3 dB
//
// Fourths and sixths plateau near -20 dB rather than going fully transparent:
// their composite period (1309) is more than twice OLA_L, so the correlation
// window cannot see a whole one and partially locks to a single note instead.
// Audible but vastly better than the total breakdown at +-450. The fifth gives
// up ~9 dB (a wider window admits more distant, less well matched candidates)
// and stays far below audibility.
//
// The real price is latency: OLA_ANCHOR grows with OLA_LAG, so this costs ~12 ms
// on top of everything else. Widening further (major seconds need +-1745) is not
// worth it. Revert to 450 for a low-latency, power-chords-only build.
#define OLA_LAG      1092

// Nominal (lag == 0) frame start is in_write - OLA_W - OLA_LAG, so the newest
// sample any candidate frame can touch is in_write - 1 even at lag = +OLA_LAG.
// Buys causality, and costs OLA_LAG samples of latency -- 22.8 ms at OLA_LAG
// = 1092, which is the dominant term in the whole signal chain. Scales 1:1
// with OLA_LAG, so this is what you pay back by narrowing the search.
#define OLA_ANCHOR   (OLA_W + OLA_LAG)

// One-pole DC blocker on the OLA input. fc = fs*(1-R)/(2*pi) ~= 7.6 Hz @ 48 kHz.
// The ADC path has no DC removal beyond the fixed "- 2048" in process_block, and
// V_mid tolerance + op-amp offset leave a residual DC term that dominates the
// cross-correlation once a note decays. See notes in OLA_Process.
#define OLA_DC_R     0.999f

// Confidence gate. The winning score is divided by ||target|| to turn it into a
// true correlation coefficient in [-1, 1]; below OLA_MIN_CORR the peak carries no
// real information (a note decayed into the pickup/ADC noise floor) and the search
// would just pick a different arbitrary lag every frame. Falling back to lag == 0
// trades a random, modulating artifact for a fixed one, which is far less audible.
// OLA_SILENCE is the same idea in absolute terms, as mean square per sample: one
// 12-bit LSB of dither is ~2e-8, so this sits a bit above the quantization floor.
// Measured on synthetic plucks: sustained notes score 0.88-0.99 all the way down
// to -60 dBFS, pure quantization noise scores 0.00, and a hard pluck attack dips
// to ~0.34 for a single frame (the target predates the attack, so nothing matches
// -- gating that frame to lag 0 is the right answer anyway).
//
// 0.3 sits well clear of the noise end rather than splitting the difference,
// because the two failure modes are not symmetric: gating too late just leaves
// the old behaviour, while gating too eagerly forces lag 0 on real notes and
// reintroduces seam buzz. Real playing (vibrato, inharmonicity, pick noise) will
// score lower than this synthetic material, so this is the number most worth
// re-tuning by ear -- raise it and the tail goes dry sooner, lower it and the
// effect stays alive longer into the decay.
#define OLA_MIN_CORR 0.3f
#define OLA_SILENCE  1e-7f   // ~ -70 dBFS RMS; a true no-input guard, not a gate

_Static_assert(OLA_W + 2 * OLA_LAG <= OLA_IN_SIZE,
               "input ring must hold the deepest search read");

// out_read trails out_write by OLA_SS. A sample is final once the frame starting
// on it has been overlap-added, and frames are added before the read loop runs,
// so trailing by one synthesis hop is safe with one frame of margin left over.
// (This was OLA_W; the extra 512 samples were pure latency. It only became safe
// once the correlation target moved off out_buffer -- the old target read exactly
// the region this now zeroes.)
//
// This assumes the caller passes exactly OLA_SA samples per call, which is what
// process_block does. out_write then advances OLA_SS per call and out_read
// 2*buffer_size, so the two stay locked with one frame of margin. Larger blocks
// still average out but let the gap breathe, and the margin is only one frame now.
_Static_assert(OLA_SS >= OLA_SA, "read pointer would overrun the write frontier");

// The correlation target is the previous frame's natural continuation,
// in[prev_start + OLA_SS .. + OLA_L). Worst case prev_start sits at
// frame_end_prev - OLA_ANCHOR + OLA_LAG and frame_end advances by OLA_SA, so
// staying causal needs  OLA_L <= OLA_SA + OLA_ANCHOR - OLA_SS - OLA_LAG  (768).
_Static_assert(OLA_L <= OLA_SA + OLA_ANCHOR - OLA_SS - OLA_LAG,
               "correlation target would read past the newest input sample");
_Static_assert(OLA_L % 8 == 0, "coarse stage decimates OLA_L by 8");
_Static_assert(OLA_LAG % 4 == 0, "coarse stage indexes lag/8 with lag step 8");

typedef struct _OLA_State {
    float in_buffer[OLA_IN_SIZE];
    float out_buffer[OLA_OUT_SIZE];
    float hann[OLA_W];
    size_t in_write;
    size_t in_count;
    size_t out_write;
    size_t out_read;
    size_t prev_start;      // in_buffer index of the previous frame's first sample
    float  dc_x1, dc_y1;    // DC blocker state
} OLA_State;

void OLA_Init(OLA_State *instance);
void OLA_Process(OLA_State* ola_state, const float* input_buffer, float* output_buffer, size_t buffer_size);

#endif //end OLA_H ifndef
