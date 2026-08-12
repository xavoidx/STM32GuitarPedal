# GuitarPedal v1.0: Phase Vocoder / WSOLA Pitch-shifter

Pitch-shift guitar pedal on an STM32F446RE development board. 

## Design Overview

The pedal itself contains 2nd order analog Sallen-Key anti-aliasing and reconstruction filters for the ADC and DAC running at 192kHz. The I/O are downsampled/upsampled x4, meaning the processing is done at 48kHz. Interpolation and Decimation are done with the polyphase FIR filter functions from the CMSIS-DSP library. The board contains analog tone and volume pots, as well as a digital dry/wet knob. A standard 2x buffer with half-complete and full-complete DMA callbacks is used for processing.

The repo has two different implementations for the pitch shifter; a Phase Vocoder suited for polyphonic instruments, and a WSOLA approach suited for monophonic instruments. The phase vocoder was entirely written by me, whereas the WSOLA implementation had some major optimizations and improvements made by pairing with Claude. I was able to learn a lot of optimizations tricks/idioms from the LLM that were crucial for meeting the audio callback deadline--which is tight on the F446 (~5.33 ms with 256-sample buffer). The start and end of the callback toggle a GPIO pin on the board that can be monitored with a scope and allows gauging of callback headroom. 

It is living on a breadboard currently, but I am working on a simple PCB for it, as the breadboard is unsurprisingly noisy. 

## Two implementations of an octave-up pitch shifting algorithm

Both approaches use some form of OLA, or Overlap-Add. By setting an analysis hop H_a and synthesis hop H_s, we can take windowed overlapping frames of the input every H_a, re-accumulate with them spaced apart by H_s, and be left with a new waveform (H_s / H_a) times as long. We can then play this new waveform back at a sped-up rate, effectively increasing the pitch. The OLA is necessary because speeding up playback requires the read pointer to eat samples quicker than are being written to the input, and thus new samples must be 'created'. The overlapping frames must be windowed such that there is constant gain. This is determined by how many times a window is applied (1 for WSOLA, 2 for PV), and the ratio H_s / FRAME_SIZE. 

Basic OLA creates phase problems that make the simple form of the algorithm unusable. When the frames are recombined, the overlap-add causes frequencies--which have now been shifted in time--to combine back with themselves at the overlapping seams, but now at different phases. This leads to detrimental cancellations.

To fix this, some form of phase correction algorithm must be applied before recombining the frames:

### 1) Phase Vocoder    (*/Core/Src/phase_vocoder.c*)

By performing an FFT on each OLA frame (STFT), the phase of a frequency bin can be compared to the phase of the previous frame at that bin to determine the true value of a partial. Using this, we can rotate the phase of that bin to the phase that frequency will be at when it is shifted ahead by k * H_s. This way, we get frequency specific phase alignment for OLA reconstruction.

This method has the consequence of smearing transients, as frequencies are shifted independently of one another and no longer exhibit time-domain coherence. There are many solutions of varying complexities to this problem, but the one I implemented is a simple global phase reset that occurs when the energy difference between a current frame's bin and that bin's previous frame (the 'flux' at that bin) exceeds a threshold.

The advantage of the Phase Vocoder over time domain methods is it can handle polyphonic instruments without creating a 'warbly' sound. Most time domain methods are limited to monoponic instruments, but have better transient preservation. The phase vocoder also needs to add significant latency to perform its FFT on a large enough frame to sound good for lower frequencies. I found that FRAME_SIZE = 2048 is about the minimum for decent performance, but this introduces ~42.6 ms of latency. 

CMSIS_DSP is used here for efficient implementation of FFT.

### 2) WSOLA     (*/Core/Src/ola.c*)

WSOLA, or 'Waveform Synchronous Overlap-Add", is a different approach to maintaining phase coherence that stays in the time domain. In the analysis phase, a small lag window is scanned for an analysis candidate that correlates the strongest with the tail of the synthesis buffer that the frame will be overlapped onto. The lag window needs to be wide enough to be able to find a sufficient candidate for low frequencies that have a large wavelength. 

The problem with WSOLA is that the dominant period of the waveform will take precedence in the searching algorithm, thus making it best suited for monophonic instruments, as there will only be phase conherence for one dominant period. Longer lag windows can accomodate for simple dyad harmonic ratios like octaves, thirds, and fifths, as the total wavelength of these ratios are still relatively short, especially at higher pitches. But full chords sound warbly.

The WSOLA implementation in this repo is able to search a fairly large lag window while healthily meeting the callback deadline. It gives very respectable performance on guitar, even during phrases with 2-3 notes playing at the same time.  

## Navigating the Repo

/Notebooks/InterpolationFilter.ipynb -> Decimation and Interpolation filter design

/Core/ -> Source files

(TODO) /Demos/ -> Audio files demoing the algorithm through a guitar amp 

(TODO) /Pics/ -> Pictures of the pedal in its breadboard form 



