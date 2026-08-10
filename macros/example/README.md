# Example ROOT Macros

This directory contains small ROOT macros showing how to inspect processed data.

## trigger_reader.C

Demonstrates the header-only triggered-frame reader API:

```bash
root -l 'macros/example/trigger_reader.C("frames.root")'
```

It iterates over spills and frames and prints the number of trigger, timing, and Cherenkov hits in each frame.

## deltat.C

Reads the triggered-frame `frames` tree and fills a 2D delta-t histogram using all hits contained in each stored frame. The ROOT function name is `deltat`, matching the filename. Example:

```bash
root -l 'macros/example/deltat.C("frames.root", 9, 200, 32, -1, -1)'
```

## timing.C

Analyzes timing-scintillator hits from triggered frame files. It uses `macros/lib/trigger_reader.h` and reads only the `timing` hit collection, which is produced by `trigger.cc` for hits with:

```text
type == 1 && device == 200
```

The timing system is made of two scintillators:

```text
TIMING0: device kc705-200, FIFOs 0..3
TIMING1: device kc705-200, FIFOs 4..7
```

Each scintillator is read out by 32 SiPM channels. The macro uses this local channel index:

```cpp
channel = pixel + 4 * column
```

For each triggered frame, the macro keeps the earliest hit per timing channel, then estimates a robust mean time separately for TIMING0 and TIMING1. The robust mean is computed by:

1. calculating the median channel time;
2. rejecting channels farther than `outlier_window` from the median;
3. averaging the remaining channels;
4. requiring at least `min_channels` surviving channels.

This rejects isolated channel outliers while still requiring a reasonable number of participating SiPM channels.

The macro writes:

```text
hDeltaTiming0   channel residuals: hit.time - timing0_mean vs channel
hDeltaTiming1   channel residuals: hit.time - timing1_mean vs channel
hDelta          timing0_mean - timing1_mean
```

Example:

```bash
root -l 'macros/example/timing.C("triggered.root")'
```

or interactively:

```cpp
.L macros/example/timing.C
timing("triggered.root", "timing.root", 16, 2.0)
```

Arguments:

```cpp
void timing(const char *filename,
            const char *outfilename = "timing.root",
            int min_channels = 16,
            double outlier_window = 2.0,
            double residual_range = 10.0,
            double delta_range = 20.0)
```

The macro uses the calibrated `hit.time` value persisted in the triggered frame file. It does not recompute nominal time from `coarse`, `rollover`, and `fine`.

## timing_calib.C

Builds timing-channel calibration offsets from triggered timing data. It reads the `timing` hit collection through `trigger_reader.h` and fits 64 channel offsets:

```text
TIMING0: 32 channels from FIFOs 0..3
TIMING1: 32 channels from FIFOs 4..7
```

The local channel number is:

```cpp
channel = pixel + 4 * column
```

The fitted time model follows the calibration convention used by `calibrator.cc`:

```cpp
calibrated_time = hit.time - channel_offset
```

For every accepted frame the macro computes robust means for TIMING0 and TIMING1. It first performs damped leave-one-out iterative residual-to-scintillator-mean corrections, keeping the reference channel fixed at zero after every iteration. For this sample, too many undamped iterations overcorrect; the default is three leave-one-out iterations with damping 0.5. Those offsets are then used as the starting point for a full 63-parameter `TMinuit` minimization. The objective contains both the inter-scintillator term and intra-scintillator residual terms, so the minimizer cannot improve `timing0_mean - timing1_mean` by making individual channel residuals worse:

```cpp
sum(delta_weight * (timing0_mean - timing1_mean)^2
    + residual_weight * channel_residual^2)
```

One offset is arbitrary, so the reference channel is fixed to zero:

```text
fifo=0 column=0 pixel=0 offset=0
```

The macro writes a ROOT diagnostic file and a text calibration snippet containing a `[CHANNEL]` section that can be included or concatenated into a `calibrator.cc` calibration file. Since `calibrator.cc` matches exact `(device,fifo,column,pixel)` addresses, the same fitted scintillator-channel offset is written for every FIFO belonging to that scintillator group.

Example:

```bash
root -l 'macros/example/timing_calib.C("triggered.timing.root")'
```

or interactively:

```cpp
.L macros/example/timing_calib.C
timing_calib("triggered.timing.root", "timing_calib.root", "timing_channel_offsets.conf")
```

Arguments:

```cpp
void timing_calib(const char *filename,
                  const char *outfilename = "timing_calib.root",
                  const char *calibfilename = "timing_channel_offsets.conf",
                  int min_channels = 16,
                  double outlier_window = 2.0,
                  double delta_range = 20.0,
                  double offset_range = 20.0,
                  int pre_iterations = 3,
                  double minimizer_step = 0.01,
                  int minimizer_calls = 5000,
                  double delta_weight = 1.0,
                  double residual_weight = 1.0,
                  int max_frames = 0)
```

The ROOT output includes `hDeltaBefore`, `hDeltaAfter`, `hOffset`, `hOffsetValue`, and before/after channel-residual maps for TIMING0 and TIMING1.


Timing calibration diagnostics now also include event-by-event spread histograms:

```text
hTiming0SpreadBefore / hTiming0SpreadAfter
hTiming1SpreadBefore / hTiming1SpreadAfter
hExpectedDeltaFromSpreadBefore / hExpectedDeltaFromSpreadAfter
```

For each accepted frame, the timing spread is the RMS of selected channel times inside one scintillator. The expected contribution to `RMS(timing0_mean - timing1_mean)` from independent channel jitter is estimated as:

```cpp
sqrt(spread0^2 / n0 + spread1^2 / n1)
```

Comparing this estimate with the observed `hDelta` RMS is a useful check of whether the two timing scintillators are limited by independent channel jitter or by an event-by-event effect common to many channels.

The timing diagnostic also writes first-order event-shape correlation histograms:

```text
hDeltaShapeCorrected
hDeltaVsSpread0 / hDeltaVsSpread1
hDeltaVsSlopeX0 / hDeltaVsSlopeX1 / hDeltaVsSlopeY0 / hDeltaVsSlopeY1
hDeltaVsLeftRight0 / hDeltaVsLeftRight1
```

`hDeltaShapeCorrected` subtracts a linear diagnostic model of `timing0_mean - timing1_mean` using per-event shape variables measured inside the two scintillators: channel spread, channel-time slope, left/right asymmetry, and even/odd asymmetry. This correction is diagnostic only; it is not written into the `[CHANNEL]` calibration snippet, because it represents event topology or light-propagation information rather than fixed channel offsets.

The event-shape diagnostics use the detector-channel geometry, not raw electronics order. The electronics channel is first mapped through:

```cpp
eo2do[32] = {22,20,18,16,24,26,28,30,25,27,29,31,23,21,19,17,
             9,11,13,15,7,5,3,1,6,4,2,0,8,10,12,14}
```

The detector channel is interpreted as a 4x8 matrix with `x = detector_channel % 4` and `y = detector_channel / 4`. The diagnostic variables include spread, x/y time gradients, and left/right asymmetry in this detector coordinate system.

The diagnostic ROOT output also includes first-hit timing comparisons:

```text
hDeltaFirstBefore
hDeltaFirstAfter
```

For these histograms, each scintillator time is the earliest selected channel time in that frame, rather than the robust average of selected channels. The same outlier-selected channel set is used, so this compares the timing estimator itself rather than changing event selection.

First-hit detector position diagnostics are also written:

```text
hDeltaVsFirstX0 / hDeltaVsFirstX1
hDeltaVsFirstY0 / hDeltaVsFirstY1
hFirstPosition0 / hFirstPosition1
hDeltaShapePositionCorrected
```

The first hit is interpreted as a beam-position proxy inside the 4x8 detector matrix. `hDeltaShapePositionCorrected` extends the diagnostic shape correction with the first-hit detector `(x,y)` coordinates for both scintillators.

First-hit channel correlations are stored as:

```text
hFirstDetectorChannel0Vs1
hFirstElectronicsChannel0Vs1
```

`hFirstDetectorChannel0Vs1` uses the `eo2do` detector-channel mapping and is the preferred plot for checking whether the earliest-light position proxy is correlated between TIMING0 and TIMING1. `hFirstElectronicsChannel0Vs1` is kept as a mapping cross-check.
