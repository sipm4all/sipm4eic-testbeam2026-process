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

## timing_pairs.C

Writes one delta-t histogram for every independent pair of DO channels inside each TIMING scintillator. It reads only the triggered-frame `timing` collection and uses the calibrated `hit.time` value stored in the file.

For every frame, the macro keeps the earliest hit in each DO channel. It then fills:

```cpp
deltat = time_A - time_B
```

for all independent pairs `A < B` among the 32 DO channels, separately for:

```text
TIMING0: FIFOs 0..3
TIMING1: FIFOs 4..7
```

The electronics channel is:

```cpp
eoch = pixel + 4 * column
```

and it is converted to detector-output channel with:

```cpp
eo2do[32] = {22,20,18,16,24,26,28,30,25,27,29,31,23,21,19,17,
             9,11,13,15,7,5,3,1,6,4,2,0,8,10,12,14}
```

Histograms are named:

```text
deltat_A_B_timingC
```

where `A` and `B` are two-digit DO channel numbers and `C` is `0` or `1`, for example:

```text
deltat_00_01_timing0
deltat_12_27_timing1
```

Example:

```bash
root -l 'macros/example/timing_pairs.C("triggered.timing.root")'
```

or interactively:

```cpp
.L macros/example/timing_pairs.C
timing_pairs("triggered.timing.root", "timing_pairs.root")
```

Arguments:

```cpp
void timing_pairs(const char *filename,
                  const char *outfilename = "timing_pairs.root",
                  int nbins = 2048,
                  double range = 32.)
```

## timing_event_display.C

Interactive event display for the TIMING scintillators. It reads triggered-frame data with `trigger_reader.h`, applies a `[CHANNEL]` timing-offset file such as the Cherenkov-860 calibration, and draws one frame at a time as two 4 by 8 DO-channel maps plus the corresponding per-scintillator time distributions:

```text
TIMING0: FIFOs 0..3
TIMING1: FIFOs 4..7
```

The macro keeps the earliest hit per DO channel. By default it plots:

```text
calibrated_time - earliest_calibrated_TIMING_time_in_frame
```

so the color scale shows the event time pattern rather than the large absolute clock. The absolute reference time is printed in the terminal.

Run interactively from a directory containing the calibration file:

```cpp
.L macros/example/timing_event_display.C
timing_event_display("triggered.timing.root",
                     "timing_offsets_from_cherenkov860.conf")
```

Press Enter to move to the next event, or type `q` then Enter to stop.

Arguments:

```cpp
void timing_event_display(const char *filename,
                          const char *calibfilename = "timing_offsets_from_cherenkov860.conf",
                          int first_event = 0,
                          int max_events = 0,
                          bool relative_time = true,
                          bool require_trigger = false)
```

## timing0_minuit.C

Fits a diagnostic set of TIMING0 channel offsets using only relative time differences between neighbouring DO channels. The macro is intentionally limited to TIMING0 while we study whether local channel-to-channel constraints are a better calibration handle than global scintillator means.

For every frame, it keeps the earliest hit in each TIMING0 electronics channel:

```cpp
eoch = pixel + 4 * column
```

The electronics channel is converted to DO channel using the same `eo2do` table as `timing_pairs.C`. DO channel coordinates are:

```cpp
x = DO % 4
y = DO / 4
```

The fit uses all channel pairs whose Manhattan distance in this DO grid is at most one:

```cpp
abs(xA - xB) + abs(yA - yB) <= 1
```

For each selected pair, the residual entering the `TMinuit` objective is:

```cpp
dt_corrected =
    (timeA - offsetA)
  - (timeB - offsetB)
```

and the minimized objective is the mean of `dt_corrected^2` over all selected neighbour-pair entries. There are 32 electronics-channel offsets, with `EOCH 0` fixed to zero.

Example:

```bash
root -l 'macros/example/timing0_minuit.C("triggered.timing.root")'
```

or interactively:

```cpp
.L macros/example/timing0_minuit.C
timing0_minuit("triggered.timing.root",
               "timing0_minuit.root",
               "timing0_offsets.conf")
```

Arguments:

```cpp
void timing0_minuit(const char *filename,
                    const char *outfilename = "timing0_minuit.root",
                    const char *calibfilename = "timing0_offsets.conf",
                    int max_calls = 5000,
                    double delta_range = 32.)
```

The ROOT output contains summary histograms `hDeltaBefore`, `hDeltaAfter`, `hOffsetEo`, `hOffsetDo`, `hPairCountDo`, plus before/after histograms for every fitted neighbouring DO pair. The optional text output is a `[CHANNEL]` calibration snippet for TIMING0 FIFOs 0..3.

## timing_calib.C

Only frames containing at least one `trigger` hit are used. Frames without a trigger hit are counted and skipped before timing-channel selection.

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
hFirstPositionDistance
hDeltaVsFirstPositionDistance
hFirstPositionDxDy
hDeltaSameFirstDetectorChannel
hDeltaSameFirstDetectorChannel00 ... hDeltaSameFirstDetectorChannel31
```

`hFirstDetectorChannel0Vs1` uses the `eo2do` detector-channel mapping and is the preferred plot for checking whether the earliest-light position proxy is correlated between TIMING0 and TIMING1. `hFirstElectronicsChannel0Vs1` is kept as a mapping cross-check.

For the position-distance diagnostics, the detector channel is converted to matrix coordinates with `x = detector_channel % 4` and `y = detector_channel / 4`. `hFirstPositionDistance` stores `sqrt((x0 - x1)^2 + (y0 - y1)^2)` between the TIMING0 and TIMING1 first-hit cell centers. `hFirstPositionDxDy` stores the signed displacement `(x0 - x1, y0 - y1)`, and `hDeltaVsFirstPositionDistance` checks whether the timing difference depends on that separation.

`hDeltaSameFirstDetectorChannel` and the 32 per-channel histograms `hDeltaSameFirstDetectorChannel00` through `hDeltaSameFirstDetectorChannel31` use only events where TIMING0 and TIMING1 have the same first-hit detector channel. These plots isolate the residual `timing0_mean - timing1_mean` separately for each matched detector cell.

The spread-sign check is written as:

```text
hTiming0RelativeVsSpread0
hTiming1RelativeVsSpread1
```

Here the event reference is `0.5 * (timing0_mean + timing1_mean)`. The plotted values are `timing0_mean - reference` and `timing1_mean - reference` versus the corresponding scintillator spread. This tests directly whether each scintillator becomes later when its own internal channel-time spread grows.
