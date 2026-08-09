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
channel = pixel + 4 * column + 32 * (fifo % 4)
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
