# Example ROOT Macros

This directory contains small ROOT macros showing how to inspect processed data.

## trigger_reader.C

Demonstrates the header-only triggered-frame reader API:

```bash
root -l 'macros/example/trigger_reader.C("frames.root")'
```

It iterates over spills and frames and prints the number of trigger, timing, and Cherenkov hits in each frame.

## display.C

Interactive Cherenkov event display for triggered frames. It draws one frame at a time using the stored hit `x/y` positions and one square per channel, rather than a regular `TH2D` bin grid. This is important because the channel centres are not assumed to lie on a perfectly regular ROOT histogram grid.

Run:

```bash
root -l 'macros/example/display.C("triggered.root")'
```

It displays only Cherenkov hits. The drawn square size defaults to 3.2 mm and can be changed with the optional `pixel_size` argument. The display axes are fixed to:

```text
x = [-100, 100] mm
y = [-100, 100] mm
```

The colour of each square is:

```text
hit.time - minimum hit.time in the displayed frame
```

so very large absolute calibrated times do not dominate the colour scale. Press Enter to advance to the next frame, or `q` then Enter to quit.

The colour reference can also be selected explicitly, using the same selector style as `deltat.C`. For example, subtract the trigger tag time from device 200:

```cpp
display("triggered.root", channel_selector_t(9, 200));
```

For hit selectors, the first matching hit in the current frame is used as the reference, scanning trigger, then timing, then Cherenkov hits.

or subtract a fully specified hit:

```cpp
display("triggered.root", field_selector_t(9, 200, 32, -1, -1));
```

For files augmented by `process/bin/timing`, the display can subtract the trained TIMING estimator:

```cpp
display("triggered.timing.root", timing_selection_t("T"));
```

Valid timing selections are `"T"`, `"T0"`, and `"T1"`.

A specific frame can be drawn directly:

```cpp
display("triggered.root", spill_id, frame_index);
display("triggered.root", spill_id, frame_index, channel_selector_t(9, 200));
display("triggered.timing.root", spill_id, frame_index, timing_selection_t("T"));
```

## xymap.C

Integrated Cherenkov occupancy map for triggered frames. It loops over all available frames, counts entries per physical channel using the stored `x/y` positions, and draws one square per channel. The drawn square size defaults to 3.2 mm and can be changed with the optional `pixel_size` argument. The display axes are fixed to:

```text
x = [-100, 100] mm
y = [-100, 100] mm
```

The colour is:

```text
entries in channel / maximum entries in any displayed channel
```

The function is named `draw_map` because a global ROOT macro function named `map` collides with C++ standard-library names in interactive ROOT.

Run:

```bash
root -l -b -q -e '.L macros/example/xymap.C' -e 'draw_map("triggered.root", "xymap.root")'
```

The output ROOT file contains the canvas `cMap` and axis object `hMapAxis`.

## deltat.C

Reads the triggered-frame `frames` tree and fills 2D delta-t histograms using hits contained in each stored frame. Selector values of `-1` are wildcards. All delta-t histograms use 2048 bins over `[-32, 32]`.

The original trigger-relative form finds one reference trigger hit and fills `hDeltaT` for all other hits:

```bash
root -l 'macros/example/deltat.C("frames.root", 9, 200, 32, -1, -1)'
```

The macro also defines selector types with fixed constructor arity:

```cpp
field_selector_t(type, device, fifo, column, pixel)
channel_selector_t(type, channel_or_trigger_device)
```

The exact constructor arity is intentional: missing selector fields are not silently zero-initialised.

Field selectors work for all stored hit categories, including trigger tags:

```cpp
deltat("frames.root",
       field_selector_t(1, 192, -1, -1, -1),
       field_selector_t(1, 192, 0, 0, 0),
       "deltat.root");
```

Channel selectors use the same global channel index as the `hDeltaT` x axis for ALCOR hits:

```cpp
channel_selector_t(1, channel)
```

where:

```cpp
channel = pixel + 4 * column + 32 * (fifo / 4)
        + 256 * (device - 192);
```

For trigger tags, the compact convention is:

```cpp
channel_selector_t(9, device)
```

which is interpreted as:

```cpp
field_selector_t(9, device, 32, -1, -1)
```

The target/reference overloads compute:

```text
delta_t = target.time - reference.time
```

and fill:

```text
hDeltaT       delta_t vs target channel
hDeltaT_spill delta_t vs spill id, 100 bins over [0,100] using the `frames.id` value written by `trigger`; this id is the DAQ START_SPILL counter.
hDeltaT_tdc0  delta_t vs target fine for target hits with tdc == 0
hDeltaT_tdc1  delta_t vs target fine for target hits with tdc == 1
hDeltaT_tdc2  delta_t vs target fine for target hits with tdc == 2
hDeltaT_tdc3  delta_t vs target fine for target hits with tdc == 3
```

For example, compare all ALCOR target channels against trigger tags from device 200:

```cpp
deltat("triggered.root",
       channel_selector_t(1, -1),
       channel_selector_t(9, 200),
       "deltat.root");
```

Files augmented by `process/bin/timing` can also use the trained TIMING estimator as the reference. The timing-reference overloads compute:

```text
delta_t = target.time - timing_reference
```

where `timing_reference_t("T")` uses the combined estimate, `timing_reference_t("T0")` uses TIMING0, and `timing_reference_t("T1")` uses TIMING1. Only frames with `timing_valid == 1` are used. The output filename remains the last argument:

```cpp
deltat("triggered.with_timing.root",
       channel_selector_t(1, -1),
       timing_reference_t("T"),
       "deltat.timing_reference.root");
```

The timing-reference overload also supports optional Cherenkov ring selection. It fits one circle per frame with deterministic RANSAC and keeps only Cherenkov target hits within the configured radial tolerance:

```cpp
deltat("triggered.timing.timing.root",
       channel_selector_t(1, -1),
       timing_reference_t("T"),
       ring_selection_t(8, 64, 256, 3.5, 1., 200.),
       "deltat.ring.root");
```

The `ring_selection_t` arguments are minimum inliers, maximum inliers, RANSAC iterations, radial tolerance in mm, minimum radius, and maximum radius. A frame is accepted only when the fitted ring has an inlier count within the inclusive minimum/maximum range. Ring selection requires `cherenkov_x` and `cherenkov_y` branches, normally added by `coordinator`. The default overload remains unchanged and does not apply ring selection. Ring diagnostics are written as `hRingRadius`, `hRingInliers`, `hRingResidual`, and `hRingCenter`. The latter is the fitted `(x,y)` centre distribution. `hRingResidual` plots the exact selection variable `abs(distance(hit, fitted_center) - fitted_radius)`; hits at or below the configured tolerance are ring-compatible. Since three non-collinear points always define a circle, the inlier and residual distributions should be inspected to reject accidental RANSAC circles in dense or unstructured frames.

The old positional target/reference form remains available as a wrapper:

```cpp
deltat(filename,
       target_type, target_device, target_fifo, target_column, target_pixel,
       reference_type, reference_device, reference_fifo, reference_column, reference_pixel,
       outfilename);
```


`trigger_reader.C` demonstrates the high-level triggered-frame reader. It prints each spill, the number of frames, the number of participating sources from `spill_participation`, and then loops over the trigger/timing/Cherenkov hit collections for each frame.
