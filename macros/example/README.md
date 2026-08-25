# Example ROOT Macros

This directory contains small ROOT macros showing how to inspect processed data.

## trigger_reader.C

Demonstrates the header-only triggered-frame reader API:

```bash
root -l 'macros/example/trigger_reader.C("frames.root")'
```

It iterates over spills and frames and prints the number of trigger, timing, and Cherenkov hits in each frame.

## display.C

Interactive Cherenkov event display for triggered frames. It draws one frame at a time using the stored hit `x/y` positions and one square per channel, rather than a regular `TH2D` bin grid. This is important because the channel centres are not assumed to lie on a perfectly regular ROOT histogram grid. The complete selection and display API is documented in [DISPLAY.md](DISPLAY.md).

Run:

```bash
root -l 'macros/example/display.C("triggered.root")'
```

It displays only Cherenkov hits. The drawn square size is controlled by the
top-level `pixel_size` constant in `display.C` (currently 3.0 mm). The display
axes are fixed to:

```text
x = [-100, 100] mm
y = [-100, 100] mm
```

If the input file contains the `ring` tree produced by `process/bin/ring-finder`,
`display.C` overlays the found circles or ellipses and can mark the associated
hits. Files without a ring tree remain supported and show only the Cherenkov
hit display.

The colour of each square is:

```text
hit.time - minimum hit.time in the displayed frame
```

so very large absolute calibrated times do not dominate the colour scale. In
the interactive loop, press Enter to advance, `s` to save a PNG, or `q` then
Enter to quit. See [DISPLAY.md](DISPLAY.md) for reference selectors, ring
cuts, and the optional trigger-hit requirement.

The colour reference can also be selected explicitly. References and filters
are polymorphic objects; there is at most one reference and all filters are
combined with logical AND. For example, subtract the trigger tag time from
device 200:

```cpp
display("triggered.root",
        std::make_shared<channel_reference_t>(
            channel_selector_t(9, 200)));
```

For hit selectors, the first matching hit in the current frame is used as the reference, scanning trigger, then timing, then Cherenkov hits.

or subtract a fully specified hit:

```cpp
display("triggered.root",
        std::make_shared<field_reference_t>(
            field_selector_t(9, 200, 32, -1, -1)));
```

For files augmented by `process/bin/timing`, the display can subtract the trained TIMING estimator:

```cpp
display("triggered.timing.root",
        std::make_shared<timing_reference_t>("T"));
```

Valid timing selections are `"T"`, `"T0"`, and `"T1"`.

To require a trigger from a specific device, pass a `trigger_selection_t`:

```cpp
display("triggered.timing.root",
        std::make_shared<timing_reference_t>("T"),
        {std::make_shared<trigger_selection_t>(200)});
```

Use `trigger_selection_t(-1)` to require a trigger from any device. The
trigger selector is optional and omitted by default.

A specific frame can be drawn directly by setting the final two arguments to
the target spill and frame. The first four optional arguments are the
reference, filter vector, starting spill, and starting frame:

```cpp
display("triggered.root", nullptr, {},
        std::numeric_limits<int>::min(), 0,
        spill_id, frame_index);
display("triggered.timing.root",
        std::make_shared<timing_reference_t>("T"),
        {std::make_shared<trigger_selection_t>(200)},
        std::numeric_limits<int>::min(), 0,
        spill_id, frame_index);
```

To display every frame without a reference or filters, use either
`display("ring.root")` or `display("ring.root", {})`.

## ring_analysis.C

`ring_analysis.C` performs a non-interactive analysis of the stored rings. It
uses the same polymorphic `reference_ptr_t` and `selection_ptr_t` arguments as
`display.C`; all supplied selections are combined with logical AND. For every
ring in a selected frame it compares every finite Cherenkov hit with the ring
time and geometry, filling:

```text
hDeltaTRing             hit time - ring time, in ns
hRingDistance           signed ellipse radial residual, in mm
hDeltaTRingVsDistance   two-dimensional time/signed-residual correlation
hDeltaTRingMatched      same time distribution for display-style matched hits
hRingDistanceMatched    same signed residual distribution for matched hits
hDeltaTRingUnmatched    time distribution for hits not matched to the ring
hRingDistanceUnmatched  signed residual distribution for hits not matched to the ring
hDeltaTRingVsDistanceUnmatched  two-dimensional unmatched-hit correlation
```

The matched-hit definitions use a spatial tolerance of 5 mm and a time
window of +/-2 native time units, matching the ring overlay in `display.C`.
The event-level ring selection can be combined with a timing reference or a
trigger requirement:

```cpp
ring_analysis("triggered.timing.rings.root",
              std::make_shared<timing_reference_t>("T"),
              {std::make_shared<trigger_selection_t>(200),
               std::make_shared<ring_selection_t>(-100, 100, -100, 100,
                                                  20, 45, 1, 2)},
              "ring_analysis.root");
```

The reference is used as an additional frame requirement; `#Delta t` itself
is always measured relative to each stored `ring_time`.

The residual is signed: positive values are outside the fitted ellipse and
negative values are inside it. The matching decision still uses its absolute
value. The unmatched histograms are the useful independent diagnostic: the ring
finder uses the matched hits to fit the ring, so a peak at zero in the matched
residual is expected by construction. Unmatched hits are not accepted by the
same spatial-and-time association (`5` mm and `+/-2` native units).

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

Reads synchronized triggered frames and fills fixed-range delta-t histograms.
`hDeltaT` uses 2048 bins over `[-32,32]`. The spill and TDC diagnostics use
128 bins over `[-2,2]`, preserving the same `0.03125` native-unit bin width.

The macro uses the shared polymorphic API from
`macros/lib/frame_selection.h`:

- one reference object
- one or more target objects
- zero or more event selections

All target objects are combined with logical AND for each hit. All event
selections are also combined with logical AND for each frame. The output
filename remains the last argument.

The basic API is:

```cpp
deltat(filename,
       targets,
       reference,
       selections,
       outfilename);
```

For the common case of one target, the vector can be created internally:

```cpp
deltat(filename,
       std::make_shared<field_target_t>(device, fifo, column, pixel),
       reference,
       selections,
       outfilename);
```

Targets include:

```cpp
field_target_t(device, fifo, column, pixel)
fifo_target_t(device, fifo)
channel_target_t(channel)
trigger_target_t(device)
ring_target_t(min_inliers, max_inliers, tolerance,
              min_radius, max_radius)
```

For example, compare all type-1 ALCOR hits against trigger tags from device
200:

```cpp
deltat("triggered.root",
       {
           std::make_shared<channel_target_t>(-1)
       },
       std::make_shared<trigger_reference_t>(200),
       {},
       "deltat.trigger.root");
```

Trigger targets use the explicit type-9 target:

```cpp
trigger_target_t(device)
```

which selects trigger tags from that device. The target constructors imply
the type: `channel_target_t`, `fifo_target_t`, and `field_target_t` select
type-1 ALCOR hits, while `trigger_target_t` selects type-9 hits.

`fifo_target_t(device, fifo)` selects all columns and pixels from one type-1
FIFO. As with the field target, `-1` can be used as a wildcard address value.

References include:

```cpp
field_reference_t(field_selector_t(...))
channel_reference_t(channel_selector_t(...))
trigger_reference_t(device)
timing_reference_t("T")
timing_reference_t("T0")
timing_reference_t("T1")
```

For a timing reference:

```cpp
deltat("triggered.timing.root",
       {
           std::make_shared<channel_target_t>(-1)
       },
       std::make_shared<timing_reference_t>("T"),
       {},
       "deltat.timing.root");
```

Only frames with a valid requested timing estimate are used.

A target vector is an AND expression. This selects type-1 hits that are also
inside an accepted stored ring:

```cpp
deltat("triggered.timing.rings.root",
       {
           std::make_shared<channel_target_t>(-1),
           std::make_shared<ring_target_t>(8, 64, 3.5, 1., 200.)
       },
       std::make_shared<timing_reference_t>("T"),
       {},
       "deltat.ring.root");
```

The ring target reads the `ring` tree produced by `ring-finder`. It does
not run RANSAC or any other ring reconstruction. Its arguments are:

```text
minimum inliers
maximum inliers
spatial hit-to-ring tolerance in mm
minimum radius
maximum radius
```

The stored ring geometry, including ellipse parameters and ring time, is used
for hit matching.

Event-level filters are passed separately. For example, require a trigger from
device 200 and a frame-level ring cut:

```cpp
deltat("triggered.timing.rings.root",
       {
           std::make_shared<channel_target_t>(-1),
           std::make_shared<ring_target_t>(8, 64, 3.5, 1., 200.)
       },
       std::make_shared<timing_reference_t>("T"),
       {
           std::make_shared<trigger_selection_t>(200),
           std::make_shared<ring_selection_t>(
               -100., 100., -100., 100., 10., 80., 1, 2)
       },
       "deltat.selected.root");
```

The output contains:

```text
hDeltaT       delta_t versus target global channel
hDeltaT_spill delta_t versus spill id, with 100 bins over [0,100] and
              delta-t range [-2,2]
hDeltaT_tdc0  delta_t versus target fine for TDC 0, range [-2,2]
hDeltaT_tdc1  delta_t versus target fine for TDC 1, range [-2,2]
hDeltaT_tdc2  delta_t versus target fine for TDC 2, range [-2,2]
hDeltaT_tdc3  delta_t versus target fine for TDC 3, range [-2,2]
hDeltaT_spill_deviceX_fifoY
              delta_t versus spill for every encountered type-1
              `(device,fifo)` source
```

The per-FIFO spill histograms are filled by default for all type-1 hits in
frames with a valid reference time. They are created lazily, so absent sources
do not produce empty histograms. Their normalization is independent of the
target histogram and uses the number of frames with a valid reference. This
is useful for scanning clock transitions without first selecting suspicious
FIFOs manually.

When a `ring_target_t` is present, ring diagnostics are also written:
`hRingRadius`, `hRingInliers`, `hRingResidual`, and `hRingCenter`.
These describe the stored rings and their hit residuals; they are not results
of a new fit performed by `deltat.C`.


`trigger_reader.C` demonstrates the high-level triggered-frame reader. It prints each spill, the number of frames, the number of participating sources from `spill_participation`, and then loops over the trigger/timing/Cherenkov hit collections for each frame.
