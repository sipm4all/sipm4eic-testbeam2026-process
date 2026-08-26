# Cherenkov Ring Finding

This document describes the time-aware Cherenkov ring finder used after the
`timing` stage. The current implementation is `ring-finder-hough`. It searches
for circular or elliptical ring patterns in the Cherenkov hits of each
triggered frame while using both spatial and timing information.

The input must contain synchronized `frames` and `cherenkov` trees. The
Cherenkov tree must provide the per-hit `x`, `y`, and trigger-relative `time`
branches. The output preserves the input trees and adds one `ring` tree entry
for every frame.

## Workflow

The usual workflow is:

```text
decoder -> checker -> process -> trigger -> timing -> ring-finder-hough
```

Run the executable directly with:

```bash
process/bin/ring-finder-hough \
    --input triggered.timing.root \
    --output triggered.timing.ring.root \
    --max-events 1000 \
    --gpu
```

For a processed run, use the workflow script:

```bash
process/scripts/ring-finder.sh \
    --run 20260623-185238 \
    --trigger timing \
    --gpu
```

By default, the script reads:

```text
/data/2026-testbeam/process/<run>/trigger/triggered.<tag>.timing.root
```

and writes:

```text
/data/2026-testbeam/process/<run>/trigger/triggered.<tag>.timing.ring.root
```

The script also supports split-spill processing:

```bash
process/scripts/ring-finder.sh \
    --run RUN \
    --trigger TAG \
    --parallel-spills \
    --jobs 8 \
    --gpu
```

The split files are processed independently and combined with `hadd`. They are
kept by default. Pass `--clean-ring-spills` to remove them after a successful
merge. Existing outputs are not overwritten unless `--overwrite` is supplied.
Use `--input-stage trigger` when the timing-estimator stage has intentionally
been skipped.

## Overall Strategy

Ring finding is deliberately split into localization and parameter estimation:

1. A time-aware RANSAC pass finds approximate circular ring seeds.
2. Each seed defines an independent local Hough search region.
3. The Hough scan evaluates a six-dimensional ellipse model.
4. The strongest maxima are extracted, refined, validated, and deduplicated.
5. Accepted candidates are written to the `ring` tree.

RANSAC is retained because it makes the Hough search practical. It identifies
candidate ring regions, but its fitted circle is not the final ring result.
Different seeds are scanned independently. Their candidate lists are combined
only after the scans, then sorted by score and deduplicated.

With the default RANSAC settings, a frame without a valid RANSAC seed produces
no Hough candidates. Set `--ransac-iterations 0` to disable the localization
stage and scan the configured global bounds instead. That mode is much more
expensive and is mainly useful as a reference or when the seed stage is known
to be inadequate.

## RANSAC Seeding

RANSAC repeatedly chooses three time-compatible Cherenkov hits and constructs
a circular hypothesis. A hypothesis is evaluated using spatial distance from
the circle and time difference from the provisional ring time.

The RANSAC spatial cut is `--ransac-tolerance`. The RANSAC time compatibility
and the local Hough time half-width are controlled by
`--ransac-time-window`. The default values are:

```text
--ransac-iterations       128
--ransac-tolerance           5 mm
--ransac-center-window      10 mm
--ransac-radius-window      10 mm
--ransac-time-window         5 native time units
```

The center and radius windows define the local Hough region around each seed:

```text
x0     = seed_x0     +/- 10 mm
y0     = seed_y0     +/- 10 mm
radius = seed_radius +/- 10 mm
time   = seed_time   +/- 5 native time units
```

The RANSAC tolerance is independent of the Hough spatial resolution. The
RANSAC pass uses circles because it is a fast way to localize a possible ring;
the subsequent Hough scan can select a nonzero eccentricity.

## Six-Dimensional Hough Scan

For every RANSAC seed, the Hough search scans:

```text
(x0, y0, radius, ring_time, eccentricity, phi)
```

The default local grid is:

```text
x0:           21 bins, seed_x0 +/- 10 mm, step 1 mm
y0:           21 bins, seed_y0 +/- 10 mm, step 1 mm
radius:       21 bins, seed_radius +/- 10 mm, step 1 mm
time:         11 bins, seed_time +/- 5 native units, step 1 native unit
eccentricity: 10 bins, 0..0.9, step 0.1
phi:          20 bins, 0..pi, step pi/19
```

Thus the normal per-seed scan contains approximately:

```text
21 x 21 x 21 x 11 x 10 x 20 = 20,404,200 Hough cells
```

The first four ranges are localized around the RANSAC seed. The eccentricity
and orientation ranges are scanned for every localized seed. The command-line
options controlling the ellipse dimensions are:

```text
--min-e       minimum eccentricity
--max-e       maximum eccentricity
--e-step      eccentricity step
--min-phi     minimum orientation in radians
--max-phi     maximum orientation in radians
--phi-step    orientation step in radians
```

Eccentricity must be in `[0,1)`. Orientations are equivalent modulo `pi`
because an ellipse has no directed major axis. The default `ring_e=0` model is
exactly a circle; `ring_phi` is then physically irrelevant.

## Ellipse Model

`ring_r` is the semi-major axis. The semi-minor axis is calculated as:

```text
b = ring_r * sqrt(1 - ring_e^2)
```

The ellipse is rotated by `ring_phi` radians. For a hit, the implementation
calculates the radial residual between the hit's distance from `(x0,y0)` and
the ellipse radius at that direction. The Hough model therefore uses the same
ellipse residual in the CPU and CUDA paths.

The Hough score is the sum of a spatial Gaussian and a time Gaussian over the
event hits. Conceptually, each hit contributes:

```text
spatial_weight = Gaussian(ellipse_radial_residual / spatial_resolution)
time_weight    = Gaussian((hit_time - ring_time) / time_resolution)
vote           = spatial_weight * time_weight
```

The Gaussian is evaluated only within four standard deviations. The controls
are:

```text
--spatial-resolution   Gaussian width in mm, default 1.5
--time-resolution      Gaussian width in native time units, default 1
```

These are measurement widths, not Hough-grid steps. The grid steps are
controlled independently by `--x0-step`, `--y0-step`, `--radius-step`, and
`--t-step`.

## GPU Implementation

With `--gpu`, the CUDA backend assigns one GPU thread to each complete
six-dimensional Hough cell. Each thread loops over the event hits and writes a
single score, so the vote calculation does not require accumulator atomics.

The coordinate maps and accumulator are allocated once and reused when
successive RANSAC seeds have the same grid dimensions. Different seed origins
update the coordinate maps without requiring a new allocation. The complete
accumulator is copied back once after the vote kernel. The common C++ code
keeps only a bounded top-score candidate pool rather than sorting every Hough
cell, then performs peak suppression, candidate merging, interpolation, and
validation.

The CUDA memory requirement is proportional to the number of local Hough cells.
Reducing the RANSAC windows or increasing grid steps reduces both memory use and
runtime. A full global scan with `--ransac-iterations 0` can be much larger
than the normal localized scan.

## Peak Extraction and Refinement

Hough maxima are selected with local suppression in all six scanned dimensions.
This prevents a broad maximum from producing many neighboring copies of the
same candidate. Maxima from separate RANSAC seed scans are then combined.

The existing sub-grid interpolation refines the four continuous coordinates:

```text
x0, y0, radius, ring_time
```

It fits a local quadratic surface to neighboring accumulator values. The
eccentricity and `phi` values remain the corresponding discrete Hough maximum.
If the local quadratic is not well formed, the discrete maximum is retained.

## Candidate Validation

Every Hough candidate is checked against the original event hits. A hit is an
inlier only if both conditions hold:

```text
abs(ellipse_radial_residual) <= 4 * spatial_resolution
abs(hit_time - ring_time)   <= 4 * time_resolution
```

The candidate is accepted only if it has at least `--min-inliers` inliers. The
default is 8. The stored `ring_ninliers` is the number of validated inliers.

The candidate list is processed in descending Hough score. Hits are not
removed after accepting a ring, so multiple rings may share hits. A new ring
is rejected only when its shared-hit count with an already accepted ring is
larger than:

```text
floor(--max-shared-fraction *
      min(number_of_new_ring_inliers, number_of_previous_ring_inliers))
```

The default `--max-shared-fraction` is `0.5`. This means that no more than half
of the smaller ring's inlier set may be shared. Previously accepted rings are
never removed because a later candidate overlaps them.

## Output Tree

The output contains one `ring` entry for every input frame. The branch layout is:

```text
nring                         UChar_t
ring_x0[nring]                float
ring_y0[nring]                float
ring_r[nring]                 float
ring_e[nring]                 float
ring_phi[nring]               float
ring_time[nring]              float
ring_ninliers[nring]          unsigned short
```

The entry index is the frame index shared by the `frames`, `trigger`, `timing`,
and `cherenkov` trees. For ring `i`:

```text
(ring_x0[i], ring_y0[i])       ellipse center in mm
ring_r[i]                      semi-major axis in mm
ring_e[i]                      eccentricity
ring_phi[i]                    orientation in radians
ring_time[i]                   ring time in native units
ring_ninliers[i]               validated inlier count
```

Multiple rings may be stored in one frame. `nring=0` is a valid result when no
candidate passes validation.

## Useful Commands

The default Hough parameters can be inspected with:

```bash
process/bin/ring-finder-hough --help
```

Use `--max-events N` to process only the first `N` frames. The default is
`-1`, which processes all frames. This is useful for smoke tests and timing
benchmarks without changing the input file.

A circle-only reference run can be made by collapsing the ellipse dimensions:

```bash
process/bin/ring-finder-hough \
    --input triggered.timing.root \
    --output triggered.circle.root \
    --min-e 0 --max-e 0 --e-step 1 \
    --min-phi 0 --max-phi 0 --phi-step 1 \
    --gpu
```

For a faster diagnostic test, reduce the local seed windows and process a
small input file or split spill. For production, use the default RANSAC-local
ellipse search unless the Hough ranges and resolutions have been deliberately
changed.
