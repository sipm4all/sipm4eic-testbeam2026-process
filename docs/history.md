# Development Notes and Design History

This file records the main design decisions made while organizing the SiPM4EIC test-beam processing repository. It is not intended as user-facing documentation; the top-level `README.md` remains the entry point for users.

## Repository Organization

The repository was organized into three main areas:

```text
process/    compiled processing programs, scripts, and configs
macros/     ROOT macro examples and reusable macro headers
docs/       development notes and historical context
```

The processing directory is further split into:

```text
process/source/              C++ sources and CMake build file
process/scripts/             workflow scripts
process/config/calibration/  calibration configuration files
process/config/trigger/      trigger configuration files
```

The macro directory is split into:

```text
macros/example/              ROOT example macros
macros/lib/                  reusable ROOT macro helper headers
```

## Files Intentionally Not Imported

The following legacy or local files were intentionally not copied into the new repository layout:

```text
read.C
sorter.C
deltat.C
deltat_trigger.C
```

`deltat_trigger_frames.C` was copied as:

```text
macros/example/deltat.C
```

because it analyzes the current triggered-frame output format.

Generated files were also excluded:

```text
*.root
build/
compiled binaries
ROOT ACLiC artifacts
backup files
logs
plots
```

## Processing Pipeline

The intended pipeline is:

```text
decode -> calibrate -> sort -> after-pulse suppress -> merge -> trigger/frame -> analysis macros
```

A dedicated calibration stage was introduced before sorting so downstream programs can use one calibrated time coordinate.

## Calibration Stage

`calibrator` writes an authoritative `time` branch into the `alcor` tree.

For ALCOR hits:

```text
phase = off + iif * fine
time = coarse + 32768 * rollover - phase - channel_offset
```

For trigger tags:

```text
time = coarse + 32768 * rollover - trigger_offset
```

The calibration file format is whitespace-separated with sections:

```text
[TDC]
[CHANNEL]
[TRIGGER]
```

Address fields support exact integer values or `*` wildcards. The most-specific matching rule wins. Equal-specificity ambiguities are errors.

## Shared Compiled-Code Headers

Two shared headers were introduced for readability and consistency:

```text
process/source/data_word.h
process/source/calibration.h
```

`data_word.h` centralizes the ROOT `alcor` word structure, branch binding, word-type helpers, and nominal fallback timing.

`calibration.h` centralizes calibration parsing, rule matching, ambiguity detection, and cached concrete lookup.

## Downstream Time Handling

Compiled downstream programs were updated to prefer an existing calibrated `time` branch. If no `time` branch is present, they fall back to the older nominal formula:

```text
coarse + 32768 * rollover - 0.0157 * fine
```

This fallback is only for legacy uncalibrated files.

## Trigger Configuration

Trigger configuration files remain declarative and are kept under:

```text
process/config/trigger/
```

When using calibrated files, selector `time_offset` values should normally be `0.` to avoid applying the same physical timing correction twice.

The trigger system supports:

```text
seeded triggers
self-coincidence triggers
require/veto conditions
hit or distinct-channel multiplicity
selector wildcard/exact/range/set matching
```

## Triggered Frame Output Reader

The triggered-frame ROOT output stores one tree entry per spill, with common `nframes` and three split flattened collections:

```text
trigger
timing
cherenkov
```

The first version stored only raw hit fields. After introducing the calibration stage, calibrated hit `time` was also persisted in each category as `trigger_time`, `timing_time`, and `cherenkov_time`, because it is now a calibrated data product rather than a temporary helper.

To hide this layout from analysis users, a header-only reader was added:

```text
macros/lib/trigger_reader.h
```

A small example macro is provided:

```text
macros/example/trigger_reader.C
```
