# Main Processing Workflow

This document describes the intended end-to-end workflow driven by
`process/scripts/main.sh`. The workflow is spill-preserving: each stage keeps
its output split by spill, and no final `hadd` file is required.

## Stage order

For each run, the processing stages are:

```text
decode
  -> calibrate
  -> coordinate
  -> sort
  -> after-pulse suppression
  -> merge
  -> trigger
  -> timing estimator
  -> ring finder
  -> filter
```

The first per-FIFO stages can run in parallel. Device and run merging happen
only after the required input FIFOs have completed. Trigger, timing, ring, and
filter processing can then run independently for each spill, subject to the
workflow's configured parallelism.

## File layout and names

The expected files are:

```text
<run>/
├── decoded/
│   └── alcdaq.fifo_<fifo>.root
│       kept: decoded input for reproducibility
│
├── process/
│   ├── calibrated.fifo_<fifo>.root
│   │   deleted after coordinated output is complete
│   ├── coordinated.fifo_<fifo>.root
│   │   deleted after sorted output is complete
│   ├── sorted.fifo_<fifo>.root
│   │   deleted after after-pulse suppression is complete
│   ├── aps.fifo_<fifo>.root
│   │   deleted after merged spill output is complete
│   └── merged.<device>.spill_<spill>.root
│       kept: device-level merged spill input for triggering
│
└── trigger/
    ├── triggered.<tag>.spill_<spill>.root
    │   deleted after successful timing production
    ├── timing.<tag>.spill_<spill>.root
    │   deleted after successful ring production
    ├── rings.<tag>.spill_<spill>.root
    │   kept: reconstructed ring information
    └── filtered.<tag>.spill_<spill>.root
        kept: final selected analysis data
```

The exact directory prefix is configured by the processing scripts. The
filenames above describe the stage identity and are deliberately retained in
the output names so that a file can be identified without inspecting its
contents.

## Cleanup policy

Cleanup is performed only after the downstream output has been created
successfully. By default:

| File | Default action |
| --- | --- |
| `decoded/alcdaq.fifo_<fifo>.root` | Keep |
| `process/calibrated.fifo_<fifo>.root` | Delete after coordination |
| `process/coordinated.fifo_<fifo>.root` | Delete after sorting |
| `process/sorted.fifo_<fifo>.root` | Delete after after-pulse suppression |
| `process/aps.fifo_<fifo>.root` | Delete after merging |
| `process/merged.<device>.spill_<spill>.root` | Keep |
| `trigger/triggered.<tag>.spill_<spill>.root` | Delete after timing estimation |
| `trigger/timing.<tag>.spill_<spill>.root` | Delete after ring finding |
| `trigger/rings.<tag>.spill_<spill>.root` | Keep |
| `trigger/filtered.<tag>.spill_<spill>.root` | Keep |

An interrupted or failed stage must not delete its input. Existing outputs are
also skipped unless `--overwrite` is explicitly requested, so a stopped run
can be resumed without destroying completed work.

## Main configuration

At this stage `main.sh` accepts a run name or a newline-separated run-list as
its positional argument. Calibration, clock-correction, trigger, filter,
window, run type, and GPU settings are configured near the top of the script.
The intended invocation is therefore:

```bash
process/scripts/main.sh 20260623-185238
process/scripts/main.sh run.list
```

Individual stages can be selected with repeatable `--step` options. This is
useful when decoding, checking, or rerunning one later stage independently:

```bash
process/scripts/main.sh RUN --step decoder
process/scripts/main.sh RUN --step checker
process/scripts/main.sh RUN --step process
process/scripts/main.sh RUN --step trigger
process/scripts/main.sh RUN --step timing
process/scripts/main.sh RUN --step ring
process/scripts/main.sh RUN --step filter
process/scripts/main.sh RUN --step timing --step ring --step filter
```

With no `--step`, all stages run in the order listed above. Selecting a stage
does not implicitly run its prerequisites; those inputs must already exist.
`--overwrite` can be added when an existing output should be regenerated.

The individual workflows remain available when a stage must be rerun or
debugged separately:

```text
decoder.sh       decoded per-FIFO ROOT files
checker.sh       consistency reports
process.sh       calibrated/coordinated/sorted/APS and merged spill files
trigger.sh       triggered spill files
timing.sh        timing-estimator spill files
ring-finder.sh   ring spill files
filter.sh        filtered spill files
```

All stages must preserve spill isolation. A spill from one input must never be
combined with a different spill during processing or analysis.
