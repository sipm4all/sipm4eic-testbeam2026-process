# Processing Pipeline

The `process/` directory contains the compiled data-processing chain, workflow scripts, and configuration files.

Recommended order:

```text
raw per-FIFO DAT file
  -> decoder
  -> decoded ROOT file
  -> optional checker
  -> cleaner
  -> calibrator
  -> sorter
  -> after-pulse-suppressor
  -> merger
  -> trigger
  -> triggered frame ROOT files
```

Most compiled programs operate on ROOT files; `decoder` is the raw `.dat` to ROOT entry point. The common raw input tree name is `alcor`. Triggered output files contain a `frames` tree with one entry per spill and flattened frame contents split into trigger, timing, and Cherenkov collections. The spill entry `id` is copied from the DAQ START_SPILL word `counter`, so it preserves the original spill number even when spills are processed as separate files.

Directory contents:

```text
source/                 C++ programs and CMake build file
scripts/                shell scripts for larger batch workflows
config/calibration/     timing calibration files
config/trigger/         trigger-definition files
```

The `cleaner` stage removes malformed stream words before calibration, currently ALCOR hits with invalid `fifo`/`column` combinations. The key design point is that low-level timing calibration happens once, early, in `calibrator`. Downstream tools prefer the calibrated `time` branch when present and only use the old nominal timing expression for legacy uncalibrated files.

The optional `checker` program scans decoded per-FIFO files before processing and writes ASCII `.check` reports with word counts and spill-counter consistency checks. It is read-only: it does not create ROOT data products and does not affect the normal processing chain.

The `decoder` stage is intentionally conservative. It first finds START_SPILL and then requires a matching END_SPILL with the same counter. ALCOR payload words are accepted only if their decoded column is valid for the FIFO: `2 * (fifo % 4)` or `2 * (fifo % 4) + 1`. If a completed spill has more decoding errors than `--allowed-spill-errors`, the output still contains the START/END spill markers but the spill payload is omitted. Suppressed spills are tagged by setting `fine = 1` on both START_SPILL and END_SPILL; normal spill markers have `fine = 0`.
