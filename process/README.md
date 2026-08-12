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

Most compiled programs operate on ROOT files; `decoder` is the raw `.dat` to ROOT entry point. The common raw input tree name is `alcor`. Triggered output files contain a `frames` tree with one entry per spill and flattened frame contents split into trigger, timing, and Cherenkov collections. The spill entry `id` is copied from the DAQ START_SPILL word `counter`, so it preserves the original spill number even when spills are processed as separate files. The merger enforces spill alignment by requiring all input streams to have the same START/END spill counters at every boundary.


Workflow scripts use `/data/2026-testbeam/process/<run>/` as the common run workspace. Device-local products are separated by stage:

```text
<device>/decoded/   decoded ROOT files and decoder/checker reports
<device>/dcalib/    TDC-calibration products from dcalib.sh
<device>/process/   calibrated, sorted, AP-suppressed, and device-merged processing products
```

`decoder.sh` is the only workflow that reads raw files from `/data/2026-testbeam/actual/<run-type>/<run>/...`. `checker.sh`, `dcalib.sh`, and `process.sh` read decoded ROOT files from `<device>/decoded/`.

Directory contents:

```text
source/                 C++ programs and CMake build file
scripts/                shell scripts for larger batch workflows
config/calibration/     timing calibration files
config/trigger/         trigger-definition files
```

The `cleaner` stage removes malformed stream words before calibration, currently ALCOR hits with invalid `fifo`/`column` combinations. The key design point is that low-level timing calibration happens once, early, in `calibrator`. Downstream tools prefer the calibrated `time` branch when present and only use the old nominal timing expression for legacy uncalibrated files.

The optional `checker` program scans decoded per-FIFO files before processing and writes ASCII `.check` reports with word counts and spill-counter consistency checks. It is read-only: it does not create ROOT data products and does not affect the normal processing chain.

The `decoder` stage is intentionally conservative. It first finds START_SPILL and then requires a matching END_SPILL with the same counter. ALCOR payload words are accepted only if their decoded column is valid for the FIFO: `2 * (fifo % 4)` or `2 * (fifo % 4) + 1`. Normal spill markers use `fine = 0`. If a completed spill has more decoding errors than `--allowed-spill-errors`, the output still contains the START/END spill markers but the spill payload is omitted and both markers use `fine = 1`. DAQ-readout-suppressed records marked by `0xdeadbeef` are counted in the decoder summary and skipped; they are never written as fake spill markers.
