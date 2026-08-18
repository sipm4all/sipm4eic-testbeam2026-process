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
  -> coordinator
  -> sorter
  -> after-pulse-suppressor
  -> merger
  -> trigger
  -> triggered frame ROOT files
```

Most compiled programs operate on ROOT files; `decoder` is the raw `.dat` to ROOT entry point. The common raw input tree name is `alcor`. Triggered output files contain a `frames` tree with one entry per spill and flattened frame contents split into trigger, timing, and Cherenkov collections. They also propagate the input `spill_participation` tree so analyses can recover which `(device,fifo)` sources contributed to each triggered spill. The spill entry `id` is copied from the DAQ START_SPILL word `counter`, so it preserves the original spill number even when spills are processed as separate files. The merger aligns spills by the START_SPILL/END_SPILL counter. If one input stream is missing a suppressed spill, that stream is absent from that merged spill rather than forcing all streams to fail synchronization. The output still collapses duplicate spill markers in the `alcor` tree, and also writes a `spill_participation` tree with one entry per merged spill and the contributing `(device,fifo)` sources.


Workflow scripts use `/data/2026-testbeam/process/<run>/` as the common run workspace. Device-local products are separated by stage:

```text
<device>/decoded/   decoded ROOT files and decoder summaries
<device>/check/     checker reports for decoded FIFO files
<device>/dcalib/    TDC-calibration products from dcalib.sh
<device>/process/   calibrated, sorted, AP-suppressed, and device-merged processing products
check/              run-level checker reports and FIFO diagnostic lists
process/            run-level merged split-spill ROOT files from process.sh
trigger/            triggered-frame ROOT files from trigger.sh
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

`coordinator` is the analogous geometry stage. It preserves the `alcor` tree and adds spatial `x` and `y` branches for ALCOR hits. TIMING hits from device 200 use the TIMING `eo2do` mapping and a 3.5 mm pitch. Cherenkov hits use the TESTBEAM2026 PDU/matrix mapping inherited from the analysis `mapping.h` logic, with the decoded electronics channel `pixel + 4 * column`.

The optional `checker` program scans decoded per-FIFO files before processing and writes ASCII `.check` reports with word counts and spill-counter consistency checks. It is read-only with respect to ROOT data products, but it also writes run-level `*.good-fifos.list` and `*.bad-fifos.list` diagnostic files that show which FIFOs are consistent with the run-level spill structure.

The `decoder` stage is intentionally conservative. It first finds START_SPILL and then requires a matching END_SPILL with the same counter before writing that spill to the output ROOT tree. ALCOR payload words are accepted only if their decoded column is valid for the FIFO: `2 * (fifo % 4)` or `2 * (fifo % 4) + 1`. Spill markers written by the decoder use `fine = 0`. If a completed spill has more decoding errors than `--allowed-spill-errors`, the entire spill is suppressed and reported in the decoder summary. If EOF arrives before the matching END_SPILL, the buffered incomplete spill is discarded and reported in the decoder summary. DAQ-readout-suppressed records marked by `0xdeadbeef` are counted in the decoder summary and skipped; they are never written as fake spill markers.
