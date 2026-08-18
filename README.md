# SiPM4EIC Test Beam 2026 Processing

This repository contains ROOT-based tools for processing SiPM4EIC 2026 test-beam data after decoding. The processing chain is organized as:

```text
raw decode -> optional check -> clean -> calibrate -> coordinate -> sort -> merge -> trigger/frame -> timing estimator -> analysis macros
```

The repository is intentionally split between compiled processing programs and ROOT analysis macros:

```text
process/source/          compiled C++ programs and build files
process/scripts/         batch workflow scripts
process/config/          calibration and trigger configuration files
macros/example/          small ROOT macro examples
macros/lib/              reusable ROOT macro helper headers
```

## Quick Build

From the repository root:

```bash
cmake -S process/source -B process/build
cmake --build process/build -j
cmake --install process/build
```

After install, this creates the processing executables in `process/bin/`:

```text
decoder
dcalib
cleaner
checker
calibrator
coordinator
sorter
after-pulse-suppressor
merger
trigger
timing
```


## Processing Layout

The workflow scripts use one run directory under `/data/2026-testbeam/process` and keep device-local products separated by stage:

```text
/data/2026-testbeam/process/<run>/<device>/decoded/   decoded ROOT files from decoder.sh
/data/2026-testbeam/process/<run>/<device>/check/     checker reports for decoded FIFO files
/data/2026-testbeam/process/<run>/<device>/dcalib/    TDC-calibration products from dcalib.sh
/data/2026-testbeam/process/<run>/<device>/process/   calibrated/sorted/AP-suppressed processing products
/data/2026-testbeam/process/<run>/check/              run-level checker reports
/data/2026-testbeam/process/<run>/process/            run-level merged split-spill ROOT files
/data/2026-testbeam/process/<run>/trigger/            triggered-frame ROOT files and timing-estimator augmented files
```

Downstream workflows read the decoded ROOT files from the `decoded/` directory. They no longer read decoded files directly from `/data/2026-testbeam/actual`.

## Main Components

- `decoder`: converts raw per-FIFO `.dat` files into decoded ROOT `alcor` trees with strict spill validation.
- `dcalib`: derives TDC `off/iif` calibration rows from decoded data.
- `cleaner`: clones the `alcor` tree while dropping malformed ALCOR hits, currently invalid `fifo`/`column` combinations.
- `checker`: scans a per-FIFO `alcor` tree and writes an ASCII `.check` sanity report.
- `calibrator`: creates or updates the calibrated `time` branch in the decoded `alcor` tree.
- `coordinator`: creates or updates the spatial `x` and `y` branches in the calibrated `alcor` tree.
- `sorter`: sorts each single-lane stream in calibrated time.
- `after-pulse-suppressor`: suppresses close repeated ALCOR hits per channel.
- `merger`: merges sorted lane/device streams while preserving spill boundaries and writes spill participation metadata.
- `trigger`: creates triggered frame trees using declarative trigger configurations and propagates spill participation metadata.
- `timing`: reads triggered frame trees and adds per-frame TIMING estimator branches.
- `macros/lib/trigger_reader.h`: header-only helper for reading triggered frame output.

Run `process/scripts/checker.sh` first when you want a non-destructive sanity pass over decoded per-FIFO files before launching calibration or full processing. The checker writes run-level `*.good-fifos.list` and `*.bad-fifos.list` diagnostic files that show which FIFOs are consistent with the run-level spill structure. Merged ROOT files contain the usual `alcor` tree plus a `spill_participation` tree with one entry per merged spill and the `(device,fifo)` sources that contributed to that spill.

## User Guides

- `docs/reconstruction.md`: standard physics reconstruction workflow, from raw decoding through triggered-frame output.
- `docs/calibration.md`: TDC calibration notes and calibration-check workflows.
- `docs/todo.md`: design notes and deferred technical improvements.

Generated ROOT files, logs, build products, and local binaries are ignored by git.
