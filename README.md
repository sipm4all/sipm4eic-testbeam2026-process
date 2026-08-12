# SiPM4EIC Test Beam 2026 Processing

This repository contains ROOT-based tools for processing SiPM4EIC 2026 test-beam data after decoding. The processing chain is organized as:

```text
raw decode -> optional check -> clean -> calibrate -> sort -> merge -> trigger/frame -> analysis macros
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
sorter
after-pulse-suppressor
merger
trigger
```


## Processing Layout

The workflow scripts use one run directory under `/data/2026-testbeam/process` and keep device-local products separated by stage:

```text
/data/2026-testbeam/process/<run>/<device>/decoded/   decoded ROOT files from decoder.sh
/data/2026-testbeam/process/<run>/<device>/dcalib/    TDC-calibration products from dcalib.sh
/data/2026-testbeam/process/<run>/<device>/process/   calibrated/sorted/AP-suppressed processing products
```

Downstream workflows read the decoded ROOT files from the `decoded/` directory. They no longer read decoded files directly from `/data/2026-testbeam/actual`.

## Main Components

- `decoder`: converts raw per-FIFO `.dat` files into decoded ROOT `alcor` trees with strict spill validation.
- `dcalib`: derives TDC `off/iif` calibration rows from decoded data.
- `cleaner`: clones the `alcor` tree while dropping malformed ALCOR hits, currently invalid `fifo`/`column` combinations.
- `checker`: scans a per-FIFO `alcor` tree and writes an ASCII `.check` sanity report.
- `calibrator`: creates or updates the calibrated `time` branch in the decoded `alcor` tree.
- `sorter`: sorts each single-lane stream in calibrated time.
- `after-pulse-suppressor`: suppresses close repeated ALCOR hits per channel.
- `merger`: merges sorted lane/device streams while preserving spill boundaries.
- `trigger`: creates triggered frame trees using declarative trigger configurations.
- `macros/lib/trigger_reader.h`: header-only helper for reading triggered frame output.

Run `process/scripts/checker.sh` first when you want a non-destructive sanity pass over decoded per-FIFO files before launching the full processing workflow.

Generated ROOT files, logs, build products, and local binaries are ignored by git.
