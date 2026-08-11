# SiPM4EIC Test Beam 2026 Processing

This repository contains ROOT-based tools for processing SiPM4EIC 2026 test-beam data after decoding. The processing chain is organized as:

```text
decode -> clean -> calibrate -> sort -> merge -> trigger/frame -> analysis macros
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
dcalib
cleaner
calibrator
sorter
after-pulse-suppressor
merger
trigger
```

## Main Components

- `dcalib`: derives TDC `off/iif` calibration rows from decoded data.
- `cleaner`: clones the `alcor` tree while dropping malformed ALCOR hits, currently invalid `fifo`/`column` combinations.
- `calibrator`: creates or updates the calibrated `time` branch in the decoded `alcor` tree.
- `sorter`: sorts each single-lane stream in calibrated time.
- `after-pulse-suppressor`: suppresses close repeated ALCOR hits per channel.
- `merger`: merges sorted lane/device streams while preserving spill boundaries.
- `trigger`: creates triggered frame trees using declarative trigger configurations.
- `macros/lib/trigger_reader.h`: header-only helper for reading triggered frame output.

Generated ROOT files, logs, build products, and local binaries are ignored by git.
