# ROOT Macro Helper Library

This directory contains reusable headers for ROOT macro analysis.

## trigger_reader.h

`trigger_reader.h` provides a header-only `trigger_reader_t` API for reading triggered frame ROOT files:

```cpp
#include "macros/lib/trigger_reader.h"

trigger_reader_t reader;
if (!reader.open("frames.root"))
    return;

while (reader.next_spill()) {
    while (reader.next_frame()) {
        const auto &trigger   = reader.trigger_hits();
        const auto &timing    = reader.timing_hits();
        const auto &cherenkov = reader.cherenkov_hits();
    }
}
```

The returned vectors contain `hit_t` objects with the original stored fields plus `double time`, the calibrated time value persisted by `trigger.cc`. The vectors remain valid until the next call to `open()`, `next_spill()`, or `next_frame()`.
