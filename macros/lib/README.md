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

The optional second argument to `open` selects the ring-result tree and
defaults to `"ring"`:

```cpp
if (!reader.open("frames.root", "ring.ellipse"))
    return;
```

This is useful when several ring-finder configurations are stored in the same
ROOT file. The selected tree must have the standard ring branches.

The returned vectors contain `hit_t` objects with the original stored fields plus `double time`, `double x`, and `double y`, the calibrated time and spatial coordinates persisted by `trigger.cc`. The vectors remain valid until the next call to `open()`, `next_spill()`, or `next_frame()`.


## Spill participation

Triggered ROOT files may contain a `spill_participation` tree propagated from the merger. `trigger_reader_t` reads it automatically when present. For the current spill:

```cpp
std::cout << reader.nsources() << std::endl;

for (const auto &source : reader.sources()) {
    std::cout << source.device << " " << source.fifo << std::endl;
}
```

The returned source vector lists only the `(device,fifo)` sources that contributed to the current spill. Files without the metadata tree simply return an empty source vector.

## frame_selection.h

`frame_selection.h` contains the shared polymorphic selection API used by
`display.C` and `deltat.C`.

The API separates three roles:

```text
target_t     selects individual hits
reference_t  calculates one reference time per frame
selection_t  accepts or rejects a complete frame
```

Targets are processed once per frame and combined with logical AND. This is
used, for example, to select type-1 hits that also belong to a stored ring:

```cpp
std::vector<target_ptr_t> targets = {
    std::make_shared<channel_target_t>(-1),
    std::make_shared<ring_target_t>(8, 64, 3.5, 1., 200.)
};
```

The available hit-target constructors are:

```cpp
field_target_t(device, fifo, column, pixel)
fifo_target_t(device, fifo)
channel_target_t(channel)
trigger_target_t(device)
```

The first three select type-1 hits. `trigger_target_t` selects type-9 hits;
`-1` is a wildcard for address fields where applicable.

`ring_target_t` consumes the rings already stored in the input `ring` tree. It
does not run RANSAC. `ring_selection_t` is a separate frame-level filter used
to select frames by stored ring centre, radius, and ring count.
