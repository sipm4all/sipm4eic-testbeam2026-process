# Trigger Configuration Language

Trigger configurations are consumed by `trigger`:

```bash
build/trigger \
  --input merged.root \
  --output frames.root \
  --config process/config/trigger/trigger_range.conf \
  --window 256
```

A file contains named selectors and one trigger block.

## Selectors

```text
selector TIMING {
    type        = 1
    device      = 200
    fifo        = [0,7]
    column      = *
    pixel       = *
    time_offset = 0.
}
```

Supported integer fields include:

```text
type device fifo column pixel counter tdc
```

Supported matching syntax:

```text
*               wildcard
3               exact match
[0,7]           inclusive range
{0,2,4,6}       explicit set
```

`time_offset` is a selector-level corrected-time offset:

```text
corrected_time = raw_time - selector.time_offset
```

When using files produced by `calibrator`, these offsets should normally be `0.` because low-level timing alignment has already been applied.

## Seeded Trigger Mode

```text
trigger {
    seed TRIGGER

    require TIMING {
        dt  = [-10,+10]
        min = 32
    }
}
```

The seed establishes the event time. `require` blocks are ANDed. `veto` blocks require zero matching hits in their `dt` interval.

## Self-Coincidence Trigger Mode

```text
trigger {
    coincidence TIMING {
        dt  = [-10,+10]
        min = 32
    }
}
```

This produces one frame per coincidence cluster, not one frame per hit. By default `min` counts hits.

Distinct-channel multiplicity is supported:

```text
trigger {
    coincidence TIMING {
        dt      = [-10,+10]
        min     = 32
        unique  = channel
    }
}
```

Here `channel` means the tuple `(device,fifo,column,pixel)`.

## Frame Center

If no `frame_center` is specified, frames are centered on event time. A trigger can instead center on the expected timing of a selector:

```text
trigger {
    seed TRIGGER
    frame_center DETECTOR
}
```

This uses `DETECTOR.time_offset` as a timing calibration reference; it does not search for a particular detector hit.

The examples in this directory demonstrate single-channel, coincidence, multiplicity, veto, range, set, and self-coincidence triggers.
