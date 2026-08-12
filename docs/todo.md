# TODO / Design Notes

This file records design discussions that are not being implemented immediately, but that may be useful for future development decisions.

## Compact decoded ROOT branch types

### Context

The decoded stream tree is currently a simple ROOT `TTree` named `alcor`. It stores one decoded DAQ word per tree entry. This stream-oriented representation is used by the early processing chain:

```text
decoder
  -> checker
  -> cleaner
  -> calibrator
  -> sorter
  -> after-pulse-suppressor
  -> merger
  -> trigger/frame builder
```

The current decoded tree branches are:

```text
device
fifo
type
counter
column
pixel
tdc
rollover
coarse
fine
```

At present all of them are written as ROOT `int` branches, i.e. 32-bit signed integer leaves:

```text
device/I
fifo/I
type/I
counter/I
column/I
pixel/I
tdc/I
rollover/I
coarse/I
fine/I
```

This is robust and simple. It also matches the existing `data_word.h` helper and all current processing programs. However, most fields have a much smaller natural range than 32 bits, so we discussed whether the decoder should already write a more compact representation.

### Measured example

A representative decoded file was inspected:

```text
/home/preghenella/CODEX/sorting/decoded.fifo_13.strict.root
```

ROOT reported approximately:

```text
entries:                 40,913,563
uncompressed tree data:   1.64 GB
compressed ROOT file:    98 MB
compression factor:      about 16.7
```

The observed min/max values in that file were:

```text
field     observed min  observed max
-----     ------------  ------------
device    199           199
fifo      13            13
type      1             15
counter   0             99
column    0             3
pixel     0             3
tdc       0             3
rollover  0             10386
coarse    0             32767
fine      0             104
```

The observed word-type counts were:

```text
type 1  ALCOR hits:    40,913,363
type 7  START_SPILL:          100
type 15 END_SPILL:            100
```

No trigger tags were present in this particular file.

### Current storage cost

With ten 32-bit integer branches, the nominal uncompressed payload is:

```text
10 fields * 4 bytes = 40 bytes per decoded word
```

For 40,913,563 entries this is about:

```text
40,913,563 * 40 bytes = 1,636,542,520 bytes
```

which is consistent with the ROOT uncompressed tree size.

The compressed branch sizes in the inspected file showed that some branches compress extremely well because they are constant or nearly constant in one FIFO file:

```text
device    very small compressed size; constant value
fifo      very small compressed size; constant value
type      very small compressed size; mostly type 1
counter   very small compressed size; mostly 0 on hits, spill counter on markers
rollover  modest compressed size
pixel     moderate compressed size
column    moderate compressed size
tdc       larger compressed size
coarse    larger compressed size
fine      largest compressed size
```

The largest compressed contributors were approximately:

```text
fine      35.7 MB
tdc       22.1 MB
coarse    15.9 MB
column    11.1 MB
pixel      5.7 MB
rollover   3.6 MB
```

The exact compressed size depends on ROOT compression, basket layout, data entropy, and whether the file contains one FIFO or merged data.

### Natural field ranges

The physical/raw ranges are smaller than 32 bits:

```text
field     expected range / meaning
-----     ------------------------
device    board/device IDs such as 192..200; uint16_t is safe
fifo      FIFO ID, currently 0..31 for ALCOR and 32/99 for trigger streams; uint8_t is safe for current IDs
type      word type: 1, 7, 9, 15; uint8_t is safe
counter   spill counter is encoded in 12 bits; trigger counter may use up to 16 bits; uint16_t is safe
column    ALCOR column, 0..7; uint8_t is safe
pixel     ALCOR pixel, 0..3; uint8_t is safe
tdc       TDC index, 0..3; uint8_t is safe
rollover  ALCOR/control rollover can require up to 25 bits; uint32_t is required
coarse    coarse clock is 15 bits, 0..32767; uint16_t is safe
fine      ALCOR fine is 9 bits; decoder also uses fine on spill markers as a status tag; uint16_t is safe
```

A compact but still straightforward branch representation could therefore be:

```text
field     candidate C++ type  bytes
-----     ------------------  -----
device    uint16_t            2
fifo      uint8_t             1
type      uint8_t             1
counter   uint16_t            2
column    uint8_t             1
pixel     uint8_t             1
tdc       uint8_t             1
rollover  uint32_t            4
coarse    uint16_t            2
fine      uint16_t            2
```

Total compact payload:

```text
2 + 1 + 1 + 2 + 1 + 1 + 1 + 4 + 2 + 2 = 17 bytes per decoded word
```

Compared with the current 40 bytes per word, this is:

```text
1 - 17 / 40 = 0.575
```

or about:

```text
57.5% less uncompressed payload
```

For the inspected file, the uncompressed payload would be roughly:

```text
40,913,563 * 17 bytes = 695,530,571 bytes
```

instead of about 1.64 GB, a nominal uncompressed reduction of about 0.94 GB for that file.

### Expected compressed-size gain

The compressed ROOT file is already much smaller than the uncompressed tree:

```text
uncompressed tree data: about 1.64 GB
compressed ROOT file:   about 98 MB
```

Therefore the on-disk saving will not scale linearly with the 57.5% uncompressed payload reduction.

Constant branches such as `device` and `fifo` already compress extremely well even as 32-bit ints, so changing those branch types may save little on disk. Higher-entropy branches such as `fine`, `tdc`, `coarse`, `column`, and `pixel` would likely benefit more.

The practical compressed-size gain is expected to be noticeable but probably much smaller than the raw 57.5% number. It may be on the order of tens of MB per large per-FIFO file, but this should be measured before committing to a format migration.

### Alternative structural savings

The biggest possible structural saving would be to avoid storing fields that are constant in one per-FIFO file, especially:

```text
device
fifo
```

Those could be stored once as file/tree metadata rather than per entry. However, this would be a much deeper data-model change:

- the tree would no longer be self-contained per entry;
- downstream code would need to recover `device` and `fifo` from metadata;
- merged files would need to reintroduce or reconstruct those fields;
- current analysis macros expect these fields in every hit.

For now this is not recommended.

### Interaction with placeholder values

We also discussed placeholder values in the new strict decoder.

The decision for the new decoder is:

```text
Do not use negative placeholders.
```

This keeps the format compatible with a possible future move to unsigned branch types.

Current convention in the strict decoder:

- Control words such as START_SPILL, END_SPILL, and trigger tags use zero for non-applicable hit fields:

```text
column = 0
pixel  = 0
tdc    = 0
```

- ALCOR hits use:

```text
counter = 0
```

because the `counter` field is not meaningful for ALCOR hit payloads.

- START_SPILL and END_SPILL use the `fine` field as a spill-status tag:

```text
fine = 0   normal spill marker
fine = 1   spill payload was suppressed by the decoder
```

If a spill payload is suppressed because its decoding error count exceeds `--allowed-spill-errors`, both the START_SPILL and END_SPILL markers are still written, and both get:

```text
fine = 1
```

This preserves the spill boundary structure while tagging the spill as empty because of decoder errors.

### Why not change the branch types now

The current decision is:

```text
Do not change the main decoded ROOT branch types yet.
```

Reasons:

1. The strict raw decoder is new and should be validated first on problematic FIFOs and full runs.

2. The entire early processing chain currently assumes ROOT `int` branches through `data_word.h`:

```text
decoder
checker
cleaner
calibrator
sorter
after-pulse-suppressor
merger
trigger
macros
```

3. Changing ROOT branch leaf types requires coordinated updates across all readers/writers. A mismatch between branch type and C++ storage type can lead to silent corruption or hard-to-debug ROOT I/O problems.

4. Old decoded ROOT files already exist with `int` branches. If the format changes, the shared reader layer would need to support both legacy and compact files.

5. The actual on-disk saving must be measured. ROOT compression already reduces file size strongly, so the benefit may be smaller than expected.

6. This optimization is independent from the current urgent problem: robust raw decoding in the presence of garbage/stale-memory blocks in the raw stream.

### Possible future implementation plan

If data volume becomes a real bottleneck, the suggested migration path is:

1. Keep the current `int` tree as the stable default until the strict decoder is validated.

2. Add a controlled compact-output option rather than changing the default immediately, for example:

```bash
decoder --input alcdaq.fifo_13.dat --output alcdaq.fifo_13.root --compact
```

or implement a separate converter:

```text
compact.cc
```

that reads the standard decoded `int` tree and writes a compact equivalent.

3. Update `data_word.h` or introduce a more flexible shared reader helper that can detect branch types and read both formats safely.

4. Add tests that compare event counts, spill counts, and field values between the legacy and compact outputs.

5. Measure actual compressed file-size reduction on representative files:

```text
single FIFO with many hits
single FIFO with trigger tags
merged device file
full merged board/run file
```

6. Only after that, decide whether the compact format should become the default.

### Candidate ROOT leaf-list mapping

If implemented, the likely ROOT leaf-list mapping would be something like:

```text
device/s     UShort_t
fifo/b       UChar_t
type/b       UChar_t
counter/s    UShort_t
column/b     UChar_t
pixel/b      UChar_t
tdc/b        UChar_t
rollover/i   UInt_t
coarse/s     UShort_t
fine/s       UShort_t
```

The exact ROOT leaf-list codes must be checked carefully before implementation. This is one of the reasons the migration should be done as a focused task rather than mixed into decoder validation.

