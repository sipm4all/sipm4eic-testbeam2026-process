# Example ROOT Macros

This directory contains small ROOT macros showing how to inspect processed data.

## trigger_reader.C

Demonstrates the header-only triggered-frame reader API:

```bash
root -l 'macros/example/trigger_reader.C("frames.root")'
```

It iterates over spills and frames and prints the number of trigger, timing, and Cherenkov hits in each frame.

## deltat.C

Reads the triggered-frame `frames` tree and fills 2D delta-t histograms using hits contained in each stored frame. Selector values of `-1` are wildcards. All delta-t histograms use 2048 bins over `[-32, 32]`.

The original trigger-relative form finds one reference trigger hit and fills `hDeltaT` for all other hits:

```bash
root -l 'macros/example/deltat.C("frames.root", 9, 200, 32, -1, -1)'
```

The macro also defines selector types with fixed constructor arity:

```cpp
field_selector_t(type, device, fifo, column, pixel)
channel_selector_t(type, channel_or_trigger_device)
```

The exact constructor arity is intentional: missing selector fields are not silently zero-initialised.

Field selectors work for all stored hit categories, including trigger tags:

```cpp
deltat("frames.root",
       field_selector_t(1, 192, -1, -1, -1),
       field_selector_t(1, 192, 0, 0, 0),
       "deltat.root");
```

Channel selectors use the same global channel index as the `hDeltaT` x axis for ALCOR hits:

```cpp
channel_selector_t(1, channel)
```

where:

```cpp
channel = pixel + 4 * column + 32 * (fifo / 4)
        + 256 * (device - 192);
```

For trigger tags, the compact convention is:

```cpp
channel_selector_t(9, device)
```

which is interpreted as:

```cpp
field_selector_t(9, device, 32, -1, -1)
```

The target/reference overloads compute:

```text
delta_t = target.time - reference.time
```

and fill:

```text
hDeltaT       delta_t vs target channel
hDeltaT_spill delta_t vs spill id, 100 bins over [0,100]
hDeltaT_tdc0  delta_t vs target fine for target hits with tdc == 0
hDeltaT_tdc1  delta_t vs target fine for target hits with tdc == 1
hDeltaT_tdc2  delta_t vs target fine for target hits with tdc == 2
hDeltaT_tdc3  delta_t vs target fine for target hits with tdc == 3
```

For example, compare all ALCOR target channels against trigger tags from device 200:

```cpp
deltat("triggered.root",
       channel_selector_t(1, -1),
       channel_selector_t(9, 200),
       "deltat.root");
```

The old positional target/reference form remains available as a wrapper:

```cpp
deltat(filename,
       target_type, target_device, target_fifo, target_column, target_pixel,
       reference_type, reference_device, reference_fifo, reference_column, reference_pixel,
       outfilename);
```
