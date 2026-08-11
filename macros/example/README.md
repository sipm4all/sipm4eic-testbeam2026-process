# Example ROOT Macros

This directory contains small ROOT macros showing how to inspect processed data.

## trigger_reader.C

Demonstrates the header-only triggered-frame reader API:

```bash
root -l 'macros/example/trigger_reader.C("frames.root")'
```

It iterates over spills and frames and prints the number of trigger, timing, and Cherenkov hits in each frame.

## deltat.C

Reads the triggered-frame `frames` tree and fills 2D delta-t histograms using hits contained in each stored frame. Selector values of `-1` are wildcards.

The original trigger-relative form finds one reference trigger hit and fills `hDeltaT` for all other hits:

```bash
root -l 'macros/example/deltat.C("frames.root", 9, 200, 32, -1, -1)'
```

A second overload compares two selected hit sets directly:

```cpp
deltat(filename,
       target_type, target_device, target_fifo, target_column, target_pixel,
       reference_type, reference_device, reference_fifo, reference_column, reference_pixel,
       outfilename)
```

It fills:

```text
hDeltaT       delta_t = target.time - reference.time vs target channel
hDeltaT_tdc0  same quantity for target hits with tdc == 0
hDeltaT_tdc1  same quantity for target hits with tdc == 1
hDeltaT_tdc2  same quantity for target hits with tdc == 2
hDeltaT_tdc3  same quantity for target hits with tdc == 3
```

For example, compare LASER-side target hits on device 192 FIFOs `16..31` against the test-pulse reference channel `fifo=0,column=0,pixel=0`:

```bash
root -l 'macros/example/deltat.C("triggered.root", 1, 192, -1, -1, -1, 1, 192, 0, 0, 0, "deltat.root")'
```
