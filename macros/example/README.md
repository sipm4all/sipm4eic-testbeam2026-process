# Example ROOT Macros

This directory contains small ROOT macros showing how to inspect processed data.

## trigger_reader.C

Demonstrates the header-only triggered-frame reader API:

```bash
root -l 'macros/example/trigger_reader.C("frames.root")'
```

It iterates over spills and frames and prints the number of trigger, timing, and Cherenkov hits in each frame.

## deltat.C

Reads the triggered-frame `frames` tree and fills a 2D delta-t histogram using all hits contained in each stored frame. The ROOT function name is `deltat`, matching the filename. Example:

```bash
root -l 'macros/example/deltat.C("frames.root", 9, 200, 32, -1, -1)'
```
