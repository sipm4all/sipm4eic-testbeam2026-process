# `display.C`

`display.C` is an interactive Cherenkov event display for the synchronized
triggered-data format. It reads the `frames`, `trigger`, `timing`,
`cherenkov`, and optional ring-result tree through `trigger_reader.h`.

The display draws Cherenkov hits as squares at their stored `(x,y)` positions.
The square colour represents time. Absolute hit times are not used for the
colour scale: the display uses either the minimum frame time or a selected
reference time. The axes are fixed to `x = [-100,100] mm` and
`y = [-100,100] mm`.

## Starting the display

```cpp
.L macros/example/display.C
display("/home/preghenella/CODEX/sorting/aps.rings.multi.time2.root");
```

The square size is configured near the top of `display.C`:

```cpp
constexpr double pixel_size = 3.0;
```

During interactive display:

```text
Enter  show the next selected frame
s      save the current Cherenkov canvas as spill_<spill>_frame_<frame>.png
q      quit
```

The PNG is written in the current working directory.

## References and Filters

The display has one optional reference and any number of filters. References
derive the time subtracted from Cherenkov hits; filters only accept or reject
the current frame. All filters are combined with logical AND.

The polymorphic base classes are:

```cpp
class reference_t {
    virtual bool process(const trigger_reader_t &) = 0;
    virtual double reference_time() const = 0;
};

class selection_t {
    virtual bool is_selected(const trigger_reader_t &) const = 0;
};
```

### No selections

Both of these display every frame:

```cpp
display("ring.root");
display("ring.root", {});
```

The second form is an empty reference pointer. It is unambiguous because the
old overload-based API has been removed.

### Reference objects

Use one of the following reference objects:

```cpp
display("frames.root",
        std::make_shared<field_reference_t>(
            field_selector_t(9, 200, 32, -1, -1)));

display("frames.root",
        std::make_shared<channel_reference_t>(
            channel_selector_t(1, -1)));

display("timing.root",
        std::make_shared<timing_reference_t>("T"));
display("timing.root",
        std::make_shared<timing_reference_t>("T0"));
display("timing.root",
        std::make_shared<timing_reference_t>("T1"));
```

`field_reference_t` uses the first matching hit. `-1` is a wildcard, and the
scan order is trigger, timing, then Cherenkov hits. For type-1 hits,
`channel_selector_t` uses the global ALCOR channel index. For type-9, its
second value is the trigger device.

### Filters

Filters are passed as a vector of `std::shared_ptr<selection_t>`:

```cpp
display("timing.root",
        std::make_shared<timing_reference_t>("T"),
        {
            std::make_shared<trigger_selection_t>(200),
            std::make_shared<ring_selection_t>(
                -100, 100, -100, 100, 10, 80, 1, 2)
        });
```

`trigger_selection_t(200)` requires a type-9 trigger from device 200.
`trigger_selection_t(-1)` accepts a trigger from any device. The default
`trigger_selection_t()` is disabled and always accepts the frame.

## Ring Display and Cuts

If the input contains the optional `ring` tree, each fitted circle or ellipse
is drawn with its fitted time colour. Cuts can be applied to the fitted
centre, radius, and number of rings:

```cpp
ring_selection_t selection(-100., 100.,
                           -100., 100.,
                           10., 80.,
                           1, 2);
display("rings.root",
        std::make_shared<timing_reference_t>("T"),
        {
            std::make_shared<ring_selection_t>(selection),
            std::make_shared<trigger_selection_t>(200)
        });
```

The ring-result tree name is the final optional argument to `display` and
defaults to `"ring"`. This allows a file containing several reconstructions to
be inspected selectively:

```cpp
display("rings.root", nullptr, {},
        std::numeric_limits<int>::min(), 0,
        std::numeric_limits<int>::min(), -1,
        "ring.ellipse");
```

The same selection and reference arguments are used; only the ring tree read
by the reader changes.

The constructor arguments are:

```text
min_x0, max_x0,
min_y0, max_y0,
min_radius, max_radius,
min_n_rings, max_n_rings
```

Ring-associated hits receive an empty circle overlay. The display-side
association uses the constants near the top of `display.C`:

```cpp
constexpr double ring_match_tolerance = 3.5;
constexpr double ring_match_time_window = 2.;
constexpr bool draw_matched_ring_hits = true;
```

The time window is in native timing units (`1 unit = 3.125 ns`). Set
`draw_matched_ring_hits` to `false` to hide the hit overlays without changing
ring drawing or selection.

## Display Details

- The eight PDU outlines are drawn using the configured Cherenkov geometry.
- The palette range is fixed by the constants near the top of the macro.
- The lower canvas shows the Cherenkov time distribution relative to the
  selected reference, when one is used.
- Ring and hit data are read from the current synchronized frame; the old
  flattened `*_frame_start` representation is not required.
