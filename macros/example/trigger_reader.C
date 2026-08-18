#include "../lib/trigger_reader.h"

#include <iostream>

void trigger_reader(const char *filename = "frames.root")
{
  trigger_reader_t reader;

  if (!reader.open(filename))
    return;

  while (reader.next_spill()) {
    std::cout << "spill " << reader.spill_id()
              << " nframes=" << reader.nframes()
              << " nsources=" << reader.nsources()
              << std::endl;

    int nprinted = 0;
    for (const auto &source : reader.sources()) {
      std::cout << "  source"
                << " device=" << source.device
                << " fifo=" << source.fifo
                << std::endl;
      if (++nprinted >= 4)
        break;
    }

    while (reader.next_frame()) {
      std::cout << "  frame " << reader.frame_index()
                << " trigger=" << reader.trigger_hits().size()
                << " timing=" << reader.timing_hits().size()
                << " cherenkov=" << reader.cherenkov_hits().size()
                << std::endl;

      for (const auto &hit : reader.cherenkov_hits()) {
        std::cout << "    first cherenkov hit:"
                  << " device=" << hit.device
                  << " fifo=" << hit.fifo
                  << " column=" << hit.column
                  << " pixel=" << hit.pixel
                  << " time=" << hit.time
                  << " x=" << hit.x
                  << " y=" << hit.y
                  << std::endl;
        break;
      }
    }
  }
}
