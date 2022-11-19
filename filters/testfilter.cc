// This Source Code Form is licensed MPL-2.0: http://mozilla.org/MPL/2.0

#include "skfilter.hh"

#include <sys/time.h>
#include <string>

using std::string;

inline double
get_time()
{
  /* return timestamp in seconds as double */
  timeval tv;
  gettimeofday (&tv, 0);

  return tv.tv_sec + tv.tv_usec / 1000000.0;
}

int
main (int argc, char **argv)
{
  if (argc < 2)
    {
      printf ("too few args\n");
      return 1;
    }

  string cmd = argv[1];

  if (argc == 2 && cmd == "perf")
    {
      SKFilter filter (/* oversample */ 4);
      filter.set_params (0, 0.95);

      const int block_size = 512;
      float left[block_size], right[block_size], freq[block_size];

      for (int i = 0; i < block_size; i++)
        {
          left[i] = right[i] = ((i % 100) - 50) / 50;
          freq[i] = 440 + i;
        }

      double start_t = get_time();

      const int blocks = 10 * 1000;
      for (int b = 0; b < blocks; b++)
        filter.process_block (block_size, left, right, freq);

      double end_t = get_time();
      double ns_per_sec = 1e9;
      double ns_per_sample = ns_per_sec * (end_t - start_t) / (blocks * block_size);

      printf ("# ns/sample %f\n",  ns_per_sample);
      printf ("# bogopolyphony = %f\n", 1e9 / (ns_per_sample * 48000));
    }

}
