#include <sys/time.h>
#include "saturationdsp.hh"
#include <string>

using namespace Ase;

using std::vector;
using std::string;

using PandaResampler::Resampler2;

inline double
get_time()
{
  /* return timestamp in seconds as double */
  timeval tv;
  gettimeofday (&tv, 0);

  return tv.tv_sec + tv.tv_usec / 1000000.0;
}

float global = 0;

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
      SaturationDSP dist;

      for (int mode = 0; mode <= 3; mode++)
        {
          for (int const_params = 0; const_params < 2; const_params++)
            {
              dist.reset (48000);
              dist.set_mode (static_cast<SaturationDSP::Mode> (mode));
              dist.set_mix (50, false);

              const int block_size = 512;
              float left[block_size], right[block_size];

              for (int i = 0; i < block_size; i++)
                {
                  left[i] = right[i] = ((i % 100) - 50) / 50;
                }
              double start_t = get_time();
              const int blocks = 50 * 1000;
              for (int b = 0; b < blocks; b++)
                {
                  dist.process<true> (left, right, left, right, block_size);
                  if (!const_params)
                    {
                      // ensure that smoothing needs to be used all the time
                      if (b & 1)
                        dist.set_drive (0, false);
                      else
                        dist.set_drive (36, false);
                    }
                }
              double end_t = get_time();
              double ns_per_sec = 1e9;
              double ns_per_sample = ns_per_sec * (end_t - start_t) / (blocks * block_size);
              printf ("with mode=%d, const=%d: ns/sample %f\n", mode, const_params, ns_per_sample);
              printf ("                    bogopolyphony = %f\n\n", 1e9 / (ns_per_sample * 48000));
              // no optimization
              for (int i = 0; i < block_size; i++)
                global += left[i] + right[i];
            }
        }
    }
  if (argc == 2 && cmd == "sin")
    {
      SaturationDSP dist;
      for (float f = -5; f < 5; f += 0.001)
        {
          //printf ("%f %f %f\n", f, dist.lookup_table (f), tanh_restricted (f));
        }
      int SR = 48000;
      vector<float> input (64), output (64);
      int x = 0;
      while (x < SR * 5)
        {
          for (size_t i = 0; i < input.size(); i++)
            input[i] = sin (x++ * 4400. * 2 * M_PI / SR);
          dist.set_drive (48.0 * x / SR / 5, true);
          dist.process<false> (input.data(), nullptr, output.data(), nullptr, 64);
          for (auto f : output)
            printf ("%f\n", f * 0.5);
        }
    }
}
