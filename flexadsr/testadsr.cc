#include <cstdio>
#include <cmath>
#include <algorithm>
#include <string>
#include <vector>

#include <sys/time.h>

using std::string;

#include "flexadsr.hh"
#include "envelope.hh"

inline double
get_time()
{
  /* return timestamp in seconds as double */
  timeval tv;
  gettimeofday (&tv, 0);

  return tv.tv_sec + tv.tv_usec / 1000000.0;
}

float no_optimize = 0;

int
main (int argc, char **argv)
{
  ADSR env;
  auto set_shape = [&env] (string s) {
    if (s == "LIN")
      env.set_shape (ADSR::Shape::LINEAR);
    else if (s == "EXP")
      env.set_shape (ADSR::Shape::EXPONENTIAL);
    else if (s == "FLEX")
      env.set_shape (ADSR::Shape::FLEXIBLE);
    else
      {
        fprintf (stderr, "bad shape\n");
        exit (1);
      }
  };

  int t = 0;
  auto process = [&] (size_t n) {
    std::vector<float> samples (n);
    env.process (&samples[0], samples.size());
    for (auto s : samples)
      printf ("%f %f\n", t++ / 48000., s);
  };

  if (string (argv[1]) == "plot")
    {
      env.set_attack (atof (argv[2]));
      env.set_attack_slope (atof (argv[3]));
      env.set_decay (atof (argv[4]));
      env.set_decay_slope (atof (argv[5]));
      env.set_sustain (atof (argv[6]) * 100);
      env.set_release (atof (argv[7]));
      env.set_release_slope (atof (argv[8]));
      string s = argv[10];

      if (s == "LIN")
        env.set_shape (ADSR::Shape::LINEAR);
      else if (s == "EXP")
        env.set_shape (ADSR::Shape::EXPONENTIAL);
      else if (s == "FLEX")
        env.set_shape (ADSR::Shape::FLEXIBLE);

      env.set_rate (48000);
      env.start();

      process (std::clamp<double> (atof (argv[9]), 0, 3) * 48000);

      env.stop();

      process (48000);
    }
  if (string (argv[1]) == "achange")
    {
      set_shape (argv[2]);

      env.set_attack (1);
      env.set_attack_slope (-0.5);
      env.set_rate (48000);
      env.start();

      process (12000);

      env.set_attack (0.5);

      process (48000);
    }
  if (string (argv[1]) == "dchange")
    {
      set_shape (argv[2]);
      env.set_attack (0.25);
      env.set_attack_slope (1);
      env.set_rate (48000);
      env.set_decay (3);
      env.set_decay_slope (1);
      env.start();

      process (48000);

      env.set_decay (0.1);

      process (48000);
    }
  if (string (argv[1]) == "rchange")
    {
      set_shape (argv[2]);
      env.set_attack (0.1);
      env.set_decay (0.1);
      env.set_sustain (75);
      env.set_release (2);
      env.set_release_slope (-0.5);
      env.start();
      process (48000);
      env.stop();
      process (48000);
      env.set_release (3);
      process (96000);
    }
  if (string (argv[1]) == "perf") // <shape>
    {
      env.set_attack (0.3);
      env.set_decay (0.3);
      set_shape (argv[2]);
      const int BS = 64;
      const int BLOCKS = 50;
      const int REPS = 10'000;
      double start_time = get_time();
      for (int r = 0; r < REPS; r++)
        {
          env.start();
          for (int i = 0; i < BLOCKS; i++)
            {
              float freq_in[BS];
              env.process (freq_in, BS);
              for (int i = 0; i < BS; i++)
                {
                  freq_in[i] = exp2f (freq_in[i]);
                  no_optimize += freq_in[i];
                }
            }
        }
      double t = get_time() - start_time;
      printf ("%f ns/sample\n", t * 1e9 / (REPS * BLOCKS * BS));
      printf ("%f streams @ 48kHz\n", 1.0 / (t * 48000 / (REPS * BLOCKS * BS)));
    }
  if (string (argv[1]) == "lperf") // <shape>
    {
      Envelope envelope;

      envelope.set_attack (0.3);
      envelope.set_decay (0.3);
      string s = argv[2];
      if (s == "LIN")
        envelope.set_shape (Envelope::Shape::LINEAR);
      else if (s == "EXP")
        envelope.set_shape (Envelope::Shape::EXPONENTIAL);
      else
        {
          fprintf (stderr, "bad shape\n");
          exit (1);
        }
      const int BS = 64;
      const int BLOCKS = 50;
      const int REPS = 10'000;
      double start_time = get_time();
      for (int r = 0; r < REPS; r++)
        {
          envelope.start (48000);
          for (int i = 0; i < BLOCKS; i++)
            {
              float freq_in[BS];
              for (int i = 0; i < BS; i++)
                {
                  freq_in[i] = exp2f (envelope.get_next());
                  no_optimize += freq_in[i];
                }
            }
        }
      double t = get_time() - start_time;
      printf ("%f ns/sample\n", t * 1e9 / (REPS * BLOCKS * BS));
      printf ("%f streams @ 48kHz\n", 1.0 / (t * 48000 / (REPS * BLOCKS * BS)));
    }
}
