#include <cstdio>
#include <cmath>
#include <algorithm>
#include <string>
#include <vector>

using std::string;

#include "flexadsr.hh"

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
}
