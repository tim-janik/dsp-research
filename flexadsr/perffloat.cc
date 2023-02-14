#include <string>

struct ADSR_F
{
  float a_ = 0;
  float b_ = 0;
  float c_ = 0;
  float rate_ = 44100;

  void
  init_abc (float time_s, float slope)
  {
    bool positive = slope > 0;
    slope = std::abs (slope);

    const float t1y = 0.5f + 0.25f * slope;

    a_ = slope * ( 1.0135809670870777f + slope * (-1.2970447050283254f + slope *   7.2390617313972063f));
    b_ = slope * (-5.8998946320566281f + slope * ( 5.7282487210570903f + slope * -15.525953208626062f));
    c_ = 1 - (t1y * a_ + b_) * t1y;

    if (!positive)
      {
        c_ += a_ + b_;
        b_ = -2 * a_ - b_;
      }

    const float time_factor = 1 / (rate_ * time_s);
    a_ *= time_factor;
    b_ *= time_factor;
    c_ *= time_factor;
  }
};

struct ADSR_D
{
  float a_ = 0;
  float b_ = 0;
  float c_ = 0;
  float rate_ = 44100;

  void
  init_abc (float time_s, float slope)
  {
    bool positive = slope > 0;
    slope = std::abs (slope);

    const float t1y = 0.5 + 0.25 * slope;

    a_ = slope * ( 1.0135809670870777 + slope * (-1.2970447050283254 + slope *   7.2390617313972063));
    b_ = slope * (-5.8998946320566281 + slope * ( 5.7282487210570903 + slope * -15.525953208626062));
    c_ = 1 - (t1y * a_ + b_) * t1y;

    if (!positive)
      {
        c_ += a_ + b_;
        b_ = -2 * a_ - b_;
      }

    const float time_factor = 1 / (rate_ * time_s);
    a_ *= time_factor;
    b_ *= time_factor;
    c_ *= time_factor;
  }
};

float dont_optimize = 0;

int
main (int argc, char **argv)
{
  std::string m = argv[1];

  int runs = 100'000'000;
  if (m == "float")
    {
      ADSR_F f;
      for (int i = 0; i < runs; i++)
        {
          float x = float (i) / runs;
          f.init_abc (x, x);
          dont_optimize += f.a_ + f.b_ + f.c_;
        }
    }
  else if (m == "double")
    {
      ADSR_D d;
      for (int i = 0; i < runs; i++)
        {
          float x = float (i) / runs;
          d.init_abc (x, x);
          dont_optimize += d.a_ + d.b_ + d.c_;
        }
    }
}
