// This Source Code Form is licensed MPL-2.0: http://mozilla.org/MPL/2.0

#include "skfilter.hh"

#include <sys/time.h>
#include <string>
#include <complex>

using std::string;
using std::vector;

static inline float
db_from_complex (float re, float im, float min_dB)
{
  float abs2 = re * re + im * im;

  if (abs2 > 0)
    {
      constexpr float log2_db_factor = 3.01029995663981; // 10 / log2 (10)

      // glibc log2f is a lot faster than glibc log10
      return log2f (abs2) * log2_db_factor;
    }
  else
    return min_dB;
}

inline double
get_time()
{
  /* return timestamp in seconds as double */
  timeval tv;
  gettimeofday (&tv, 0);

  return tv.tv_sec + tv.tv_usec / 1000000.0;
}

vector<float>
gen()
{
  vector<float> samples;
  double phase = 0;
  double l = 48000 * 5;
  double factor = pow (24000 / 20., (1./l));
  double vol = 0;
  for (double f = 20; f < 24000; f *= factor)
    {
      samples.push_back (sin (phase) * vol);
      samples.push_back (cos (phase) * vol);
      phase += f / 48000 * 2 * M_PI;
      vol += 1. / 500; /* avoid click at start */
      if (vol > 1)
        vol = 1;
    }
  return samples;
}

void
test (const vector<float>& samples)
{
  double l = 48000 * 5;
  double factor = pow (24000 / 20., (1./l));
  size_t i = 0;
  for (double f = 20; f < 24000; f *= factor)
    {
      i++;
      if (i > 500 && i * 2 < samples.size())
        {
          printf ("%f %f\n", f, db_from_complex (samples[i * 2], samples[i * 2 + 1], -144));
          // printf ("%f\n", samples[i * 2]);
        }
    }
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
  if (argc == 6 && cmd == "sweep")
    {
      vector<float> left;
      vector<float> right;
      vector<float> freq;
      vector<float> samples = gen();

      const float xfreq = atof (argv[2]);
      for (size_t i = 0; i < samples.size() / 2; i++)
        {
          left.push_back (samples[i * 2]);
          right.push_back (samples[i * 2 + 1]);
          freq.push_back (xfreq);
        }

      SKFilter filter (atoi (argv[3]));

      /* test how the filter behaves as a linear filter (without distortion) */
      float pre_scale = 0.001;
      filter.set_scale (pre_scale, 1 / pre_scale);

      filter.set_params (atoi (argv[4]), atof (argv[5]));
      filter.process_block (left.size(), left.data(), right.data(), freq.data());

      vector<float> out;
      for (size_t i = 0; i < left.size(); i++)
        {
          out.push_back (left[i]);
          out.push_back (right[i]);
        }

      test (out);
    }
  if (argc == 6 && cmd == "4pole")
    {
      // roots of 4th order butterworth in s, left semi-plane
      double sq2inv = 1 / sqrt (2);
      std::complex<double> r1 {-sq2inv,  sq2inv};
      std::complex<double> r2 {-sq2inv, -sq2inv};

      // R must be in interval [0:1]
      std::complex<double> R = 1 - atof (argv[5]);
      std::complex<double> a1 = R + std::sqrt (R * R - 1.0);

      // roots with resonance
      std::complex<double> r1res = r1 * std::sqrt (a1);
      std::complex<double> r2res = r2 * std::sqrt (a1);

      auto R1 = -r1res.real();
      auto R2 = -r2res.real();
      printf ("# HH(s)=1/((s*s+2*s*%f+1)*(s*s+2*s*%f+1))\n", R1, R2);

      vector<float> left;
      vector<float> right;
      vector<float> freq;
      vector<float> samples = gen();

      const float xfreq = atof (argv[2]);
      for (size_t i = 0; i < samples.size() / 2; i++)
        {
          left.push_back (samples[i * 2]);
          right.push_back (samples[i * 2 + 1]);
          freq.push_back (xfreq);
        }

      SKFilter filter (atoi (argv[3]));
      SKFilter filter2 (atoi (argv[3]));

      /* test how the filter behaves as a linear filter (without distortion) */
      float pre_scale = 0.001;
      filter.set_scale (pre_scale, 1 / pre_scale);
      filter2.set_scale (pre_scale, 1 / pre_scale);

      filter.set_params (atoi (argv[4]), 1 - R1);
      filter2.set_params (atoi (argv[4]), 1 - R2);
      filter.process_block (left.size(), left.data(), right.data(), freq.data());
      filter2.process_block (left.size(), left.data(), right.data(), freq.data());

      vector<float> out;
      for (size_t i = 0; i < left.size(); i++)
        {
          out.push_back (left[i]);
          out.push_back (right[i]);
        }

      test (out);
    }
}
