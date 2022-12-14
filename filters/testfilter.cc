// This Source Code Form is licensed MPL-2.0: http://mozilla.org/MPL/2.0

#include "skfilter.hh"
#include "laddervcf.hh"

#include <sys/time.h>
#include <string>
#include <complex>
#include <iostream>

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

  if (argc == 3 && cmd == "perf")
    {
      SKFilter filter (/* oversample */ 4);
      filter.set_params (atoi (argv[2]), 0.95);

      const int block_size = 512;
      float left[block_size], right[block_size], freq[block_size], reso[block_size];

      for (int i = 0; i < block_size; i++)
        {
          left[i] = right[i] = ((i % 100) - 50) / 50;
          freq[i] = 440 + i;
          reso[i] = i / double (block_size);
        }

      double start_t = get_time();

      const int blocks = 10 * 1000;
      for (int b = 0; b < blocks; b++)
        filter.process_block (block_size, left, right, freq, reso);

      double end_t = get_time();
      double ns_per_sec = 1e9;
      double ns_per_sample = ns_per_sec * (end_t - start_t) / (blocks * block_size);

      printf ("# ns/sample %f\n",  ns_per_sample);
      printf ("# bogopolyphony = %f\n", 1e9 / (ns_per_sample * 48000));
    }
  if (argc == 6 && cmd == "sweep") // <freq> <over> <mode> <res>
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
  if (argc == 6 && cmd == "4pole") // <freq> <over> <mode> <res>
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
      printf ("R %f %f %f\n", atof (argv[5]), R1, R2);
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
  if (argc == 5 && cmd == "ldsweep")
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

      SpectMorph::LadderVCFNonLinear filter;
      filter.set_mode (SpectMorph::LadderVCFMode (atoi (argv[3])));

      /* test how the filter behaves as a linear filter (without distortion) */
      float pre_scale = 0.001;
      filter.set_scale (pre_scale, 1 / pre_scale);

      const float *inputs[2] = { left.data(), right.data() };
      float *outputs[2] = { left.data(), right.data() };
      filter.run_block (left.size(), 0, atof (argv[4]), inputs, outputs, true, true, freq.data(), nullptr);

      vector<float> out;
      for (size_t i = 0; i < left.size(); i++)
        {
          out.push_back (left[i]);
          out.push_back (right[i]);
        }

      test (out);
    }
  if (argc == 4 && cmd == "roots") // <order> <res>
    {
      int order = atoi (argv[2]);

      // R must be in interval [0:1]
      const double R = 1 - atof (argv[3]);
      const double r_alpha = std::acos (R) / (order / 2);

      vector<double> Rn;
      for (int i = 0; i < order / 2; i++)
        {
          /* butterworth roots in s, left semi plane */
          const double bw_s_alpha = M_PI * (4 * i + order + 2) / (2 * order);
          /* add resonance */
          Rn.push_back (-cos (bw_s_alpha + r_alpha));
        }

      std::sort (Rn.begin(), Rn.end(), std::greater<double>());

      printf ("H%d(s)=1/(", order);
      for (size_t i = 0; i < Rn.size(); i++)
        {
          if (i)
            printf ("*");
          printf ("(s*s+2*s*%.17g+1)", Rn[i]);
        }
      printf (")\n");
      for (size_t i = 0; i < Rn.size(); i++)
        {
          printf ("H%d_%zd(s)=1/(s*s+2*s*%.17g+1)\n", order, i, Rn[i]);
        }
    }
  if (argc == 3 && cmd == "rsweep") // <mode>
    {
      SKFilter filter (/* oversample */ 4);

      const int block_size = 8192;
      const int len = block_size * 30;
      float left[len], right[len], freq[len];

      for (int i = 0; i < len; i++)
        {
          const auto saw_len = 500;
          const auto saw_len1 = 505;
          const auto sl2 = saw_len * 0.5;
          const auto sl21 = saw_len1 * 0.5;
          left[i] = right[i] = ((i % saw_len) - sl2) / sl2 + ((i % saw_len1) - sl21) / sl21;
          freq[i] = 1000;
        }

      double last_reso = 0;
      filter.set_params (atoi (argv[2]), last_reso);

      for (int b = 0; b + block_size <= len; b += block_size)
        {
          double reso = ((float)rand() / RAND_MAX);
          reso = 1 - reso * reso;
          filter.set_scale (0.05, 1);

          float reso_in[block_size];
          for (int i = 0; i < block_size; i++)
            reso_in[i] = last_reso + (reso - last_reso) * i / double (block_size);

          last_reso = reso;
          filter.process_block (block_size, left + b, right + b, freq + b, reso_in);
        }
      float mx = 0;
      for (int i = 0; i < len; i++)
        mx = std::max (std::abs (left[i]), mx);
      for (int i = 0; i < len; i++)
        printf ("%f %f\n", left[i] / mx, right[i] / mx);
    }
  if (argc == 3 && cmd == "lmax") // <res>
    {
      /* numerically find maximum of HLP (freq, res)
       * there is probably a direct solution to this problem as well
       */
      auto HLP = [] (double freq, double res) {
        auto s = std::complex (0.0, freq);

        return std::abs (1.0/(s*s + 2*(1-res)*s + 1.0));
      };
      double res = atof (argv[2]);
      double x = 0.1, dx = 0.1;
      double best = -1;
      while (dx > 1e-14)
        {
          double v;
          v = HLP (x + dx, res);
          if (v > best)
            {
              x += dx;
              best = v;
            }
          else
            {
              v = HLP (x - dx, res);
              if (v > best)
                {
                  x -= dx;
                  best = v;
                }
              else
                {
                  dx *= 0.5;
                }
            }
        }
      printf ("%f\n", best);
    }
}
