// This Source Code Form is licensed MPL-2.0: http://mozilla.org/MPL/2.0

#include "skfilter.hh"
#include "laddervcf.hh"

#include <fftw3.h>
#include <sys/time.h>
#include <string>
#include <complex>
#include <iostream>

using std::string;
using std::vector;
using std::min;

using namespace SpectMorph; // LadderVCF

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
window_blackman_harris_92 (double x)
{
  if (fabs (x) > 1)
    return 0;

  const double a0 = 0.35875, a1 = 0.48829, a2 = 0.14128, a3 = 0.01168;

  return a0 + a1 * cos (M_PI * x) + a2 * cos (2.0 * M_PI * x) + a3 * cos (3.0 * M_PI * x);
}

void
fft (const uint n_values, float *r_values_in, float *ri_values_out)
{
  auto plan_fft = fftwf_plan_dft_r2c_1d (n_values, r_values_in, (fftwf_complex *) ri_values_out, FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);

  fftwf_execute_dft_r2c (plan_fft, r_values_in, (fftwf_complex *) ri_values_out);

  // usually we should keep the plan, but for this simple test program, the fft
  // is only computed once, so we can destroy the plan here
  fftwf_destroy_plan (plan_fft);
}

double
db (double x)
{
  return 20 * log10 (std::max (x, 0.00000001));
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

vector<float>
filter_partials_test (SKFilter& filter)
{
  const int blocks = 2;
  const int block_size = 8192;
  const int len = block_size * blocks;
  vector<float> in (len), freq (len);

  for (int i = 0; i < len; i++)
    {
      in[i] = ((i % 240) - 120) / 120.0; // 200 Hz Saw
      freq[i] = 5000;
    }
  filter.process_block (len, in.data(), nullptr, freq.data());

  static constexpr int ZERO_PAD = 8;
  vector<float> fft_in (block_size * ZERO_PAD), fft_out (block_size * ZERO_PAD);
  vector<float> fft_out_sum (block_size * ZERO_PAD / 2);

  double window_weight = 0;
  for (int i = 0; i < block_size; i++)
    {
      const double wsize_2 = block_size / 2;
      const double w = window_blackman_harris_92 ((i - wsize_2) / wsize_2);
      window_weight += w;
      fft_in[i] = in[block_size + i] * w;
    }
  for (auto& x : fft_in)
    x *= 2 / window_weight;
  fft (block_size * ZERO_PAD, fft_in.data(), fft_out.data());
  for (size_t i = 0; i < block_size * ZERO_PAD; i += 2)
    fft_out_sum[i / 2] = std::abs (std::complex (fft_out[i], fft_out[i + 1]));

  vector<float> partials;

  for (size_t i = 1; i < fft_out_sum.size() - 1; i++)
    {
      if ((fft_out_sum[i] > fft_out_sum[i - 1]) && (fft_out_sum[i] > fft_out_sum[i + 1]))
        {
          auto norm_freq = 48000 / 2.0 * i / fft_out_sum.size();
          auto norm_height = fft_out_sum[i];
          int partial = (norm_freq + 100) / 200;
          auto fdiff = norm_freq - partial * 200;

          if (abs (fdiff) < 10)
            {
              partials.resize (partial + 1);
              partials[partial] = norm_height;
            }
        }
    }
  return partials;
}

SKFilter::Mode
skmode (const std::string& arg)
{
  vector<string> modes = { "lp1", "lp2", "lp3", "lp4", "lp6", "lp8", "bp2", "bp4", "bp6", "bp8", "hp1", "hp2", "hp3", "hp4", "hp6", "hp8" };
  for (size_t i = 0; i < modes.size(); i++)
    {
      if (arg == modes[i])
        return SKFilter::Mode (i);
    }
  printf ("unsupported filter mode: %s\n", arg.c_str());
  exit (1);
}

struct LadderVCFRef
{
  struct Channel
  {
    struct LP {
      double a = 0;
      double b = 0;

      void
      tick (double in, double g)
      {
        double tmp = 1 / 1.3 * in + a * 0.3 / 1.3;
        b = b + (tmp - b) * g;
        a = in;
      }
    };

    LP lps[4];

    double
    run (double input, double fc, double res)
    {
      fc = M_PI * fc;
      const double g = 0.9892 * fc - 0.4342 * fc * fc + 0.1381 * fc * fc * fc - 0.0202 * fc * fc * fc * fc;
      res *= 1.0029 + 0.0526 * fc - 0.0926 * fc * fc + 0.0218 * fc * fc * fc;

      double g_comp = 0.5;
      lps[0].tick (input - (lps[3].b - g_comp * input) * res * 4, g);
      lps[1].tick (lps[0].b, g);
      lps[2].tick (lps[1].b, g);
      lps[3].tick (lps[2].b, g);
      return lps[3].b;
    }
  } channels[2];

  float freq_ = 440;
  float reso_ = 0;

  void
  set_freq (float freq)
  {
    freq_ = freq;
  }
  void
  set_reso (float reso)
  {
    reso_ = reso;
  }
  void
  run_block (size_t         n_samples,
             float         *left,
             float         *right,
             const float   *freq_in = nullptr,
             const float   *reso_in = nullptr)
  {
    double fc = freq_ / 24000.;
    double res = reso_;
    for (size_t i = 0; i < n_samples; i++)
      {
        left[i]  = channels[0].run (left[i], fc, res);
        right[i] = channels[1].run (right[i], fc, res);
      }
  }
};

template<class Ladder> void
lsweep (Ladder& vcf, double freq, double res)
{
  vector<float> samples = gen();
  vcf.set_freq (freq);
  vcf.set_reso (res);
  for (size_t i = 0; i < samples.size(); i += 2)
    vcf.run_block (1, &samples[i], &samples[i + 1]);
  test (samples);
}

template<class Ladder> void
ladder_perf (const char *label,
           bool          stereo = true,
           const float  *freq_in = nullptr,
           const float  *reso_in = nullptr,
           const float  *drive_in = nullptr)
{
  vector<float> lsamples (48000 * 5);
  for (auto& s : lsamples)
    s = (rand() * (1.0 / RAND_MAX) - 0.5) * 2;
  vector<float> rsamples = lsamples;

  constexpr int    n_runs = 20;
  const double start_t = get_time();
  for (int runs = 0; runs < n_runs; runs++)
    {
      Ladder vcf;
      vcf.set_freq (500);
      vcf.set_reso (0.75);
      size_t i = 0;
      while (i < lsamples.size())
        {
          size_t n_samples = min<size_t> (lsamples.size() - i, 1024);
          vcf.run_block (n_samples, lsamples.data(), stereo ? rsamples.data() : nullptr, freq_in, reso_in, drive_in);
          i += n_samples;
        }
    }
  const double t = (get_time() - start_t) / n_runs;
  printf ("%s: %f ns/sample\n", label, t / lsamples.size() * 1e9);
  printf ("%s: %f stereo streams @ 48kHz\n", label, 1.0 / (t * 48000 / lsamples.size()));
  printf ("\n");
}

void
ladder_perf_streams()
{
  vector<float> freq_in (1024, 440.0);
  vector<float> reso_in (1024, 0.9);
  vector<float> drive_in (1024, 10);

  ladder_perf<LadderVCFNonLinear> ("nl-none");
  ladder_perf<LadderVCFNonLinear> ("nl-f",          true,  &freq_in[0], nullptr);
  ladder_perf<LadderVCFNonLinear> ("nl-f+r",        true,  &freq_in[0], &reso_in[0]);
  ladder_perf<LadderVCFNonLinear> ("nl-f+r+d",      true,  &freq_in[0], &reso_in[0], &drive_in[0]);
  ladder_perf<LadderVCFNonLinear> ("nl-mono-none",  false);
  ladder_perf<LadderVCFNonLinear> ("nl-mono-f",     false, &freq_in[0], nullptr);
  ladder_perf<LadderVCFNonLinear> ("nl-mono-f+r",   false, &freq_in[0], &reso_in[0]);
  ladder_perf<LadderVCFNonLinear> ("nl-mono-f+r+d", false, &freq_in[0], &reso_in[0], &drive_in[0]);
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
      for (int test = 0; test < 6; test++)
        {
          int r, f, mono = test / 3;
          switch (test % 3)
            {
              case 0: r = 0; f = 0; break;
              case 1: r = 0; f = 1; break;
              case 2: r = 1; f = 1; break;
            }

          SKFilter filter (/* oversample */ 4);
          filter.set_mode (skmode (argv[2]));
          filter.set_reso (0.95);

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
            filter.process_block (block_size, left, mono ? nullptr : right, f ? freq : nullptr, r ? reso : nullptr);

          double end_t = get_time();
          double ns_per_sec = 1e9;
          double ns_per_sample = ns_per_sec * (end_t - start_t) / (blocks * block_size);

          printf ("with mono=%d freq=%d reso=%d: ns/sample %f\n", mono, f, r, ns_per_sample);
          printf ("                              bogopolyphony = %f\n\n", 1e9 / (ns_per_sample * 48000));
        }
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

      filter.set_test_linear (true);
      filter.set_mode (skmode (argv[4]));
      filter.set_reso (atof (argv[5]));
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

      filter.set_test_linear (true);
      filter.set_mode (skmode (argv[4]));
      filter.set_reso (1 - R1);

      filter2.set_test_linear (true);
      filter2.set_mode (skmode (argv[4]));
      filter2.set_reso (1 - R2);

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
  if (argc == 5 && cmd == "ldsweep") // <freq> <mode> <reso>
    {
      vector<float> left;
      vector<float> right;
      vector<float> samples = gen();

      for (size_t i = 0; i < samples.size() / 2; i++)
        {
          left.push_back (samples[i * 2]);
          right.push_back (samples[i * 2 + 1]);
        }

      SpectMorph::LadderVCFNonLinear filter;
      filter.set_mode (SpectMorph::LadderVCFMode (atoi (argv[3])));
      filter.set_freq (atof (argv[2]));
      filter.set_reso (atof (argv[4]));

      /* test how the filter behaves as a linear filter (without distortion) */
      float pre_scale = 0.001;
      filter.set_scale (pre_scale, 1 / pre_scale);
      filter.run_block (left.size(), left.data(), right.data());

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
      filter.set_mode (skmode (argv[2]));
      filter.set_reso (last_reso);

      for (int b = 0; b + block_size <= len; b += block_size)
        {
          double reso = ((double) rand() / RAND_MAX);
          reso = 1 - reso * reso;

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
  if (argc == 5 && cmd == "nfilt") // <mode> <res> <vol>
    {
      SKFilter filter (/* oversample */ 4);

      filter.set_mode (skmode (argv[2]));
      filter.set_reso (atof (argv[3]));
      filter.set_drive (atof (argv[4]));

      auto partials = filter_partials_test (filter);
      for (size_t p = 0; p < partials.size(); p++)
        {
          if (partials[p])
            printf ("%zd %f\n", p * 200, db (partials[p]));
        }
    }
  if (argc == 4 && cmd == "nfscan") // <mode> <drive>
    {
      for (int reso = 0; reso < 125; reso++)
        {
          SKFilter filter (/* oversample */ 4);

          filter.set_mode (skmode (argv[2]));
          filter.set_drive (atof (argv[3]));
          filter.set_reso (reso * 0.01);

          auto partials = filter_partials_test (filter);
          printf ("%d %f\n", reso, db (partials[5000 / 200]));
        }
    }
  if (argc == 4 && strcmp (argv[1], "ls-ref") == 0) // <freq> <res>
    {
      LadderVCFRef vcf;
      lsweep (vcf, atof (argv[2]), atof (argv[3]));
    }
  if (argc == 4 && strcmp (argv[1], "ls") == 0) // <freq> <res>
    {
      LadderVCFLinear vcf;
      lsweep (vcf, atof (argv[2]), atof (argv[3]));
    }
  if (argc == 4 && strcmp (argv[1], "ls-nl") == 0) // <freq> <res>
    {
      LadderVCFNonLinear vcf;
      vcf.set_scale (0.01, 100);
      lsweep (vcf, atof (argv[2]), atof (argv[3]));
    }
  if (argc == 2 && strcmp (argv[1], "lperf") == 0)
    {
      ladder_perf<LadderVCFRef>("LadderVCFRef");
      ladder_perf<LadderVCFLinear>("LadderVCFLinear");
      ladder_perf<LadderVCFNonLinear>("MoogVCFNonLinear");
      ladder_perf<LadderVCFNonLinearCheap>("MoogVCFNonLinearCheap");
    }
  if (argc == 2 && strcmp (argv[1], "lperfs") == 0)
    {
      ladder_perf_streams();
    }
}
