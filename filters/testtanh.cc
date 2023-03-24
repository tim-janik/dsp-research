#include <cmath>
#include <cstdio>
#include <vector>
#include <complex>
#include <algorithm>

#include <sys/time.h>

#include <fftw3.h>

#define PANDA_RESAMPLER_HEADER_ONLY
#include "pandaresampler.hh"

using PandaResampler::Resampler2;


using std::vector;
using std::string;

inline double
window_blackman_harris_92 (double x)
{
  if (fabs (x) > 1)
    return 0;

  const double a0 = 0.35875, a1 = 0.48829, a2 = 0.14128, a3 = 0.01168;

  return a0 + a1 * cos (M_PI * x) + a2 * cos (2.0 * M_PI * x) + a3 * cos (3.0 * M_PI * x);
}

template<typename WFunc>
void
apply_normalized_window (WFunc wfunc, std::vector<float>& fft_in)
{
  double window_weight = 0;
  for (size_t i = 0; i < fft_in.size(); i++)
    {
      const double wsize_2 = fft_in.size() / 2;
      const double w = window_blackman_harris_92 ((i - wsize_2) / wsize_2);
      window_weight += w;
      fft_in[i] *= w;
    }
  for (auto& x : fft_in)
    x *= 2 / window_weight;
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

float
distort (float x)
{
  /* shaped somewhat similar to tanh() and others, but faster */
  x = std::clamp (x, -1.0f, 1.0f);

  return x - x * x * x * (1.0f / 3);
}

// https://www.musicdsp.org/en/latest/Other/238-rational-tanh-approximation.html
float
cheap_tanh (float x)
{
  x = std::clamp (x, -3.0f, 3.0f);

  return (x * (27.0f + x * x) / (27.0f + 9.0f * x * x));
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

enum TFunc
{
  TF_TANH,
  TF_CHEAP_TANH,
  TF_DISTORT,
  TF_CLIP,
  TF_COPY
};

void
upsample (vector<float>& block, int factor)
{
  Resampler2 res (Resampler2::UP, factor, Resampler2::PREC_72DB);
  vector<float> up_block (block.size() * factor);
  res.process_block (block.data(), block.size(), up_block.data());
  block = up_block;
}

void
downsample (vector<float>& block, int factor)
{
  Resampler2 res (Resampler2::DOWN, factor, Resampler2::PREC_72DB);
  vector<float> down_block (block.size() / factor);
  res.process_block (block.data(), block.size(), down_block.data());
  block = down_block;
}

template<TFunc F>
void test(int argc, char **argv)
{
  auto tanh_f = [] (float f) {
    switch (F)
      {
        case TF_TANH:       return std::tanh (f);
        case TF_CHEAP_TANH: return cheap_tanh (f);
        case TF_DISTORT:    return distort (f);
        case TF_CLIP:       return std::clamp (f, -1.0f, 1.0f);
        case TF_COPY:       return f;
      }
  };
  string cmd = argv[2];
  if (cmd == "spect") // <factor>
    {
      float amp = atof (argv[3]);

      const double freq = 876;
      int BS = 4096;
      vector<float> block;
      for (int i = 0; i < 4096; i++)
        {
          block.push_back (sin (i * 2 * M_PI / 48000 * freq));
        }
      upsample (block, 8);
      for (size_t i = 0; i < block.size(); i++)
        block[i] = tanh_f (block[i] * amp);
      downsample (block, 8);
      apply_normalized_window (window_blackman_harris_92, block);
      vector<float> out (BS + 2);
      fft (block.size(), block.data(), out.data());

      for (size_t i = 0; i < out.size(); i += 2)
        {
          printf ("%f\n", db (std::abs (std::complex (out[i], out[i + 1]))));
        }
    }
  if (cmd == "perf")
    {
      vector<float> block;
      vector<float> out;

      for (int i = 0; i < 512; i++)
        block.push_back ((rand() * (1.0 / RAND_MAX) - 0.5) * 2);

      out.resize (block.size());

      {
        double t = get_time();
        int REPS = 200'000;
        for (int j = 0; j < REPS; j++)
          for (size_t i = 0; i < block.size(); i++)
            out[i] = tanh_f (block[i]);
        printf ("flat: %f ns/sample\n", (get_time() - t) / (REPS * block.size()) * 1e9);
      }
      {
        double t = get_time();
        int REPS = 200'000;
        float last = 0;
        for (int j = 0; j < REPS; j++)
          for (size_t i = 0; i < block.size(); i++)
            {
              last = tanh_f (block[i] + last);
              out[i] = last;
            }
        printf ("rec:  %f ns/sample\n", (get_time() - t) / (REPS * block.size()) * 1e9);
      }
    }
  if (cmd == "plot")
    {
      for (double x = -5; x <= 5; x += 0.01)
        printf ("%f %f %f\n", x, tanh (x), tanh_f (x));
    }
  if (cmd == "tplot")
    {
      float amp = atof (argv[3]);
      const double freq = 876;
      for (int i = 0; i < 48000; i++)
        {
          printf ("%f %f\n", tanh (sin (i * 2 * M_PI / 48000 * freq) * amp), tanh_f (sin (i * 2 * M_PI / 48000 * freq) * amp));
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

    if (cmd == "tanh")
      test<TF_TANH> (argc, argv);
    if (cmd == "cheap_tanh")
      test<TF_CHEAP_TANH> (argc, argv);
    if (cmd == "distort")
      test<TF_DISTORT> (argc, argv);
    if (cmd == "clip")
      test<TF_CLIP> (argc, argv);
    if (cmd == "copy")
      test<TF_COPY> (argc, argv);
  }
