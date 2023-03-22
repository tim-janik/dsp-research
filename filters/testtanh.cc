#include <cmath>
#include <cstdio>
#include <vector>
#include <complex>
#include <algorithm>

#include <sys/time.h>

#include <fftw3.h>

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

template<typename F>
void test(F f, int argc, char **argv)
{
  string cmd = argv[2];
  if (cmd == "spect") // <factor>
    {
      float amp = atof (argv[3]);

      const double freq = 1000;
      int BS = 4096;
      vector<float> block;
      for (int i = 0; i < 4096; i++)
        {
          block.push_back (f (sin (i * 2 * M_PI / 48000 * freq) * amp));
        }
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
            out[i] = f (block[i]);
        printf ("flat: %f ns/sample\n", (get_time() - t) / (REPS * block.size()) * 1e9);
      }
      {
        double t = get_time();
        int REPS = 200'000;
        float last = 0;
        for (int j = 0; j < REPS; j++)
          for (size_t i = 0; i < block.size(); i++)
            {
              last = f (block[i] + last);
              out[i] = last;
            }
        printf ("rec:  %f ns/sample\n", (get_time() - t) / (REPS * block.size()) * 1e9);
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
    test (tanh, argc, argv);
  if (cmd == "cheap_tanh")
    test (cheap_tanh, argc, argv);
  if (cmd == "distort")
    test (distort, argc, argv);
#if 0
  for (int i = 0; i < 48000; i++)
    {
      //printf ("%f\n", tanh (sin (i * 2 * M_PI / 48000 * freq) * 2));
    }
#endif
}