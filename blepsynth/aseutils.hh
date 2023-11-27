// This Source Code Form is licensed MPL-2.0: http://mozilla.org/MPL/2.0
#pragma once

#include <fftw3.h>

#include <stdarg.h>
#include <cassert>
#include <cmath>
#include <complex>
#include <string>
#include <vector>

// detect compiler
#if __clang__
  #define ASEUTILS_COMP_CLANG 1
#elif __GNUC__ > 2
  #define ASEUTILS_COMP_GCC 1
#else
  #error "unsupported compiler"
#endif

#if ASEUTILS_COMP_GCC
  #define ASEUTILS_PRINTF(format_idx, arg_idx)      __attribute__ ((__format__ (gnu_printf, format_idx, arg_idx)))
#else
  #define ASEUTILS_PRINTF(format_idx, arg_idx)      __attribute__ ((__format__ (__printf__, format_idx, arg_idx)))
#endif

static inline void
gsl_power2_fftac (unsigned int n_values, double *ri_in, double *ri_out)
{
  assert (false);
}

static inline float
random_frange (float min, float max)
{
  return min + static_cast<float> (std::rand()) / RAND_MAX * (max - min);
}

static inline void
floatfill (float *ptr, float value, unsigned int n_values)
{
  std::fill_n (ptr, n_values, value);
}

static inline float
fast_exp2 (float v)
{
  return exp2f (v);
}

static inline int
irintf (float f)
{
  return __builtin_rintf (f);
}

static inline float
fast_voltage2hz (float f)
{
  // FIXME: remove this function after implementing freq modulation
  return f;
}

#define ASE_SIGNAL_FROM_FREQ(f) f /* FIXME: remove me */

static inline double
ase_db_from_factor (double factor, double min_dB)
{
  if (factor > 0)
    {
      double dB = log10 (factor); /* Bell */
      dB *= 20;
      return dB;
    }
  else
    return min_dB;
}

static inline double
ase_window_cos (double x) /* von Hann window */
{
  if (fabs (x) > 1)
    return 0;
  return 0.5 * cos (x * M_PI) + 0.5;
}

static inline std::string string_printf (const char *format, ...) ASEUTILS_PRINTF (1, 2);

static inline std::string
string_printf (const char *format, ...)
{
  std::vector<char> buffer;
  va_list ap;

  /* figure out required size */
  va_start (ap, format);
  int size = vsnprintf (nullptr, 0, format, ap);
  va_end (ap);

  if (size < 0)
    return format;

  /* now print with large enough buffer */
  buffer.resize (size + 1);

  va_start (ap, format);
  size = vsnprintf (&buffer[0], buffer.size(), format, ap);
  va_end (ap);

  if (size < 0)
    return format;

  return &buffer[0];
}

#define ASE_ASSERT_RETURN(expr,...) do { if (expr) break; assert (false); } while (0)

static inline void
fft (unsigned int n_values, const double *in, double *out)
{
  assert ((n_values & (n_values - 1)) == 0); // power-of-two check
  auto plan_fft = fftw_plan_dft_1d (n_values, (fftw_complex *) in, (fftw_complex *) out, FFTW_FORWARD, FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);

  fftw_execute_dft (plan_fft, (fftw_complex *) in, (fftw_complex *) out);

  // usually we should keep the plan, but for this simple test program, the fft
  // is only computed once, so we can destroy the plan here
  fftw_destroy_plan (plan_fft);
}

static inline void
ifft (unsigned int n_values, double *in, double *out)
{
  assert ((n_values & (n_values - 1)) == 0); // power-of-two check
  auto plan_fft = fftw_plan_dft_1d (n_values, (fftw_complex *) in, (fftw_complex *) out, FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);

  fftw_execute_dft (plan_fft, (fftw_complex *) in, (fftw_complex *) out);

  // usually we should keep the plan, but for this simple test program, the fft
  // is only computed once, so we can destroy the plan here
  fftw_destroy_plan (plan_fft);
}

static inline std::vector<std::complex<double>>
fft (const std::vector<std::complex<double>>& in)
{
  std::vector<std::complex<double>> out (in.size());
  fft (in.size(), reinterpret_cast<const double *> (in.data()), reinterpret_cast<double *> (out.data()));
  return out;
}
