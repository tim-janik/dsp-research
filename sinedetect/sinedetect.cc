// This Source Code Form is licensed MPL-2.0: http://mozilla.org/MPL/2.0

#include <vector>
#include <bse/bsemath.hh>
#include <bse/bsemathsignal.hh>
#include <bse/gslfft.hh>

using std::vector;

struct SineDetectPartial
{
  double freq;
  double mag;
};

namespace
{

class QInterpolator
{
  double a, b, c;

public:
  QInterpolator (double y1, double y2, double y3)
  {
    a = (y1 + y3 - 2*y2) / 2;
    b = (y3 - y1) / 2;
    c = y2;
  }
  double
  eval (double x)
  {
    return a * x * x + b * x + c;
  }
  double
  x_max()
  {
    return -b / (2 * a);
  }
};

}

vector<SineDetectPartial>
sine_detect (double mix_freq, const vector<float>& signal)
{
  /* possible improvements for this code
   *
   *  - could produce phases (odd-centric window)
   *  - could eliminate the sines to produce residual signal (spectral subtract)
   */
  vector<SineDetectPartial> partials;

  constexpr double MIN_PADDING = 4;

  size_t padded_length = 2;
  while (signal.size() * MIN_PADDING >= padded_length)
    padded_length *= 2;

  vector<float> padded_signal;
  float window_weight = 0;
  for (size_t i = 0; i < signal.size(); i++)
    {
      const float w = bse_window_cos ((i - signal.size() * 0.5) / (signal.size() * 0.5));
      window_weight += w;
      padded_signal.push_back (signal[i] * w);
    }
  padded_signal.resize (padded_length);

  vector<float> fft_values (padded_signal.size());
  gsl_power2_fftar_simple (padded_signal.size(), padded_signal.data(), fft_values.data());

  vector<float> mag_values;
  for (size_t i = 0; i < fft_values.size(); i += 2)
    mag_values.push_back (sqrt (fft_values[i] * fft_values[i] + fft_values[i + 1] * fft_values[i + 1]));

  for (size_t x = 1; x + 2 < mag_values.size(); x++)
    {
      /* check for peaks
       *  - single peak : magnitude of the middle value is larger than
       *                  the magnitude of the left and right neighbour
       *  - double peak : two values in the spectrum have equal magnitude,
         *                this must be larger than left and right neighbour
       */
      const auto [m1, m2, m3, m4] = std::tie (mag_values[x - 1], mag_values[x], mag_values[x + 1],  mag_values[x + 2]);
      if ((m1 < m2 && m2 > m3) || (m1 < m2 && m2 == m3 && m3 > m4))
        {
          size_t xs, xe;
          for (xs = x - 1; xs > 0 && mag_values[xs] < mag_values[xs + 1]; xs--);
          for (xe = x + 1; xe < (mag_values.size() - 1) && mag_values[xe] > mag_values[xe + 1]; xe++);

          const double normalized_peak_width = double (xe - xs) * signal.size() / padded_length;

          const double mag1 = bse_db_from_factor (mag_values[x - 1], -100);
          const double mag2 = bse_db_from_factor (mag_values[x], -100);
          const double mag3 = bse_db_from_factor (mag_values[x + 1], -100);
          QInterpolator mag_interp (mag1, mag2, mag3);
          double x_max = mag_interp.x_max();
          double peak_mag_db = mag_interp.eval (x_max);
          double peak_mag = bse_db_to_factor (peak_mag_db) * (2 / window_weight);

          if (peak_mag > 0.0001 && normalized_peak_width > 2.9)
            {
              SineDetectPartial partial;
              partial.freq = (x + x_max) * mix_freq / padded_length;
              partial.mag  = peak_mag;
              partials.push_back (partial);
              // printf ("%f %f %f\n", (x + x_max) * mix_freq / padded_length, peak_mag * (2 / window_weight), normalized_peak_width);
            }
        }
    }

  return partials;
}

void
add_sine (double mix_freq, vector<float>& signal, double freq, double mag)
{
  for (size_t i = 0; i < signal.size(); i++)
    signal[i] += sin (i * freq * 2 * PI / mix_freq) * mag;
}

int
main()
{
  double mix_freq = 44100;
  vector<float> test_signal (mix_freq * 0.1); /* 100ms */

  add_sine (mix_freq, test_signal, 513, 1);
  add_sine (mix_freq, test_signal, 8009, 1. / 3);
  for (auto p : sine_detect (mix_freq, test_signal))
    {
      printf ("%f Hz - %f\n", p.freq, p.mag);
    }
}
