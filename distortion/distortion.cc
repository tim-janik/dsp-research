#include <algorithm>
#include <sys/time.h>

#include "distortiondsp.hh"

#include <vector>
#include <complex>
#include <cmath>
#include <iostream>

using std::complex;

static constexpr double PI = 3.14159265358979323846;

class LatencyTester
{
  DistortionDSP *distortion_dsp;
  double sample_rate_hz;
  double freq_hz;
  double phase;
  int warmup_samples;
public:
  LatencyTester(double sample_rate_hz, double freq_hz, DistortionDSP *distortion_dsp)
    : distortion_dsp (distortion_dsp),
      sample_rate_hz(sample_rate_hz),
      freq_hz (freq_hz),
      phase(0.0),
      warmup_samples((int)(0.1 * sample_rate_hz)) // 100 ms
  {
  }

  double run_test (int total_samples, int block_size)
  {
    std::vector<float> left (block_size);
    std::vector<float> right (block_size);

    double phase_step = 2.0 * PI * freq_hz / sample_rate_hz;

    int sample_index = 0;

    double re = 0.0;
    double im = 0.0;
    int measured_samples = 0;

    while (sample_index < total_samples)
      {
        // ----------------------------
        // Generate input sine
        // ----------------------------
        for (int i = 0; i < block_size; i++)
          {
            double sample = std::sin(phase);
            phase += phase_step;

            left[i] = (float)sample;
            right[i] = (float)sample;
          }

        // ----------------------------
        // Process through Distortion
        // ----------------------------
        distortion_dsp->process_block (left.data(), right.data(), block_size);

        // ----------------------------
        // Warm-up discard
        // ----------------------------
        if (sample_index < warmup_samples)
          {
            sample_index += block_size;
            continue;
          }

        // ----------------------------
        // Lock-in phase measurement
        // ----------------------------
        for (int i = 0; i < block_size; i++)
          {
            double t = (double)(sample_index + i) / sample_rate_hz;

            double ref_cos = std::cos(2.0 * PI * freq_hz * t);
            double ref_sin = std::sin(2.0 * PI * freq_hz * t);

            double y = left[i];

            // multiply signal with complex exp function from fourier transform:
            //
            //   v * exp (-j * x) = v * (cos (x) - j * sin (x))
            re += y * ref_cos;
            im -= y * ref_sin;

            measured_samples++;
          }

          sample_index += block_size;
      }

    // ----------------------------
    // Compute output phase
    // ----------------------------
    double output_phase = std::atan2 (im, re) + 0.5 * PI;

    // input phase is 0 (because we generated clean sine reference)
    double phase_diff = -output_phase;

    // unwrap to [-pi, pi]
    while (phase_diff > PI) phase_diff -= 2.0 * PI;
    while (phase_diff < 0) phase_diff += 2.0 * PI;

    // ----------------------------
    // Convert phase → time delay
    // ----------------------------
    double latency_seconds = phase_diff / (2.0 * PI * freq_hz);
    return latency_seconds;
  }
};

inline double
get_time()
{
  /* return timestamp in seconds as double */
  timeval tv;
  gettimeofday (&tv, 0);

  return tv.tv_sec + tv.tv_usec / 1000000.0;
}

float global_f = 0;

complex<double>
expect (SVF::Output output, double f, double freq, double Q, double gain_db)
{
  double A = pow (10, gain_db / 40);

  const double fs = 44100.0;
  const double fc = freq;

  // digital rad/sample
  double w = 2.0 * M_PI * f / fs;

  // prewarped analog cutoff
  double Omega_c = 2.0 * fs * tan(M_PI * fc / fs);

  // bilinear mapping: s = j*Omega
  double Omega = 2.0 * fs * tan(w / 2.0);
  complex<double> s(0.0, Omega);

  // scaled Butterworth prototype
  s /= Omega_c;
  complex<double> H;
  switch (output)
    {
      case SVF::LP:     H = 1.0 / (s*s + s / Q + 1.0);
                        break;
      case SVF::BP:     H = (s / Q) / (s*s + s / Q + 1.0);
                        break;
      case SVF::HP:     H = s * s / (s*s + s / Q + 1.0);
                        break;
      case SVF::AP:     H = (s * s - s / Q + 1.0) /  (s*s + s / Q + 1.0);
                        break;
      case SVF::PEQ:    H = (s * s + s * A / Q + 1.0) / (s * s + s / (A * Q) + 1.0);
                        break;
      case SVF::LSH:    H = A * (s * s + std::sqrt (A) / Q * s + A) / (A * s * s + sqrt (A) / Q * s + 1.0);
                        break;
      case SVF::HSH:    H = A * (A * s * s + std::sqrt (A) / Q * s + 1.0) / (s * s + std::sqrt (A) / Q * s + A);
                        break;
      case SVF::NOTCH:  H = (s * s + 1.0) / (s * s + s / Q + 1.0);
                        break;
      default:          H = 0;
    }

  return H;
}

double
unwrap_phase (double phase, double prev)
{
  while (phase - prev > M_PI)
    phase -= 2.0 * M_PI;

  while (phase - prev < -M_PI)
    phase += 2.0 * M_PI;

  return phase;
}

int
main (int argc, char **argv)
{
  if (argc == 2 && !strcmp (argv[1], "dbg"))
    {
      DistortionDSP distortion_dsp;
      distortion_dsp.set_mode (5);
      distortion_dsp.set_oversample (1);
      distortion_dsp.set_drive (0, true);
      distortion_dsp.set_symmetry (0, true);
      distortion_dsp.set_pre_eq_params (1000, 0, 1, true);
      distortion_dsp.reset (44100);
      distortion_dsp.enable_filters (false);
      const double eps = 0.0011;
      for (double d = -12; d < 12; d += eps)
        {
          float l = d, r = d;
          distortion_dsp.process_block (&l, &r, 1);

          auto cheap_tanh = [] (float x)
            {
              x = std::clamp (x, -3.0f, 3.0f);
              return (x * (27.0f + x * x) / (27.0f + 9.0f * x * x));
            };

          printf ("%.17g %.17g %.17g %.17g\n", d - eps * 0.5, l, r, cheap_tanh (d - eps * 0.5));
        }
    }
  else if (argc == 2 && !strcmp (argv[1], "perf"))
    {
      DistortionDSP distortion_dsp;
      distortion_dsp.set_mode (0);
      distortion_dsp.set_oversample (4);
      distortion_dsp.set_drive (6, true);
      distortion_dsp.set_symmetry (0, true);
      distortion_dsp.reset (48000);

      for (int p = 0; p < 3; p++)
        {
          bool filters = p > 0;
          bool mod = p == 2;
          distortion_dsp.enable_filters (filters);
          const int block_size = 512;
          float left[block_size], right[block_size];

          for (int i = 0; i < block_size; i++)
            {
              left[i] = right[i] = ((i % 100) - 50) / 50;
            }
          double start_t = get_time();
          const int blocks = 50 * 1000;
          for (int b = 0; b < blocks; b++)
            {
              if (mod)
                {
                  const float F[2] = { 440, 1000 };
                  const float G[2] = { 6, 12 };
                  const float Q[2] = { 1, 2 };
                  const float SYM[2] = { -50, 100 };
                  const float SLEW[2] = { 25, 66 };

                  distortion_dsp.set_pre_eq_params (F[b & 1], G[b & 1], Q[b & 1], false);
                  distortion_dsp.set_post_lp (F[b & 1], false);
                  distortion_dsp.set_post_hp (F[b & 1], false);
                  distortion_dsp.set_symmetry (SYM[b & 1], false);
                  distortion_dsp.set_slew (SLEW[b & 1], false);
                }
              distortion_dsp.process_block (left, right, block_size);
            }

          double end_t = get_time();
          double ns_per_sec = 1e9;
          double ns_per_sample = ns_per_sec * (end_t - start_t) / (blocks * block_size);
          printf ("ns/sample %f %s filters, %s modulation\n", ns_per_sample, filters ? "with" : "without", mod ? "with" : "without");
          printf ("                    bogopolyphony = %f\n\n", 1e9 / (ns_per_sample * 48000));
        }
    }
  else if (argc == 2 && !strcmp (argv[1], "latency"))
    {
      DistortionDSP distortion_dsp;
      distortion_dsp.reset (44100);
      distortion_dsp.set_mode (5);
      distortion_dsp.set_oversample (4);
      distortion_dsp.set_drive (-20, true);
      distortion_dsp.set_symmetry (0, true);
      distortion_dsp.set_mix (100, true);
      distortion_dsp.set_pre_eq_params (1000, 12, 1, true);
      for (auto freq_hz : { 110.0, 220.0, 440.0, 880.0, 1000.0, 2000.0, 4000.0, 8000.0, 16000.0 })
        {
          LatencyTester ltest (44100, freq_hz, &distortion_dsp);
          double latency = ltest.run_test (44100, 1024);
          printf ("%.2f Hz -> %f samples latency\n", freq_hz, latency * 44100);
        }
    }
  else if (argc == 3 && !strcmp (argv[1], "sweep"))
    {
      int SR = 44100;
      float mix = atof (argv[2]);
      float buffer[5*SR], buffer2[5*SR], in_freq[5*SR];
      DistortionDSP distortion_dsp;
      double phase = 0;
      for (int i = 0; i < 5*SR; i++)
        {
          double freq = 20*(pow (1000,(double (i)/44100/5)));
          in_freq[i] = freq;
          buffer[i] = sin (phase) * 0.1;
          buffer2[i] = cos (phase) * 0.1;
          phase += freq * 2 * M_PI / 44100;
        }
      distortion_dsp.reset (SR);
      distortion_dsp.set_mode (0);
      distortion_dsp.set_oversample (4);
      distortion_dsp.set_drive (0, true);
      distortion_dsp.set_symmetry (0, true);
      distortion_dsp.set_pre_eq_params (1000, 12, 1, true);
      distortion_dsp.set_mix (mix, true);
      int i = 0;
      while (i < 5 * SR)
        {
          const int TODO = std::min (5 * SR - i, 1024);
          distortion_dsp.process_block (&buffer[i], &buffer2[i], TODO);
          i += TODO;
        }
      for (int i = 0; i < 5*SR; i++)
        {
          printf ("%f %.8f\n", in_freq[i], sqrt (buffer[i] * buffer[i] + buffer2[i] * buffer2[i]));
        }
    }
  else if ((argc == 5 || argc == 6) && (strcmp (argv[1], "sweep-svf") == 0 || strcmp (argv[1], "sweep-svf-mod") == 0))
    {
      int SR = 44100;
      float buffer[5*SR], buffer2[5*SR], in_freq[5*SR], in_phase[5*SR];
      double phase = 0;
      double fade_samples = 250;
      for (int i = 0; i < 5*SR; i++)
        {
          double freq = 20*(pow (1000,(double (i)/44100/5)));
          double fade_in = i < fade_samples ? i / fade_samples : 1.0;
          in_freq[i] = freq;
          in_phase[i] = phase;
          buffer[i] = sin (phase) * fade_in;
          buffer2[i] = cos (phase) * fade_in;
          phase += freq * 2 * M_PI / 44100;
        }
      SVF svf;
      svf.reset (SR);

      SVF::Output output;
      bool found = false;

      for (size_t o = 0; o < svf.output_name.size(); o++)
        if (!strcmp (svf.output_name[o], argv[2]))
          {
            output = SVF::Output (o);
            found = true;
          }
      assert (found);

      float cutoff = atof (argv[3]);
      float Q = atof (argv[4]);
      float gain_db = (argc == 6) ? atof (argv[5]) : 0;

      if (strcmp (argv[1], "sweep-svf") == 0)
        {
          svf.set_params (output, cutoff, 1 / Q, gain_db);
          svf.process_block (output, buffer, buffer2, 5 * SR);
        }
      else
        {
          float freq[5 * SR];
          float Q_inv_in[5 * SR];
          float gain_in[5 * SR];
          for (int i = 0; i < 5 * SR; i++)
            {
              freq[i] = cutoff;
              Q_inv_in[i] = 1 / Q;
              gain_in[i] = gain_db;
            }
          svf.process_mod (output, buffer, buffer2, freq, Q_inv_in, gain_in, 5 * SR);
        }
      double prev_phase = 0;
      float out_phase[5*SR];
      for (int i = 0; i < 5*SR; i++)
        {
          double S     = buffer[i];   // filtered sin sample
          double C     = buffer2[i];  // filtered cos sample
          double phase = unwrap_phase (std::atan2 (S, C) - in_phase[i], prev_phase);

          prev_phase = phase;
          out_phase[i] = phase;
        }
      auto avg_phase = [&] {
        double avg = 0;
        for (int i = 0; i < 5*SR; i++)
          avg += out_phase[i];
        return avg / (5*SR);
      };
      while (avg_phase() > M_PI)
        for (int i = 0; i < 5*SR; i++)
          out_phase[i] -= 2 * M_PI;
      while (avg_phase() < -M_PI)
        for (int i = 0; i < 5*SR; i++)
          out_phase[i] += 2 * M_PI;

      prev_phase = 0;
      for (int i = 2 * fade_samples; i < 5*SR; i++) /* skip first samples (filter fade in) */
        {
          double magnitude = std::sqrt (buffer[i] * buffer[i] + buffer2[i] * buffer2[i]);

          complex<double> H = expect (output, in_freq[i], cutoff, Q, gain_db);
          double expect_phase = unwrap_phase (std::arg (H), prev_phase);
          prev_phase = expect_phase;

          printf ("%f %.8f %.8f %.8f %.8f\n", in_freq[i], magnitude, std::abs (H), out_phase[i], expect_phase);
        }
    }
  else if (argc == 2 && !strcmp (argv[1], "svf-perf"))
    {
      SVF svf;
      svf.reset (48000);

      const int block_size = 512;
      float left[block_size], right[block_size], freq[block_size], Q_inv[block_size], gain_db[block_size];

      for (int i = 0; i < block_size; i++)
        {
          left[i] = right[i] = ((i % 100) - 50) / 50;
          freq[i] = 440 + i;
          Q_inv[i] = 0.7 + i* 0.0001;
          gain_db[i] = 20 + i * 0.0001;
        }
      for (bool modulation : { true, false })
        {
          for (size_t o = 0; o < svf.output_name.size(); o++)
            {
              SVF::Output output = (SVF::Output) o;

              svf.set_params (output, freq[0], 1 / sqrt (2), gain_db[0]);

              double start_t = get_time();
              const int blocks = 50 * 1000;
              if (modulation)
                {
                  for (int b = 0; b < blocks; b++)
                    {
                      svf.process_mod (output, left, right, freq, Q_inv, gain_db, block_size);
                    }
                }
              else
                {
                  for (int b = 0; b < blocks; b++)
                    svf.process_block (output, left, right, block_size);
                }
              global_f = left[0] + right[0]; // avoid optimization

              double end_t = get_time();
              double ns_per_sec = 1e9;
              double ns_per_sample = ns_per_sec * (end_t - start_t) / (blocks * block_size);
              printf ("ns/sample %5.2f, mode %6s, %8s modulation", ns_per_sample, svf.output_name[output], modulation ? "with" : "without");
              printf ("  bogopolyphony = %.2f\n", 1e9 / (ns_per_sample * 48000));
            }
          printf ("\n");
        }
    }
  else if (argc == 2 && !strcmp (argv[1], "svf-lfo"))
    {
      SVF svf;
      int SR = 44100;
      svf.reset (SR);

      float buffer[5*SR], buffer2[5*SR], freq[5*SR], Q_inv[5*SR], gain_db[5*SR];
      double phase = 0;
      for (int i = 0; i < 5*SR; i++)
        {
          float lfo = sin (phase);
          phase += (50.0 * i / (5*SR)) * 2 * M_PI / SR;

          buffer[i] = ((i % 300) - 150)/300.;
          buffer2[i] = buffer[i];
          freq[i] = 1000; // * exp2 (lfo * 4);
          Q_inv[i] = 1;
          gain_db[i] = 20 + lfo * 10; //(i % 50) * 0.0001;
        }
      svf.process_mod<SVF::PEQ> (buffer, buffer2, freq, Q_inv, gain_db, 5 * SR);
      for (int i = 0; i < 5*SR; i++)
        {
          //filter.process_peq_mono (buffer + i, 1000 * exp2 (lfo * 2), 1, 20, 1);
          printf ("%f\n", buffer[i] * 0.1);
        }
    }
  else if (argc == 2 && !strcmp (argv[1], "adaa-table"))
    {
      auto antiderivative_tanh = [] (double x)
        {
          double ax = std::abs(x);
          return ax + std::log1p(std::exp(-2.0 * ax)) - std::log(2.0);
        };
      struct TableRange { static constexpr float range() { return 4; } };
      ADAATable<TableRange, 32> table ([] (double x) { return tanh (x); });
      for (double x = -5; x < 5; x += 0.001)
        printf ("%f %.17g %.17g #f\n", x, table.f (x), tanh (x));
      for (double x = -5; x < 5; x += 0.001)
        printf ("%f %.17g %.17g #F\n", x, table.F (x), antiderivative_tanh (x));
      for (double x = -5; x < 5; x += 0.001)
        printf ("%f %.17g %.17g #d\n", x, (table.F (x + 0.001) - table.F (x))/0.001, tanh (x));
    }
  else if (argc == 2 && !strcmp (argv[1], "adaa-table-sin"))
    {
      struct TableRange { static constexpr float range() { return M_PI; } };
      ADAATable<TableRange, 32, true> table ([] (double x) { return sin (x) + 0.2; });
      for (double x = -15; x < 15; x += 0.001)
        printf ("%f %.17g %.17g #f\n", x, table.f (x), sin (x) + 0.2);
      for (double x = -15; x < 15; x += 0.001)
        printf ("%f %.17g %.17g #F\n", x, table.F (x), -cos (x) + 0.2 * x);
      for (double x = -15; x < 15; x += 0.001)
        printf ("%f %.17g %.17g #d\n", x, (table.F (x + 0.001) - table.F (x))/0.001, sin (x) + 0.2);
    }
  else if (argc == 2 && !strcmp (argv[1], "slew"))
    {
      int SR = 44100;
      float buffer[5*SR], buffer2[5*SR];
      double phase = 0;
      double freq = 1000;
      for (int i = 0; i < 5*SR; i++)
        {
          buffer[i] = buffer2[i] = sin (phase);
          phase += freq * 2 * M_PI / 44100;
        }
      DistortionDSP distortion_dsp;
      distortion_dsp.reset (SR);
      distortion_dsp.set_mode (3);
      distortion_dsp.set_oversample (4);
      distortion_dsp.set_drive (20, true);
      distortion_dsp.set_symmetry (0, true);
      distortion_dsp.enable_filters (false);
      distortion_dsp.set_mix (100, true);
      distortion_dsp.set_slew (50, true);
      int i = 0;
      while (i < 5 * SR)
        {
          const int TODO = std::min (5 * SR - i, 1024);
          distortion_dsp.process_block (&buffer[i], &buffer2[i], TODO);
          i += TODO;
        }
      for (int i = 0; i < 5*SR; i++)
        printf ("%.8f\n", buffer[i] * 0.25);
    }
  else
    {
      int SR = 44100;
      float buffer[5*SR], buffer2[5*SR];
      DistortionDSP distortion_dsp;
      double phase = 0;
      for (int i = 0; i < 5*SR; i++)
        {
          double freq = 20*(pow (1000,(double (i)/44100/5)));
          buffer[i] = buffer2[i] = sin (phase);
          phase += freq * 2 * M_PI / 44100;
        }
      distortion_dsp.reset (SR);
      distortion_dsp.set_mode (2);
      distortion_dsp.set_oversample (4);
      distortion_dsp.set_drive (36, true);
      distortion_dsp.set_symmetry (100, true);
      distortion_dsp.enable_filters (false);
      distortion_dsp.set_mix (100, true);
      int i = 0;
      while (i < 5 * SR)
        {
          const int TODO = std::min (5 * SR - i, 1024);
          distortion_dsp.process_block (&buffer[i], &buffer2[i], TODO);
          i += TODO;
        }
      for (int i = 0; i < 5*SR; i++)
        printf ("%.8f\n", buffer[i] * 0.25);
    }
}
