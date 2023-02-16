// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl.html
#ifndef SPECTMORPH_LADDER_VCF_HH
#define SPECTMORPH_LADDER_VCF_HH

#include "pandaresampler.hh"

#include <array>

namespace SpectMorph {

using PandaResampler::Resampler2;

enum class LadderVCFMode { LP1, LP2, LP3, LP4 };

template<bool OVERSAMPLE, bool NON_LINEAR>
class LadderVCF
{
  struct Channel {
    float x1, x2, x3, x4;
    float y1, y2, y3, y4;

    Resampler2 res_up   { Resampler2::UP,   2, Resampler2::PREC_72DB };
    Resampler2 res_down { Resampler2::DOWN, 2, Resampler2::PREC_72DB };
  };
  std::array<Channel, 2> channels;
  LadderVCFMode mode;
  float pre_scale, post_scale;
  float rate;
  float freq_ = 440;
  float reso_ = 0;
  float drive_ = 0;

public:
  LadderVCF()
  {
    reset();
    set_mode (LadderVCFMode::LP4);
    set_scale (1, 1);
    set_rate (48000);
  }
  void
  set_mode (LadderVCFMode new_mode)
  {
    mode = new_mode;
  }
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
  set_drive (float drive)
  {
    drive_ = drive;
  }
  void
  set_scale (float pre, float post)
  {
    pre_scale = pre;
    post_scale = post;
  }
  void
  set_rate (float r)
  {
    rate = r;
  }
  void
  reset()
  {
    for (auto& c : channels)
      {
        c.x1 = c.x2 = c.x3 = c.x4 = 0;
        c.y1 = c.y2 = c.y3 = c.y4 = 0;

        c.res_up.reset();
        c.res_down.reset();
      }
  }
  float
  distort (float x)
  {
    if (NON_LINEAR)
      {
        /* shaped somewhat similar to tanh() and others, but faster */
        x = std::clamp (x, -1.0f, 1.0f);

        return x - x * x * x * (1.0f / 3);
      }
    else
      {
        return x;
      }
  }
private:
  /*
   * This ladder filter implementation is mainly based on
   *
   * Välimäki, Vesa & Huovilainen, Antti. (2006).
   * Oscillator and Filter Algorithms for Virtual Analog Synthesis.
   * Computer Music Journal. 30. 19-31. 10.1162/comj.2006.30.2.19.
   */
  template<LadderVCFMode MODE, bool STEREO> inline void
  run (float *left, float *right, float fc, float res)
  {
    const float pi = M_PI;
    fc = pi * fc;
    const float g = 0.9892f * fc - 0.4342f * fc * fc + 0.1381f * fc * fc * fc - 0.0202f * fc * fc * fc * fc;
    const float b0 = g * (1 / 1.3f);
    const float b1 = g * (0.3f / 1.3f);
    const float a1 = g - 1;

    res *= 1.0029f + 0.0526f * fc - 0.0926f * fc * fc + 0.0218f * fc * fc * fc;

    constexpr uint oversample_count = OVERSAMPLE ? 2 : 1;
    for (uint os = 0; os < oversample_count; os++)
      {
        for (uint i = 0; i < (STEREO ? 2 : 1); i++)
          {
            float &value = i == 0 ? left[os] : right[os];

            Channel& c = channels[i];
            const float x = value * pre_scale;
            const float g_comp = 0.5f; // passband gain correction
            const float x0 = distort (x - (c.y4 - g_comp * x) * res * 4);

            c.y1 = b0 * x0 + b1 * c.x1 - a1 * c.y1;
            c.x1 = x0;

            c.y2 = b0 * c.y1 + b1 * c.x2 - a1 * c.y2;
            c.x2 = c.y1;

            c.y3 = b0 * c.y2 + b1 * c.x3 - a1 * c.y3;
            c.x3 = c.y2;

            c.y4 = b0 * c.y3 + b1 * c.x4 - a1 * c.y4;
            c.x4 = c.y3;

            switch (MODE)
              {
                case LadderVCFMode::LP1:
                  value = c.y1 * post_scale;
                  break;
                case LadderVCFMode::LP2:
                  value = c.y2 * post_scale;
                  break;
                case LadderVCFMode::LP3:
                  value = c.y3 * post_scale;
                  break;
                case LadderVCFMode::LP4:
                  value = c.y4 * post_scale;
                  break;
                default:
                  assert (false);
              }
          }
      }
  }
  template<LadderVCFMode MODE, bool STEREO> inline void
  do_run_block (uint          n_samples,
                float        *left,
                float        *right,
                const float  *freq_in,
                const float  *reso_in,
                const float  *drive_in)
  {
    float over_samples_left[2 * n_samples];
    float over_samples_right[2 * n_samples];
    float freq_scale = OVERSAMPLE ? 0.5 : 1.0;
    float nyquist    = rate * 0.5;

    if (OVERSAMPLE)
      {
        channels[0].res_up.process_block (left, n_samples, over_samples_left);
        if (STEREO)
          channels[1].res_up.process_block (right, n_samples, over_samples_right);
      }

    float fc = freq_ * freq_scale / nyquist;
    float res = reso_;

    for (uint i = 0; i < n_samples; i++)
      {
        float mod_fc = fc;
        float mod_res = res;

        if (freq_in)
          mod_fc = freq_in[i] * freq_scale / nyquist;

        if (reso_in)
          mod_res = reso_in[i];

        mod_fc  = std::clamp<float> (mod_fc, 0.0f, 1.0f);
        mod_res = std::clamp<float> (mod_res, 0.0f, 1.0f);

        if (OVERSAMPLE)
          {
            const uint over_pos = i * 2;

            run<MODE, STEREO> (over_samples_left + over_pos, over_samples_right + over_pos, mod_fc, mod_res);
          }
        else
          {
            run<MODE, STEREO> (left + i, right + i, mod_fc, mod_res);
          }
      }
    if (OVERSAMPLE)
      {
        channels[0].res_down.process_block (over_samples_left, 2 * n_samples, left);
        if (STEREO)
          channels[1].res_down.process_block (over_samples_right, 2 * n_samples, right);
      }
  }
  template<LadderVCFMode MODE> inline void
  run_block_mode (uint          n_samples,
                  float        *left,
                  float        *right,
                  const float  *freq_in,
                  const float  *reso_in,
                  const float  *drive_in)
  {
    if (right) // stereo?
      do_run_block<MODE, true> (n_samples, left, right, freq_in, reso_in, drive_in);
    else
      do_run_block<MODE, false> (n_samples, left, right, freq_in, reso_in, drive_in);
  }
public:
  void
  run_block (uint         n_samples,
             float       *left,
             float       *right,
             const float *freq_in = nullptr,
             const float *reso_in = nullptr,
             const float *drive_in = nullptr)
  {
    switch (mode)
    {
      case LadderVCFMode::LP4: run_block_mode<LadderVCFMode::LP4> (n_samples, left, right, freq_in, reso_in, drive_in);
                               break;
      case LadderVCFMode::LP3: run_block_mode<LadderVCFMode::LP3> (n_samples, left, right, freq_in, reso_in, drive_in);
                               break;
      case LadderVCFMode::LP2: run_block_mode<LadderVCFMode::LP2> (n_samples, left, right, freq_in, reso_in, drive_in);
                               break;
      case LadderVCFMode::LP1: run_block_mode<LadderVCFMode::LP1> (n_samples, left, right, freq_in, reso_in, drive_in);
                               break;
    }
  }
};

// fast linear model of the filter
typedef LadderVCF<false, false> LadderVCFLinear;

// slow but accurate non-linear model of the filter (uses oversampling)
typedef LadderVCF<true,  true>  LadderVCFNonLinear;

// fast non-linear version (no oversampling), may have aliasing
typedef LadderVCF<false, true>  LadderVCFNonLinearCheap;

} // SpectMorph

#endif // __BSE_DEVICES_LADDER_VCF_HH__
