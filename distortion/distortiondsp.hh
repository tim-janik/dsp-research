#define PANDA_RESAMPLER_HEADER_ONLY

#include "pandaresampler.hh"
#include "log2.hh"
#include "paramsmoother.hh"
#include "adaatable.hh"

using PandaResampler::Resampler2;

#include <array>
#include <cstddef>
#include <cstdint>

/* SVF: The Art of VA Filter Desgin 2.1.2 by Vadim Zavalishin */
class SVF
{
  float g = 0;
  float d = 1;
  float g1 = 1;

  float s1l = 0;
  float s1r = 0;
  float s2l = 0;
  float s2r = 0;

  float m_lp = 0;
  float m_bp = 0;
  float m_hp = 0;

  float cutoff_warp_factor = 0;

  /* tan(x) approximation for x in [0:pi/2] with low relative error (a lot better than 0.1 cent) */
  float
  tan_approx (float x)
  {
    const float c1 = -2.4908436646011975;
    const float c2 = 0.16995685098890287;
    if (x < float (M_PI / 4))
      {
        float x2 = x * x;
        return x * (c1 + c2 * x2) / (c1 + x2);
      }
    else /* for x > pi/4, use (1 / approx (pi/2 - x)) */
      {
        x = float (M_PI / 2) - x;

        float x2 = x * x;
        return (c1 + x2) / (x * (c1 + c2 * x2));
      }
  }

  float
  cutoff_warp (float cutoff)
  {
    return tan_approx (cutoff * cutoff_warp_factor);
  }

  void
  set_params_g_Q_inv (float Q_inv)
  {
    /* instead of using R to compute the SVF parameters, use Q_inv == 1 / Q == 2 * R */
    d = 1.f / (1 + Q_inv*g + g*g);
    g1 = Q_inv + g;
  }
public:
  enum Output {
    LP,
    BP,
    HP,
    AP,
    PEQ,
    LSH,
    HSH,
    NOTCH
  };
  static constexpr std::array<const char *, 8> output_name = { "lp", "bp", "hp", "ap", "peq", "lsh", "hsh", "notch" };
  void
  reset (float sample_rate)
  {
    s1l = s1r = 0;
    s2l = s2r = 0;

    cutoff_warp_factor = M_PI / sample_rate;
  }
  void
  set_params (Output output, float cutoff, float Q_inv, float gain_db)
  {
    if (output == LP || output == BP || output == HP || output == AP || output == NOTCH)
      {
        g = cutoff_warp (cutoff);
        set_params_g_Q_inv (Q_inv);
        if (output == BP)
          m_bp = Q_inv;
        else if (output == AP)
          m_bp = -Q_inv;
      }
    else if (output == PEQ)
      set_peq_params_A (cutoff, Q_inv, powf (10, gain_db / 40));
    else if (output == LSH)
      set_lsh_params_M (cutoff, Q_inv, powf (10, gain_db / 80));
    else if (output == HSH)
      set_hsh_params_M (cutoff, Q_inv, powf (10, gain_db / 80));
    else
      {
        assert (false);
      }
  }
  void
  set_peq_params_A (float cutoff, float Q_inv, float A)
  {
    g = cutoff_warp (cutoff);
    set_params_g_Q_inv (Q_inv / A);

    m_bp = A * Q_inv;
  }
  void
  set_lsh_params_M (float cutoff, float Q_inv, float M)
  {
    float A = M * M;

    m_lp = A * A;
    m_bp = A * Q_inv;
    g = cutoff_warp (cutoff) / M;
    set_params_g_Q_inv (Q_inv);
  }
  void
  set_hsh_params_M (float cutoff, float Q_inv, float M)
  {
    float A = M * M;

    m_bp = A * Q_inv;
    m_hp = A * A;
    g = cutoff_warp (cutoff) * M;
    set_params_g_Q_inv (Q_inv);
  }
  template<Output output>
  void
  process_mod (float *left, float *right, float *freq_in, float *Q_inv_in, float *gain_db_in, uint n_frames)
  {
    if (!n_frames)
      return;

    auto convert_gain_to_F = [gain_db_in] (uint i)
      {
        if constexpr (output == PEQ)
          return powf (10, gain_db_in[i] / 40);
        else if constexpr (output == LSH || output == HSH)
          return powf (10, gain_db_in[i] / 80);
        else
          return 0.f;
      };
    static constexpr int BS = 16;
    float next_F = convert_gain_to_F (0);
    uint i = 0;
    while (i < n_frames)
      {
        float delta_F = 0;
        float F = next_F;

        uint todo = std::min<uint> (BS, n_frames - i);
        if (n_frames - i > BS)
          {
            next_F = convert_gain_to_F (i + BS);
            delta_F = (next_F - F) / BS;
          }
        else if (todo > 1)
          {
            next_F = convert_gain_to_F (i + todo - 1);
            delta_F = (next_F - F) / (todo - 1);
          }
        for (uint j = 0; j < todo; j++)
          {
            const float freq = freq_in[i + j];
            const float Q_inv = Q_inv_in[i + j];

            if constexpr (output == PEQ)
              {
                set_peq_params_A (freq, Q_inv, F);
                F += delta_F;
              }
            else if constexpr (output == LSH)
              {
                set_lsh_params_M (freq, Q_inv, F);
                F += delta_F;
              }
            else if constexpr (output == HSH)
              {
                set_hsh_params_M (freq, Q_inv, F);
                F += delta_F;
              }
            else
              {
                g = cutoff_warp (freq);

                set_params_g_Q_inv (Q_inv);

                if (output == BP)
                  m_bp = Q_inv;
                if (output == AP)
                  m_bp = -Q_inv;
              }
            process_s<output> (left + i + j, right + i + j);
          }
        i += todo;
      }
  }
  template<Output output>
  float process (float l)
  {
    float dummy = 0;
    process_s<output> (&l, &dummy);
    return l;
  }
  template<Output output>
  void process_s (float *l, float *r) __restrict__
  {
    // hp
    float hpl, hpr;
    hpl = (*l - g1*s1l - s2l) * d;
    hpr = (*r - g1*s1r - s2r) * d;

    // first integrator
    float v1l, v1r;
    v1l = g * hpl;
    v1r = g * hpr;

    float bpl, bpr;
    bpl = v1l + s1l;
    bpr = v1r + s1r;
    s1l = bpl + v1l;
    s1r = bpr + v1r;

    // second integrator
    float v2l, v2r;
    v2l = g * bpl;
    v2r = g * bpr;

    float lpl, lpr;
    lpl = v2l + s2l;
    lpr = v2r + s2r;
    s2l = lpl + v2l;
    s2r = lpr + v2r;

    if (output == LP)
      {
        *l = lpl;
        *r = lpr;
      }
    else if (output == BP)
      {
        *l = bpl * m_bp;
        *r = bpr * m_bp;
      }
    else if (output == HP)
      {
        *l = hpl;
        *r = hpr;
      }
    else if (output == AP)
      {
        *l = lpl + bpl * m_bp + hpl;
        *r = lpr + bpr * m_bp + hpr;
      }
    else if (output == PEQ)
      {
        *l = lpl + hpl + bpl * m_bp;
        *r = lpr + hpr + bpr * m_bp;
      }
    else if (output == LSH)
      {
        *l = m_lp * lpl + bpl * m_bp + hpl;
        *r = m_lp * lpr + bpr * m_bp + hpr;
      }
    else if (output == HSH)
      {
        *l = lpl + bpl * m_bp + m_hp * hpl;
        *r = lpr + bpr * m_bp + m_hp * hpr;
      }
    else if (output == NOTCH)
      {
        *l = lpl + hpl;
        *r = lpr + hpr;
      }
    else
      {
        assert (false);
      }
  }
  template<Output output>
  void
  process_block (float *left, float *right, uint n_samples)
  {
    for (uint i = 0; i < n_samples; i++)
      {
        process_s<output> (left + i, right + i);
      }
  }
  void
  process_block (SVF::Output output, float *left, float *right, uint n_samples)
  {
    switch (output)
      {
        case LP:    process_block<LP> (left, right, n_samples);
                    break;
        case BP:    process_block<BP> (left, right, n_samples);
                    break;
        case HP:    process_block<HP> (left, right, n_samples);
                    break;
        case AP:    process_block<AP> (left, right, n_samples);
                    break;
        case PEQ:   process_block<PEQ> (left, right, n_samples);
                    break;
        case LSH:   process_block<LSH> (left, right, n_samples);
                    break;
        case HSH:   process_block<HSH> (left, right, n_samples);
                    break;
        case NOTCH: process_block<NOTCH> (left, right, n_samples);
                    break;
        default:    assert (false);
      }
  }
  void
  process_mod (SVF::Output output, float *left, float *right, float *freq_in, float *Q_inv_in, float *gain_db_in, uint n_frames)
  {
    switch (output)
      {
        case LP:    process_mod<LP> (left, right, freq_in, Q_inv_in, gain_db_in, n_frames);
                    break;
        case BP:    process_mod<BP> (left, right, freq_in, Q_inv_in, gain_db_in, n_frames);
                    break;
        case HP:    process_mod<HP> (left, right, freq_in, Q_inv_in, gain_db_in, n_frames);
                    break;
        case AP:    process_mod<AP> (left, right, freq_in, Q_inv_in, gain_db_in, n_frames);
                    break;
        case PEQ:   process_mod<PEQ> (left, right, freq_in, Q_inv_in, gain_db_in, n_frames);
                    break;
        case LSH:   process_mod<LSH> (left, right, freq_in, Q_inv_in, gain_db_in, n_frames);
                    break;
        case HSH:   process_mod<HSH> (left, right, freq_in, Q_inv_in, gain_db_in, n_frames);
                    break;
        case NOTCH: process_mod<NOTCH> (left, right, freq_in, Q_inv_in, gain_db_in, n_frames);
                    break;
        default:    assert (false);
      }
  }
};

template <size_t max_delay_samples_pow2>
class StereoDelay
{
  static_assert ((max_delay_samples_pow2 & (max_delay_samples_pow2 - 1)) == 0,
                 "max_delay_samples_pow2 must be a power of two");

public:
  StereoDelay (size_t delay_samples = 0)
  {
    set_delay (delay_samples);
    reset();
  }

  void
  set_delay (size_t delay_samples)
  {
    assert (delay_samples < max_delay_samples_pow2);
    m_delay_samples = delay_samples;
  }

  void
  reset()
  {
    m_buffer_l.fill (0.0f);
    m_buffer_r.fill (0.0f);
    m_write_index = 0;
  }

  inline void
  process_sample (float in_l, float in_r, float& out_l, float& out_r)
  {
    constexpr size_t mask = max_delay_samples_pow2 - 1;

    const size_t read_index = (m_write_index - m_delay_samples) & mask;

    m_buffer_l[m_write_index] = in_l;
    m_buffer_r[m_write_index] = in_r;

    out_l = m_buffer_l[read_index];
    out_r = m_buffer_r[read_index];

    m_write_index = (m_write_index + 1) & mask;
  }

private:
  std::array<float, max_delay_samples_pow2> m_buffer_l{};
  std::array<float, max_delay_samples_pow2> m_buffer_r{};

  size_t m_delay_samples = 0;
  size_t m_write_index = 0;
};

class DistortionDSP
{
  struct ADAATables
  {
    static constexpr int N_TABLES = 31; // odd to have the center table represent a linear function
    std::array<std::unique_ptr<ADAATable<4, 1024>>, N_TABLES> tables;
    ADAATables()
    {
      for (size_t i = 0; i < N_TABLES; i++)
        {
          auto distort_asymmetric = [&] (double x, double s)
             {
              double k  = s * 10;  // symmetry between -1 and 1 scales to k between -10 and 10
              double kx = k * x;
              if (std::abs (kx) < 0.001)
                return x;
              else
                return (x / (1 - exp (-kx)) - 1./k)*2;
            };
          float symmetry = (i / (N_TABLES - 1.0)) * 2 - 1;
          tables[i] = std::make_unique<ADAATable<4, 1024>> (
            [&] (double x) { return distort_asymmetric (tanh (x), symmetry); }
          );
        }
    }
    static ADAATables& the()
    {
      static ADAATables instance;
      return instance;
    }
  };
  ADAATables& adaa_tables { ADAATables::the() };

  std::unique_ptr<Resampler2> up_left;
  std::unique_ptr<Resampler2> up_right;
  std::unique_ptr<Resampler2> down_left;
  std::unique_ptr<Resampler2> down_right;

  static constexpr int MAX_OVERSAMPLE = 8;
  int oversample = -1;
  float mix = 1;
  int mode = 0;
  int   last_table = -1;
  float last_left = 0;
  float last_right = 0;
  float last_left_F_1 = 0;
  float last_left_F_2 = 0;
  float last_right_F_1 = 0;
  float last_right_F_2 = 0;

  int over_delay = 0;
  std::array<float, MAX_OVERSAMPLE> left_over_delay_history {};
  std::array<float, MAX_OVERSAMPLE> right_over_delay_history {};

  bool filters_enabled = true;

  int    sample_rate = 44100;

  SVF                                       pre_eq_filter;
  ParamSmoother<SmootherType::logarithmic>  pre_eq_freq_smoother  { 1000 };
  ParamSmoother<SmootherType::logarithmic>  pre_eq_Q_inv_smoother { 1 };
  ParamSmoother<SmootherType::linear>       pre_eq_gain_smoother  { 6 };


  SVF                                       post_lp_filter;
  ParamSmoother<SmootherType::logarithmic>  post_lp_freq_smoother { 20000 };

  SVF                                       post_hp_filter;
  ParamSmoother<SmootherType::logarithmic>  post_hp_freq_smoother { 20 };

  StereoDelay<64> dry_delay;

  ParamSmoother<SmootherType::linear>       symmetry_smoother     { 0 };
  ParamSmoother<SmootherType::logarithmic>  drive_factor_smoother { 1 };

  // https://www.musicdsp.org/en/latest/Other/238-rational-tanh-approximation.html
  float
  cheap_tanh (float x)
  {
    x = std::clamp (x, -3.0f, 3.0f);

    return (x * (27.0f + x * x) / (27.0f + 9.0f * x * x));
  }

  double
  cheap_tanh_antiderivative (double x)
  {
    return (1.0/18.0) * x * x + (4.0/3.0) * log (x * x + 3);
  }
  float
  cheap_tanh_antiderivative_approx (float x)
  {
    float ax = std::abs (x);
    if (ax > 3)
      {
        return ax + 0.813208866384f;
      }
    else
      {
        const float ln2 = 0.693147180559945f;
        return (1.f/18.f) * x * x + (4.f / 3.f * ln2) * fast_log2 (x * x + 3);
      }
  }
public:
  void
  reset (int sample_rate)
  {
    pre_eq_filter.reset (sample_rate);

    post_lp_filter.reset (sample_rate);
    post_hp_filter.reset (sample_rate);

    pre_eq_freq_smoother.reset (sample_rate, 0.025);
    pre_eq_Q_inv_smoother.reset (sample_rate, 0.025);
    pre_eq_gain_smoother.reset (sample_rate, 0.025);

    post_lp_freq_smoother.reset (sample_rate, 0.025);
    post_hp_freq_smoother.reset (sample_rate, 0.025);

    symmetry_smoother.reset (sample_rate, 0.025);

    this->sample_rate = sample_rate;
    left_over_delay_history.fill (0);
    right_over_delay_history.fill (0);
    dry_delay.reset();

    last_table = -1;
  }
  void
  set_pre_eq_params (float freq, float gain, float Q, bool now)
  {
    pre_eq_freq_smoother.set_target (freq, now);
    pre_eq_Q_inv_smoother.set_target (1 / Q, now);
    pre_eq_gain_smoother.set_target (gain, now);
  }
  void
  set_post_lp (float lp_freq, bool now)
  {
    post_lp_freq_smoother.set_target (lp_freq, now);
  }
  void
  set_post_hp (float hp_freq, bool now)
  {
    post_hp_freq_smoother.set_target (hp_freq, now);
  }
  void
  set_oversample (int new_oversample)
  {
    if (new_oversample != oversample)
      {
        oversample = new_oversample;
        up_left = std::make_unique<Resampler2> (Resampler2::UP, oversample, Resampler2::PREC_72DB);
        up_right = std::make_unique<Resampler2> (Resampler2::UP, oversample, Resampler2::PREC_72DB);
        down_left = std::make_unique<Resampler2> (Resampler2::DOWN, oversample, Resampler2::PREC_72DB);
        down_right = std::make_unique<Resampler2> (Resampler2::DOWN, oversample, Resampler2::PREC_72DB);

        /* delay compensation:
         *  - oversampling path     -> delay compensation to compensate for PandaResampler fractional delay
         *  - non-oversampled path  -> whole sample delay compensation
         *
         * we don't compensate for ADAA delay (0.5 samples in the oversampled path), but this is small
         * so the frequency response is still relatively flat for wet/dry mix
         */
        double delay = up_left->delay() / oversample + down_left->delay();

        int best_i = 0;
        double best_err = 1e10;
        for (int i = 0; i < oversample; i++)
          {
            double d = delay + double (i) / oversample;
            double err = std::fabs (d - std::round (d));

            if (err < best_err)
              {
                best_err = err;
                best_i   = i;
              }
          }
        over_delay = best_i;

        left_over_delay_history.fill (0);
        right_over_delay_history.fill (0);

        dry_delay.set_delay (std::round (delay + double (over_delay) / oversample));
        dry_delay.reset();
      }
  }
  void
  set_drive (float drive, bool now)
  {
    drive_factor_smoother.set_target (exp10f (drive * (1/20.f)), now);
  }
  void
  set_symmetry (float new_symmetry, bool now)
  {
    symmetry_smoother.set_target (new_symmetry, now);
  }
  void
  set_mode (int new_mode)
  {
    mode = new_mode;
  }
  void
  set_mix (float percent)
  {
    mix = std::clamp (percent * 0.01, 0.0, 1.0);
  }
  void
  enable_filters (bool enable)
  {
    filters_enabled = enable;
  }
  void
  process_block (float *left_in, float *right_in, int n_samples)
  {
    if (!n_samples)
      {
        /*
         * oversample delay compensation requires at least one sample to be processed
         * so that we can always update the delay history buffer
         */
        return;
      }

    float dry_delay_left[n_samples];
    float dry_delay_right[n_samples];
    for (int i = 0; i < n_samples; i++)
      dry_delay.process_sample (left_in[i], right_in[i], dry_delay_left[i], dry_delay_right[i]);

    if (filters_enabled)
      {
        if (pre_eq_freq_smoother.is_constant() && pre_eq_Q_inv_smoother.is_constant() && pre_eq_gain_smoother.is_constant())
          {
            /* FAST: no smoothing case */
            pre_eq_filter.set_params (SVF::PEQ, pre_eq_freq_smoother.get_next(), pre_eq_Q_inv_smoother.get_next(), pre_eq_gain_smoother.get_next());
            pre_eq_filter.process_block (SVF::PEQ, left_in, right_in, n_samples);
          }
        else
          {
            /* SLOW: at least one parameter is smoothing */
            float freq_peq[n_samples];
            float Q_inv_peq[n_samples];
            float gain_peq[n_samples];

            pre_eq_freq_smoother.process_block (freq_peq, n_samples);
            pre_eq_Q_inv_smoother.process_block (Q_inv_peq, n_samples);
            pre_eq_gain_smoother.process_block (gain_peq, n_samples);

            pre_eq_filter.process_mod (SVF::PEQ, left_in, right_in, freq_peq, Q_inv_peq, gain_peq, n_samples);
          }
      }

    if (drive_factor_smoother.is_constant())
      {
        float drive_factor = drive_factor_smoother.get_next();
        for (int i = 0; i < n_samples; i++)
          {
            left_in[i] *= drive_factor;
            right_in[i] *= drive_factor;
          }
      }
    else
      {
        for (int i = 0; i < n_samples; i++)
          {
            float drive_factor = drive_factor_smoother.get_next();
            left_in[i] *= drive_factor;
            right_in[i] *= drive_factor;
          }
      }

    float left_over_raw[over_delay + oversample * n_samples];
    float right_over_raw[over_delay + oversample * n_samples];
    float *left_over = left_over_raw + over_delay;
    float *right_over = right_over_raw + over_delay;

    up_left->process_block (left_in, n_samples, left_over);
    up_right->process_block (right_in, n_samples, right_over);

    float *left = left_over;
    float *right = right_over;

    if (mode == 0)
      {
        float symmetry[n_samples];
        symmetry_smoother.process_block (symmetry, n_samples);

        float left_pre_F[n_samples * oversample];
        float right_pre_F[n_samples * oversample];
        for (size_t i = 0; i < n_samples * oversample; i++)
          {
            left[i] = std::clamp (left[i], -10.f, 10.f);
            right[i] = std::clamp (right[i], -10.f, 10.f);
          }
        for (size_t i = 0; i < n_samples * oversample; i++)
          {
            left_pre_F[i] = cheap_tanh_antiderivative_approx (left[i]);
            right_pre_F[i] = cheap_tanh_antiderivative_approx (right[i]);
          }

        for (size_t i = 0; i < n_samples * oversample; i++)
          {
            auto adaa = [&] (float x, float last_x, float F, float last_F)
              {
                /* ADAA quotient is (F - last_F) / (x - last_x)
                 *
                 * This is problematic if F and last_F are very close, because
                 * then float cancellation will remove the significant bits, so
                 * ADAA approximation will be inaccurate.
                 *
                 * We could do everything in double precision but this would
                 * be slow.
                 *
                 * Insead, we use a rather high epsilon, because in real world
                 * signals if x and last_x are very similar then the ADAA value
                 * is close to the cheap_tanh value anyway.
                 */
                const float epsilon = 0.001f;

                float delta = x - last_x;
                if (std::abs (delta) > epsilon)
                  return (F - last_F) / delta;
                else
                  return cheap_tanh (0.5f * (x + last_x));
              };

            float l = left[i];
            float left_F = left_pre_F[i];
            left[i] = adaa (l, last_left, left_F, last_left_F_1);
            last_left = l;
            last_left_F_1 = left_F;

            float r = right[i];
            float right_F = right_pre_F[i];
            right[i] = adaa (r, last_right, right_F, last_right_F_1);
            last_right = r;
            last_right_F_1 = right_F;

            auto distort = [&] (float x, float s)
              {
                float k  = s * 10;  // symmetry between -1 and 1 scales to k between -10 and 10
                float kx = k * x;
                if (std::abs (kx) < 0.001f)
                  return x;
                else
                  return (x / (1 - exp (-kx)) - 1.f/k)*2;
              };
            left[i] = distort (left[i], symmetry[i / oversample] * 0.01f);
            right[i] = distort (right[i], symmetry[i / oversample] * 0.01f);
          }
        goto out;
      }
    if (mode == 1)
      {
        for (size_t i = 0; i < n_samples * oversample; i++)
          {
            left[i] = std::sin (left[i]);
            right[i] = std::sin (right[i]);
          }
      }
    if (mode == 2)
      {
        if (symmetry_smoother.is_constant())
          {
            process_with_symmetry (left, right, n_samples * oversample, symmetry_smoother.get_next());
          }
        else
          {
            for (size_t i = 0; i < n_samples * oversample; i += oversample)
              process_with_symmetry (left + i, right + i, oversample, symmetry_smoother.get_next());
          }
      }
#if 0
    if (mode == 6)
      {
        for (size_t i = 0; i < n_samples * oversample; i++)
          {
            left[i]  = cheap_tanh (drive_factor * left[i]);
            right[i] = cheap_tanh (drive_factor * right[i]);
          }
        goto out;
      }
    for (size_t i = 0; i < n_samples * oversample; i++)
      {
        if (mode == 0)
          {
            left[i] = left[i]>0?tanh (drive_factor * left[i]):tanh(drive_factor * left[i] / neg_scale) * neg_scale;
            right[i] = right[i]>0?tanh (drive_factor * right[i]):tanh(drive_factor * right[i] / neg_scale) * neg_scale;
          }
        else if (mode == 1)
          {
            left[i] = tanh(drive_factor * left[i]) / (1-s*tanh(drive_factor * left[i]))*(1-s);
            right[i] = tanh(drive_factor * right[i]) / (1-s*tanh(drive_factor * right[i]))*(1-s);
          }
        else if (mode == 2)
          {
            float bias = std::clamp (symmetry * 0.01f, -1.f, 1.f) * 0.5f;
            auto deriv = [] (float bias) { return (tanh(bias+0.001)-tanh(bias))/0.001; };
            float norm = 1/deriv (bias);
            left[i] = (tanh(drive_factor * left[i] + bias) - tanh (bias)) * norm;
            right[i] = (tanh(drive_factor * right[i] + bias) - tanh (bias)) * norm;
          }
        else if (mode == 3)
          {
            float sl = 0;
            float sr = 0;
            for (int j = 0; j < 10; j++)
              {
                float frac = j / 10.;

                sl += tanh ((last_left * (1 - frac) + left[i] * frac) * drive_factor);
                sr += tanh ((last_right * (1 - frac) + right[i] * frac) * drive_factor);
              }
            last_right = right[i];
            last_left = left[i];

            left[i] = sl * 0.1f;
            right[i] = sr * 0.1f;
          }
        else if (mode == 4)
          {
            auto F = [] (double x) { return log(cosh(x)); };
            auto adaa = [&] (double x, double last_x) { 
              x = std::clamp (x, -10., 10.);
              last_x = std::clamp (last_x, -10., 10.);
              if (fabs (x - last_x) > 0.00001)
                return (F(x) - F(last_x)) / (x - last_x);
              else
                return tanh(x);
            };
            float l = adaa (left[i] * drive_factor, last_left * drive_factor);
            float r = adaa (right[i] * drive_factor, last_right * drive_factor);
            last_left = left[i];
            last_right = right[i];
            left[i] = l;
            right[i] = r;
          }
      }
#endif
out:
    std::copy_n (left_over_delay_history.begin(), over_delay, left_over_raw);
    std::copy_n (right_over_delay_history.begin(), over_delay, right_over_raw);
    std::copy_n (&left[n_samples * oversample - over_delay], over_delay, left_over_delay_history.begin());
    std::copy_n (&right[n_samples * oversample - over_delay], over_delay, right_over_delay_history.begin());

    down_left->process_block (left_over_raw, oversample * n_samples, left_in);
    down_right->process_block (right_over_raw, oversample * n_samples, right_in);

    if (filters_enabled)
      {
        constexpr double BUTTERWORTH_Q = M_SQRT1_2; /* 1 / sqrt (2) */
        constexpr float  BUTTERWORTH_Q_INV = 1 / BUTTERWORTH_Q;

        float Q_inv_lp_hp[n_samples];
        if (!post_lp_freq_smoother.is_constant() || !post_hp_freq_smoother.is_constant())
          {
            std::fill_n (Q_inv_lp_hp, n_samples, BUTTERWORTH_Q_INV);
          }

        if (post_lp_freq_smoother.is_constant())
          {
            /* FAST: no smoothing case */
            post_lp_filter.set_params (SVF::LP, post_lp_freq_smoother.get_next(), BUTTERWORTH_Q_INV, 0);
            post_lp_filter.process_block (SVF::LP, left_in, right_in, n_samples);
          }
        else
          {
            /* SLOW: with lowpass frequency smoothing */
            float freq_lp[n_samples];

            post_lp_freq_smoother.process_block (freq_lp, n_samples);
            post_lp_filter.process_mod (SVF::LP, left_in, right_in, freq_lp, Q_inv_lp_hp, nullptr, n_samples);
          }

        if (post_hp_freq_smoother.is_constant())
          {
            /* FAST: no smoothing case */
            post_hp_filter.set_params (SVF::HP, post_hp_freq_smoother.get_next(), BUTTERWORTH_Q_INV, 0);
            post_hp_filter.process_block (SVF::HP, left_in, right_in, n_samples);
          }
        else
          {
            /* SLOW: with highpass frequency smoothing */
            float freq_hp[n_samples];

            post_hp_freq_smoother.process_block (freq_hp, n_samples);
            post_hp_filter.process_mod (SVF::HP, left_in, right_in, freq_hp, Q_inv_lp_hp, nullptr, n_samples);
          }
      }

    for (int i = 0; i < n_samples; i++)
      {
        left_in[i] = dry_delay_left[i] + mix * (left_in[i] - dry_delay_left[i]);
        right_in[i] = dry_delay_right[i] + mix * (right_in[i] - dry_delay_right[i]);
      }
  }

  void
  process_with_symmetry (float *left, float *right, int n_samples, float symmetry)
  {
    // map symmetry [-100..100] to table index [0..N_TABLES - 1]
    float ftable_index = (symmetry * 0.01f + 1) / 2 * (adaa_tables.N_TABLES - 1);
    ftable_index = std::max (0.f, ftable_index);
    int table_index = int (ftable_index);
    table_index = std::min (table_index, adaa_tables.N_TABLES - 2);
    float frac = ftable_index - table_index;

    auto& table_1 = *adaa_tables.tables[table_index];
    auto& table_2 = *adaa_tables.tables[table_index + 1];

    if (last_table != table_index)
      {
        last_left_F_1 = table_1.F (last_left);
        last_right_F_1 = table_1.F (last_right);
        last_left_F_2 = table_2.F (last_left);
        last_right_F_2 = table_2.F (last_right);
        last_table = table_index;
      }
    for (int i = 0; i < n_samples; i++)
      {
        float l = left[i];
        float r = right[i];

        auto adaa = [&] (float x, float last_x, float F, float last_F, auto& table)
          {
            /* ADAA quotient is (F - last_F) / (x - last_x)
             *
             * This is problematic if F and last_F are very close, because
             * then float cancellation will remove the significant bits, so
             * ADAA approximation will be inaccurate.
             *
             * We could do everything in double precision but this would
             * be slow.
             *
             * Insead, we use a rather high epsilon, because in real world
             * signals if x and last_x are very similar then the ADAA value
             * is close to the sin value anyway.
             */
            const float epsilon = 0.001f;

            float delta = x - last_x;
            if (std::abs (delta) > epsilon)
              return (F - last_F) / delta;
            else
              return table.f (0.5f * (x + last_x));
          };

        float left_F_1 = table_1.F (l);
        float left_F_2 = table_2.F (l);
        left[i] = adaa (l, last_left, left_F_1, last_left_F_1, table_1) * (1 - frac) +
                  adaa (l, last_left, left_F_2, last_left_F_2, table_2) * frac;
        last_left = l;
        last_left_F_1 = left_F_1;
        last_left_F_2 = left_F_2;

        float right_F_1 = table_1.F (r);
        float right_F_2 = table_2.F (r);
        right[i] = adaa (r, last_right, right_F_1, last_right_F_1, table_1) * (1 - frac) +
                   adaa (r, last_right, right_F_2, last_right_F_2, table_2) * frac;
        last_right = r;
        last_right_F_1 = right_F_1;
        last_right_F_2 = right_F_2;
      }
  }
};


