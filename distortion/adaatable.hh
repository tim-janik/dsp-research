#include <functional>

template<class RangeClass, int BINS, bool PERIODIC = false>
class ADAATable
{
  static constexpr float RANGE = RangeClass::range();

  std::array<float, BINS + 1> f_;
  std::array<float, BINS + 1> F_;
  float                       period_integral_;

  struct Wrap
  {
    int   ibin;
    float frac;
    int   periods;
  };

  Wrap
  periodic_wrap (float x) const
  {
    Wrap wrap;

    /* map periodicity, [-RANGE..RANGE] -> [0..1] */
    x = (x + RANGE) * (1.f / (2.f * RANGE));

    wrap.periods = std::floor (x);
    x -= wrap.periods;

    float fbin = x * BINS;

    // floating point arithmetic does not guarantee fbin < BINS at this point
    if (fbin >= BINS)
      {
        wrap.ibin = BINS - 1;
        wrap.frac = 1;
      }
    else
      {
        wrap.ibin = fbin;
        wrap.frac = fbin - wrap.ibin;
      }
    return wrap;
  }
public:
  ADAATable (std::function<double(double)> f)
  {
    // build function f(x) table
    for (size_t i = 0; i < BINS + 1; i++)
      f_[i] = f ((i / double (BINS) * 2 - 1) * RANGE);

    // build antiderivative F(x) table
    std::vector<double> F (BINS + 1);

    double Fx = 0;
    double dx = 2.0 * RANGE / BINS;
    F[0] = Fx;

    for (size_t i = 0; i < BINS; i++)
      {
        Fx += 0.5f * dx * (f_[i] + f_[i+1]);
        F[i + 1] = Fx;
      }
    period_integral_ = F[BINS];

    /* shift antiderivative so that smallest elements are close to zero
     *   -> better resolution for small floats
     */
    auto min_F = *std::min_element (F.begin(), F.end());
    for (size_t i = 0; i < BINS + 1; i++)
      F_[i] = F[i] - min_F;
  }

  float
  f (float x) const
  {
    int ibin;
    float frac;
    if constexpr (PERIODIC)
      {
        auto w = periodic_wrap (x);

        ibin = w.ibin;
        frac = w.frac;
      }
    else
      {
        float fbin = (x + float (RANGE)) * float (BINS / (RANGE * 2.0));

        if (fbin < 0)
          return f_[0];
        if (fbin >= BINS)
          return f_[BINS];

        ibin = (int) fbin;
        frac = fbin - ibin;
      }
    return f_[ibin] + frac * (f_[ibin + 1] - f_[ibin]);
  }

  float
  F (float x) const
  {
    Wrap wrap;
    int ibin;
    float frac;
    if constexpr (PERIODIC)
      {
        wrap = periodic_wrap (x);

        ibin = wrap.ibin;
        frac = wrap.frac;
      }
    else
      {
        float fbin = (x + float (RANGE)) * float (BINS / (RANGE * 2.0));

        if (fbin < 0)
          return F_[0] + (x + RANGE) * f_[0];
        if (fbin >= BINS)
          return F_[BINS] + (x - RANGE) * f_[BINS];

        ibin = (int) fbin;
        frac = fbin - ibin;
      }

    float dx = 2.0 * RANGE / BINS;

    float y = F_[ibin] + dx * frac * (f_[ibin] + 0.5f * frac * (f_[ibin+1] - f_[ibin]));

    if constexpr (PERIODIC)
      y += wrap.periods * period_integral_;
    return y;
  }

};
