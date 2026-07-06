#include <functional>

template<int RANGE, int BINS>
class ADAATable
{
  std::array<float, BINS + 1> f_;
  std::array<float, BINS + 1> F_;

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

    /* shift antiderivative so that smallest elements are close to zero
     *   -> better resolution for small floats
     */
    auto min_F = *std::min_element (F.begin(), F.end());
    for (size_t i = 0; i < BINS + 1; i++)
      F_[i] = F[i] - min_F;
  }

  float
  f (float x)
  {
    float fbin = (x + float (RANGE)) * float (BINS / (RANGE * 2.0));

    if (fbin < 0)
      return f_[0];
    if (fbin >= BINS)
      return f_[BINS];

    int ibin = (int) fbin;

    float frac = fbin - ibin;
    return f_[ibin] + frac * (f_[ibin + 1] - f_[ibin]);
  }

  float
  F (float x)
  {
    float fbin = (x + float (RANGE)) * float (BINS / (RANGE * 2.0));

    if (fbin < 0)
      return F_[0] + (x + RANGE) * f_[0];
    if (fbin >= BINS)
      return F_[BINS] + (x - RANGE) * f_[BINS];

    int ibin = (int) fbin;

    float frac = fbin - ibin;
    float dx = frac * 2.f * RANGE / BINS;
    return F_[ibin] + dx * f_[ibin] + 0.5f * dx * dx * (f_[ibin+1] - f_[ibin]) / (2.f * RANGE / BINS);;
  }

};
