#include "minsearch.hh"

using std::vector;

static constexpr double rate = 192000;

double
distq (double ax, double ay, double bx, double by)
{
  return (ax - bx) * (ax - bx) + (ay - by) * (ay - by);
}

void
debug_adsr_curve (int i, double a, double b)
{
  double slope = i / 100.;
  double t1y = 0.5 + 0.25 * slope;
  double c = 1 - t1y * t1y * a - t1y * b;
  double x = 0;
  double y = 0;

  for (int j = 0; j < rate; j++)
    {
      x += 1 / rate;
      y += (y * y * a + y * b + c) / rate;
      printf ("%d %f %f\n", i, x, y);
    }
}

/*
 * find optimal a, b, c for a given ADSR slope, such that y can be recursively
 * updated in the ADSR using
 *
 *   y += (a * y^2 + b * y + c)
 *
 */
int
main()
{
  vector<double> vec (2);
  for (int i = 0; i <= 100; i += 1)
    {
      double slope = i / 100.;
      double t1x = 0.5 - 0.25 * slope;
      double t1y = 0.5 + 0.25 * slope;

      printf ("#%f %f\n", t1x, t1y);

      auto score_func = [&](const vector<double>& vec) {
        double y = 0;
        double t1d = 1e9;
        double t2d = 1e9;
        double a = vec[0];
        double b = vec[1];
        double c = 1 - t1y * t1y * a - t1y * b;
        double x = 0;
        for (int i = 0; i < rate + 20; i++)
          {
            // ADSR curve should cross this point: (t1x, t1y)
            t1d = std::min (t1d, distq (x, y, t1x, t1y));

            // ADSR curve should reach point (1, 1)
            t2d = std::min (t2d, distq (x, y, 1, 1));

            x += 1 / rate;
            y += (y * y * a + y * b + c) / rate;
          }
        return t1d + t2d;
      };

      char str[20];
      sprintf (str, "S %.3f", slope);

      MinSearch::minimize (vec, score_func, str);

      printf("%d %f %f %.10f\n", i, vec[0], vec[1], score_func (vec));
      fflush (stdout);

      //debug_adsr_curve (i, vec[0], vec[1]);
  }
}
