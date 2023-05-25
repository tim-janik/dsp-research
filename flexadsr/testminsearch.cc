// This Source Code Form is licensed MPL-2.0: http://mozilla.org/MPL/2.0

#include <cmath>

#include "minsearch.hh"

// Rosenbrock function, minimum at (x,y) = (a, a^2)
double rosenbrock (double x, double y)
{
  const double a = 1;
  const double b = 100;

  return std::pow (a - x, 2) + b * std::pow (y - std::pow (x, 2), 2);
}

int
main()
{
  int n = 0;
  auto score_func = [&] (std::vector<double>& v) {
    n++;
    return rosenbrock (v[0], v[1]);
  };

  std::vector<double> v(2);
  MinSearch::minimize (v, score_func);
  printf ("%f %f, n=%d function evaluations\n", v[0], v[1], n);
}
