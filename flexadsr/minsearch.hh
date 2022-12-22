// This Source Code Form is licensed MPL-2.0: http://mozilla.org/MPL/2.0

#ifndef __MIN_SEARCH_HH__
#define __MIN_SEARCH_HH__

#include <vector>
#include <string>

namespace MinSearch
{

/*
 * Search a local minimum for N-dimensional vector according to a score
 * function using gradient descent algorithm.
 */
template<class ScoreFunc>
void
minimize (std::vector<double>& vec, ScoreFunc score_func, const std::string& label = "")
{
  int i = 0;
  int maxi = 50;

  double old_score = score_func (vec);
  double learning_rate = 1;
  bool trace = !label.empty();
  while (i < maxi)
    {
      std::vector<double> d;

      if (trace)
        fprintf (stderr, "[ ");

      for (size_t m = 0; m < vec.size(); m++)
        {
          double h = 1e-9;
          std::vector<double> vec1 = vec;
          vec1[m] += h;
          double deriv = (score_func (vec1) - old_score) / h;
          d.push_back (-deriv * learning_rate);

          if (trace)
            fprintf (stderr, "% 3.7f ", deriv);
        }

      if (trace)
          fprintf (stderr, "]  %25.17g %s           \r", old_score, label.c_str());

      std::vector<double> vec1 = vec;
      for (size_t m = 0; m < vec.size(); m++)
        vec1[m] += d[m];

      double new_score = score_func (vec1);
      if (new_score < old_score)
        {
          vec = vec1;
          old_score = new_score;

          learning_rate *= 1.5;
          i = 0;
        }
      else
        {
          learning_rate /= 1.5;
        }
      i++;
    }
  if (trace)
    fprintf (stderr, "\n");
}

}

#endif
