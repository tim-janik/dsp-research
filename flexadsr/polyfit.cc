#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <cctype>
#include <cmath>

#include <vector>
#include <algorithm>

#include "minsearch.hh"

using std::vector;

int
main (int argc, char **argv)
{
  struct P {
    int slope;
    double p;
  };
  vector<P> as, bs;

  char buffer[1024];

  while (fgets (buffer, 1024, stdin))
    {
      if (isdigit (buffer[0]) || buffer[0] == '-')
        {
          char *sls = strtok (buffer, " \n");
          int slope = atoi (sls);

          float a = atof (strtok (NULL, " \n"));
          float b = atof (strtok (NULL, " \n"));
          as.push_back ({ slope, a });
          bs.push_back ({ slope, b });
        }
    }
  int k = 0;
  for (auto ps : { as, bs })
    {
      auto score_func = [ps](const vector<double>& v) {
        double score = 0;
        for (auto p : ps)
          {
            double slope = p.slope / 100.;
            double est = 0;
            for (size_t i = 0; i < v.size(); i++)
              est += pow (slope, i + 1) * v[i];
            score += pow (est - p.p, 2);
          }
        return score;
      };
      vector<double> v (atoi (argv[1]));

      MinSearch::minimize (v, score_func, k ? "#B" : "#A");
      for (int i = -100; i <= 100; i++)
        {
          double slope = i / 100.;
          double est = 0;
          for (size_t i = 0; i < v.size(); i++)
            est += pow (slope, i + 1) * v[i];
          //printf ("%d %f #%c\n", i, est, "AB"[k]);
        }
      for (size_t i = 0; i < v.size(); i++)
        {
          printf ("  params.%c += pow (slope, %zd) * %.17g;\n", "ab"[k], i + 1, v[i]);
        }
      k++;
    }
}
