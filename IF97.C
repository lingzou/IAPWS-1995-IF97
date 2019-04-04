#include <cmath>

#include "IF97.h"

void a_func() {}

double
B23_p_from_T(double T)
{
  double theta = T;
  return (B23_n[0] + B23_n[1] * theta + B23_n[2] * theta * theta) * 1.0e6;
}

double
B23_T_from_p(double p)
{
  double pi = p * 1.0e-6;
  return B23_n[3] + std::sqrt((pi - B23_n[4])/B23_n[2]);
}
