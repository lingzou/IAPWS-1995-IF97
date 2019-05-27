#include <cmath>
#include <iostream>
#include "Viscosity.h"

double mu0_bar(double T_bar)
{
  double denom = MU0_H[0] + MU0_H[1] / T_bar + MU0_H[2] / (T_bar * T_bar) + MU0_H[3] / (T_bar * T_bar * T_bar);

  return 100.0 * std::sqrt(T_bar) / denom;
}

double mu1_bar(double rho_bar, double T_bar)
{
  double SecTerm[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  for (int k = 0; k < 21; k++)
  {
    // MU1_POS[k][0] -> i
    // MU1_POS[k][1] -> j
    // MU1_H[k] -> H_ij
    SecTerm[MU1_POS[k][0]] += MU1_H[k] * std::pow(rho_bar - 1.0, MU1_POS[k][1]);
  }

  double sum = 0.0;
  for (int i = 0; i < 6; i++)
    sum += std::pow(1.0 / T_bar - 1.0, i) * SecTerm[i];

  return std::exp(rho_bar * sum);
}

double viscosity(double rho, double T)
{
  double rho_bar = rho / 322.0;
  double T_bar = T / 647.096;

  return mu0_bar(T_bar) * mu1_bar(rho_bar, T_bar) * 1.0e-6;
}
