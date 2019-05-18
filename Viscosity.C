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

double viscosity(double rho, double T, bool critical_enhancement)
{
  double rho_bar = rho / 322.0;
  double T_bar = T / 647.096;

  double mu = mu0_bar(T_bar) * mu1_bar(rho_bar, T_bar) * 1.0e-6;

  if (!critical_enhancement)
    return mu;
  else
    return mu * mu2_bar(rho_bar, T_bar);
}

double zeta(double rho_bar, double T_bar)
{
  return Pcrit / R3_dp_ddelta(rho_bar, 1.0 / T_bar);
}

double correlation_length(double rho_bar, double T_bar)
{
  double dchi_bar = std::max(rho_bar * (zeta(rho_bar, T_bar) - zeta(rho_bar, 1.5) * 1.5 / T_bar), 0.0);

  return 0.13 * std::pow(dchi_bar / 0.06, 0.630 / 1.239);
}

double mu2_bar(double rho_bar, double T_bar)
{
  double Y = 0.0;
  double xi = correlation_length(rho_bar, T_bar);

  double qC_xi = xi / 1.9;
  double qD_xi = xi / 1.1;

  if (xi < 0.3817016416)
    Y = 0.2 * qC_xi * std::pow(qD_xi, 5) * (1.0 - qC_xi + qC_xi * qC_xi - 765.0 / 504.0 * qD_xi * qD_xi);
  else
  {
    double psi_D = acos(1.0 / sqrt(1.0 + qD_xi * qD_xi));
    double w = sqrt(fabs((qC_xi - 1.0) / (qC_xi + 1.0))) * tan(0.5 * psi_D);
    double Lw = (qC_xi > 1.0) ? log((1.0 + w) / (1.0 - w)) : 2.0 * atan(fabs(w));
    double qC_xi2 = qC_xi * qC_xi;

    Y = sin(3.0 * psi_D) / 12.0 - 0.25 / qC_xi * sin(2.0 * psi_D) + 1.0 / qC_xi2 * (1.0 - 1.25 * qC_xi2) * sin(psi_D)
        - 1.0 / (qC_xi * qC_xi2) * ((1.0 - 1.5 * qC_xi2) * psi_D - pow(fabs(qC_xi2 - 1.0), 1.5) * Lw);
  }

  return std::exp(0.068 * Y);
}
