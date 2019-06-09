#include "ThermalConductivity.h"
#include "Viscosity.h"

double labmda0_bar(double T_bar)
{
  double T_bar_sq = T_bar * T_bar;
  double denom = L0[0] + L0[1] / T_bar + L0[2] / T_bar_sq + L0[3] / (T_bar * T_bar_sq) + L0[4] / (T_bar_sq * T_bar_sq);

  return sqrt(T_bar) / denom;
}

double labmda1_bar(double rho_bar, double T_bar)
{
  double sum = 0.0;
  for (int i = 0; i < 5; i++)
  {
    double secSum = 0.0;
    for (int j = 0; j < 6; j++)
      secSum += L1[i][j] * pow(rho_bar - 1.0, j);
    sum += pow(1.0 / T_bar - 1.0, i) * secSum;
  }

  return exp(rho_bar * sum);
}

double labmda2_bar(double rho_bar, double T_bar, double cp, double cv, double mu_bar)
{
  /*
  double y = correlation_length_TC(rho_bar, T_bar) / 0.40;
  double zy_val = Zy(y, rho_bar, cp / cv);

  //std::cout << "y = " << y << std::endl;
  //std::cout << "zy_val = " << zy_val << std::endl;

  return 177.8514 * rho_bar * T_bar * cp * zy_val / (0.46151805e3 * mu_bar);*/
  return 0.0;
}

double zeta_R1(double p, double T)
{
  return P_CRIT / RHO_CRIT * R1_drho_dp(p, T);
}

double zeta_R2(double p, double T)
{
  return P_CRIT / RHO_CRIT * R2_drho_dp(p, T);
}

double zeta_R3(double rho_bar, double T_bar)
{
  return P_CRIT / R3_dp_ddelta(rho_bar, 1.0 / T_bar);
}

double zeta_REF(double rho_bar)
{
  int j = -1;
  if      (rho_bar <= 0.310559006)    j = 0;
  else if (rho_bar <= 0.776397516)    j = 1;
  else if (rho_bar <= 1.242236025)    j = 2;
  else if (rho_bar <= 1.863354037)    j = 3;
  else                                j = 4;

  double denom = 0.0;
  for (int i = 0; i < 6; i++)
    denom += A[i][j] * pow(rho_bar, i);

  return 1.0 / denom;
}

double correlation_length_TC(double rho_bar, double T_bar, double zeta)
{
  double dchi_bar = rho_bar * (zeta - zeta_REF(rho_bar) * 1.5 / T_bar);
  dchi_bar = (dchi_bar > 0.0) ? dchi_bar : 0.0;

  return 0.13 * std::pow(dchi_bar / 0.06, 0.630 / 1.239);
}

double Zy(double y, double rho_bar, double kappa)
{
  if (y < 1.2e-7)
    return 0.0;
  else
  {
    double val_1 = (1.0 - 1.0 / kappa) * atan(y) + y / kappa;
    double val_2 = 1.0 - exp(-1.0 / (1.0 / y + y * y / (3.0 * rho_bar * rho_bar)));

    return 2.0 / (PI * y) * (val_1 - val_2);
  }
}

double thermal_conductivity_no_enhancement(double rho, double T)
{
  double rho_bar = rho / 322.0;
  double T_bar = T / 647.096;

  return labmda0_bar(T_bar) * labmda1_bar(rho_bar, T_bar) * 1.0e-3;
}

double thermal_conductivity_R1(double p, double T)
{
  double rho_bar = 1.0 / R1_specific_volume(p, T) / 322.0;
  double T_bar = T / 647.096;
  double cp = R1_cp(p, T);
  double cv = R1_cv(p, T);

  double lambda_bar_part1 = labmda0_bar(T_bar) * labmda1_bar(rho_bar, T_bar);
  double zeta = zeta_R1(p, T);
  double xi = correlation_length_TC(rho_bar, T_bar, zeta);
  double mu_bar = mu0_bar(T_bar) * mu1_bar(rho_bar, T_bar);

  double y = xi / 0.40;
  double zy_val = Zy(y, rho_bar, cp / cv);
  double lambda2_bar = 177.8514 * rho_bar * T_bar * cp * zy_val / (0.46151805e3 * mu_bar);

  return (lambda_bar_part1 + lambda2_bar) * 1.0e-3;
}

double thermal_conductivity_R2(double p, double T)
{
  double rho_bar = 1.0 / R2_specific_volume(p, T) / 322.0;
  double T_bar = T / 647.096;
  double cp = R2_cp(p, T);
  double cv = R2_cv(p, T);

  double lambda_bar_part1 = labmda0_bar(T_bar) * labmda1_bar(rho_bar, T_bar);
  double zeta = zeta_R2(p, T); // change to R2
  double xi = correlation_length_TC(rho_bar, T_bar, zeta);
  double mu_bar = mu0_bar(T_bar) * mu1_bar(rho_bar, T_bar);

  double y = xi / 0.40;
  double zy_val = Zy(y, rho_bar, cp / cv);
  double lambda2_bar = 177.8514 * rho_bar * T_bar * cp * zy_val / (0.46151805e3 * mu_bar);

  return (lambda_bar_part1 + lambda2_bar) * 1.0e-3;
}

double thermal_conductivity_R3(double rho, double T)
{
  double rho_bar = rho / 322.0;
  double T_bar = T / 647.096;

  double lambda_bar_part1 = labmda0_bar(T_bar) * labmda1_bar(rho_bar, T_bar);

  double cp = R3_cp(rho, T);
  double cv = R3_cv(rho, T);
  double zeta_val = zeta_R3(rho_bar, T_bar);
  double xi = correlation_length_TC(rho_bar, T_bar, zeta_val);
  double mu_bar = mu0_bar(T_bar) * mu1_bar(rho_bar, T_bar);

  double y = xi / 0.40;
  double zy_val = Zy(y, rho_bar, cp / cv);
  double lambda2_bar = 177.8514 * rho_bar * T_bar * cp * zy_val / (0.46151805e3 * mu_bar);

  return (lambda_bar_part1 + lambda2_bar) * 1.0e-3;
}
