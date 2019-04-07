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

double
R1_gamma(double p, double T)
{
  double pi  = p / R1_pStar;
  double tau = R1_TStar / T;

  double gamma = 0.0;
  for (int i = 0; i < 34; i++)
    gamma += R1Coef[i][2] * std::pow(7.1 - pi, R1Coef[i][0]) * std::pow(tau - 1.222, R1Coef[i][1]);

  return gamma;
}

double
R1_gamma_pi(double p, double T)
{
  double pi  = p / R1_pStar;
  double tau = R1_TStar / T;

  double sum = 0.0;
  for (int i = 0; i < 34; i++)
    sum -= R1Coef[i][2] * R1Coef[i][0] * std::pow(7.1 - pi, R1Coef[i][0] - 1.0)
            * std::pow(tau - 1.222, R1Coef[i][1]);

  return sum;
}

double
R1_gamma_tau(double p, double T)
{
  double pi  = p / R1_pStar;
  double tau = R1_TStar / T;

  double sum = 0.0;
  for (int i = 0; i < 34; i++)
    sum += R1Coef[i][2] * std::pow(7.1 - pi, R1Coef[i][0]) * R1Coef[i][1]
            * std::pow(tau - 1.222, R1Coef[i][1] - 1.0);

  return sum;
}

double R1_gamma_pi_pi(double p, double T)
{
  double pi  = p / R1_pStar;
  double tau = R1_TStar / T;

  double sum = 0.0;
  for (int i = 0; i < 34; i++)
    sum += R1Coef[i][2] * R1Coef[i][0] * (R1Coef[i][0] - 1.0) * std::pow(7.1 - pi, R1Coef[i][0] - 2.0)
            * std::pow(tau - 1.222, R1Coef[i][1]);

  return sum;
}

double R1_gamma_tau_tau(double p, double T)
{
  double pi  = p / R1_pStar;
  double tau = R1_TStar / T;

  double sum = 0.0;
  for (int i = 0; i < 34; i++)
    sum += R1Coef[i][2] * std::pow(7.1 - pi, R1Coef[i][0]) * R1Coef[i][1] * (R1Coef[i][1] - 1.0)
            * std::pow(tau - 1.222, R1Coef[i][1] - 2.0);

  return sum;
}

double
R1_gamma_pi_tau(double p, double T)
{
  double pi  = p / R1_pStar;
  double tau = R1_TStar / T;

  double sum = 0.0;
  for (int i = 0; i < 34; i++)
    sum -= R1Coef[i][2] * R1Coef[i][0] * std::pow(7.1 - pi, R1Coef[i][0] - 1.0)
            * R1Coef[i][1] * std::pow(tau - 1.222, R1Coef[i][1] - 1.0);

  return sum;
}

double R1_specific_volume(double p, double T)
{
  double pi  = p / R1_pStar;
  return pi * R1_gamma_pi(p, T) * Rgas * T / p;
}

double R1_specific_int_energy(double p, double T)
{
  double pi  = p / R1_pStar;
  double tau = R1_TStar / T;

  return (tau * R1_gamma_tau(p, T) - pi * R1_gamma_pi(p, T)) * Rgas * T;
}

double R1_specific_entropy(double p, double T)
{
  double tau = R1_TStar / T;
  return (tau * R1_gamma_tau(p, T) - R1_gamma(p, T)) * Rgas;
}

double R1_specific_enthalpy(double p, double T)
{
  double tau = R1_TStar / T;
  return tau * R1_gamma_tau(p, T) * Rgas * T;
}

double R1_cp(double p, double T)
{
  double tau = R1_TStar / T;
  return -tau * tau * R1_gamma_tau_tau(p, T) * Rgas;
}

double R1_cv(double p, double T)
{
  double tau = R1_TStar / T;
  double val = R1_gamma_pi(p, T) - tau * R1_gamma_pi_tau(p, T);
  double val2 = -tau * tau * R1_gamma_tau_tau(p, T) + val * val / R1_gamma_pi_pi(p, T);
  return val2 * Rgas;
}

double R1_sound_speed(double p, double T)
{
  double tau = R1_TStar / T;
  double gamma_pi = R1_gamma_pi(p, T);

  double val = gamma_pi - tau * R1_gamma_pi_tau(p, T);
  double val2 = gamma_pi * gamma_pi / (val * val / (tau * tau * R1_gamma_tau_tau(p, T)) - R1_gamma_pi_pi(p, T));
  return std::sqrt(val2 * Rgas * T);
}

double R1_T_from_p_h(double p, double h)
{
  double pi = p / 1.0e6;
  double eta = h / 2500e3;

  double sum = 0.0;
  for (int i = 0; i < 20; i++)
    sum += R1_ph_Coef[i][2] * std::pow(pi, R1_ph_Coef[i][0]) * std::pow(eta + 1.0, R1_ph_Coef[i][1]);

  return sum;
}

double R1_T_from_p_s(double p, double s)
{
  double pi = p / 1.0e6;
  double sigma = s / 1.0e3;

  double sum = 0.0;
  for (int i = 0; i < 20; i++)
    sum += R1_ps_Coef[i][2] * std::pow(pi, R1_ps_Coef[i][0]) * std::pow(sigma + 2.0, R1_ps_Coef[i][1]);

  return sum;
}
