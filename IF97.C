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
R1_gamma(double pi, double tau)
{
  double gamma = 0.0;
  for (int i = 0; i < 34; i++)
    gamma += R1Coef[i][2] * std::pow(7.1 - pi, R1Coef[i][0]) * std::pow(tau - 1.222, R1Coef[i][1]);

  return gamma;
}

double
R1_gamma_pi(double pi, double tau)
{
  double sum = 0.0;
  for (int i = 0; i < 34; i++)
    sum -= R1Coef[i][2] * R1Coef[i][0] * std::pow(7.1 - pi, R1Coef[i][0] - 1.0)
            * std::pow(tau - 1.222, R1Coef[i][1]);

  return sum;
}

double
R1_gamma_tau(double pi, double tau)
{
  double sum = 0.0;
  for (int i = 0; i < 34; i++)
    sum += R1Coef[i][2] * std::pow(7.1 - pi, R1Coef[i][0]) * R1Coef[i][1]
            * std::pow(tau - 1.222, R1Coef[i][1] - 1.0);

  return sum;
}

double R1_gamma_pi_pi(double pi, double tau)
{
  double sum = 0.0;
  for (int i = 0; i < 34; i++)
    sum += R1Coef[i][2] * R1Coef[i][0] * (R1Coef[i][0] - 1.0) * std::pow(7.1 - pi, R1Coef[i][0] - 2.0)
            * std::pow(tau - 1.222, R1Coef[i][1]);

  return sum;
}

double R1_gamma_tau_tau(double pi, double tau)
{
  double sum = 0.0;
  for (int i = 0; i < 34; i++)
    sum += R1Coef[i][2] * std::pow(7.1 - pi, R1Coef[i][0]) * R1Coef[i][1] * (R1Coef[i][1] - 1.0)
            * std::pow(tau - 1.222, R1Coef[i][1] - 2.0);

  return sum;
}

double
R1_gamma_pi_tau(double pi, double tau)
{
  double sum = 0.0;
  for (int i = 0; i < 34; i++)
    sum -= R1Coef[i][2] * R1Coef[i][0] * std::pow(7.1 - pi, R1Coef[i][0] - 1.0)
            * R1Coef[i][1] * std::pow(tau - 1.222, R1Coef[i][1] - 1.0);

  return sum;
}

double R1_specific_volume(double p, double T)
{
  double pi  = p / R1_pStar;
  double tau = R1_TStar / T;

  return pi * R1_gamma_pi(pi, tau) * Rgas * T / p;
}

double R1_specific_int_energy(double p, double T)
{
  double pi  = p / R1_pStar;
  double tau = R1_TStar / T;

  return (tau * R1_gamma_tau(pi, tau) - pi * R1_gamma_pi(pi, tau)) * Rgas * T;
}

double R1_specific_entropy(double p, double T)
{
  double pi  = p / R1_pStar;
  double tau = R1_TStar / T;

  return (tau * R1_gamma_tau(pi, tau) - R1_gamma(pi, tau)) * Rgas;
}

double R1_specific_enthalpy(double p, double T)
{
  double pi  = p / R1_pStar;
  double tau = R1_TStar / T;

  return tau * R1_gamma_tau(pi, tau) * Rgas * T;
}

double R1_cp(double p, double T)
{
  double pi  = p / R1_pStar;
  double tau = R1_TStar / T;

  return -tau * tau * R1_gamma_tau_tau(pi, tau) * Rgas;
}

double R1_cv(double p, double T)
{
  double pi  = p / R1_pStar;
  double tau = R1_TStar / T;

  double val = R1_gamma_pi(pi, tau) - tau * R1_gamma_pi_tau(pi, tau);
  double val2 = -tau * tau * R1_gamma_tau_tau(pi, tau) + val * val / R1_gamma_pi_pi(pi, tau);
  return val2 * Rgas;
}

double R1_sound_speed(double p, double T)
{
  double pi  = p / R1_pStar;
  double tau = R1_TStar / T;

  double gamma_pi = R1_gamma_pi(pi, tau);

  double val = gamma_pi - tau * R1_gamma_pi_tau(pi, tau);
  double val2 = gamma_pi * gamma_pi / (val * val / (tau * tau * R1_gamma_tau_tau(pi, tau)) - R1_gamma_pi_pi(pi, tau));
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

double R2_gamma_0(double pi, double tau)
{
  double gamma_0 = std::log(pi);
  for (int i = 0; i < 9; i++)
    gamma_0 += R2Coef0[i][1] * std::pow(tau, R2Coef0[i][0]);

  return gamma_0;
}

double R2_gamma_r(double pi, double tau)
{
  double gamma_r = 0;
  for (int i = 0; i < 43; i++)
    gamma_r += R2Coefr[i][2] * std::pow(pi, R2Coefr[i][0]) * std::pow(tau - 0.5, R2Coefr[i][1]);

  return gamma_r;
}

double R2_gamma_0_pi(double pi, double /*tau*/)
{
  return 1.0 / pi;
}

double R2_gamma_0_pi_pi(double pi, double tau)
{
  return -1.0 / (pi * pi);
}
double R2_gamma_0_tau(double pi, double tau)
{
  double der = 0.0;
  for (int i = 0; i < 9; i++)
    der += R2Coef0[i][1] * R2Coef0[i][0] * std::pow(tau, R2Coef0[i][0] - 1.0);

  return der;
}

double R2_gamma_0_tau_tau(double pi, double tau)
{
  double der = 0.0;
  for (int i = 0; i < 9; i++)
    der += R2Coef0[i][1] * R2Coef0[i][0] * (R2Coef0[i][0] - 1.0) * std::pow(tau, R2Coef0[i][0] - 2.0);

  return der;
}

double R2_gamma_0_pi_tau(double /*pi*/, double /*tau*/)
{
  return 0.0;
}

double R2_gamma_r_pi(double pi, double tau)
{
  double der = 0.0;
  for (int i = 0; i < 43; i++)
    der += R2Coefr[i][2] * R2Coefr[i][0] * std::pow(pi, R2Coefr[i][0] - 1.0)
            * std::pow(tau - 0.5, R2Coefr[i][1]);

  return der;
}

double R2_gamma_r_pi_pi(double pi, double tau)
{
  double der = 0.0;
  for (int i = 0; i < 43; i++)
    der += R2Coefr[i][2] * R2Coefr[i][0] * (R2Coefr[i][0] - 1.0 )* std::pow(pi, R2Coefr[i][0] - 2.0)
            * std::pow(tau - 0.5, R2Coefr[i][1]);

  return der;
}

double R2_gamma_r_tau(double pi, double tau)
{
  double der = 0.0;
  for (int i = 0; i < 43; i++)
    der += R2Coefr[i][2] * std::pow(pi, R2Coefr[i][0]) * R2Coefr[i][1]
            * std::pow(tau - 0.5, R2Coefr[i][1] - 1.0);

  return der;
}
double R2_gamma_r_tau_tau(double pi, double tau)
{
  double der = 0.0;
  for (int i = 0; i < 43; i++)
    der += R2Coefr[i][2] * std::pow(pi, R2Coefr[i][0]) * R2Coefr[i][1] * (R2Coefr[i][1] - 1.0)
            * std::pow(tau - 0.5, R2Coefr[i][1] - 2.0);

  return der;
}

double R2_gamma_r_pi_tau(double pi, double tau)
{
  double der = 0.0;
  for (int i = 0; i < 43; i++)
    der += R2Coefr[i][2] * R2Coefr[i][0] * std::pow(pi, R2Coefr[i][0] - 1.0)
            * R2Coefr[i][1] * std::pow(tau - 0.5, R2Coefr[i][1] - 1.0);

  return der;
}

double R2_specific_volume(double p, double T)
{
  double pi = p / R2_pStar;
  double tau = R2_TStar / T;

  return pi * (R2_gamma_0_pi(pi, tau) + R2_gamma_r_pi(pi, tau)) * Rgas * T / p;
}

double R2_specific_int_energy(double p, double T)
{
  double pi = p / R2_pStar;
  double tau = R2_TStar / T;

  double val = tau * (R2_gamma_0_tau(pi, tau) + R2_gamma_r_tau(pi, tau))
              - pi * (R2_gamma_0_pi(pi, tau)  + R2_gamma_r_pi(pi, tau));

  return val * Rgas * T;
}

double R2_specific_entropy(double p, double T)
{
  double pi = p / R2_pStar;
  double tau = R2_TStar / T;

  double val = tau * (R2_gamma_0_tau(pi, tau) + R2_gamma_r_tau(pi, tau))
              - (R2_gamma_0(pi, tau)  + R2_gamma_r(pi, tau));

  return val * Rgas;
}

double R2_specific_enthalpy(double p, double T)
{
  double pi = p / R2_pStar;
  double tau = R2_TStar / T;

  return tau * (R2_gamma_0_tau(pi, tau) + R2_gamma_r_tau(pi, tau)) * Rgas * T;
}

double R2_cp(double p, double T)
{
  double pi = p / R2_pStar;
  double tau = R2_TStar / T;

  return -tau * tau * (R2_gamma_0_tau_tau(pi, tau) + R2_gamma_r_tau_tau(pi, tau)) * Rgas;
}

double R2_cv(double p, double T)
{
  double pi = p / R2_pStar;
  double tau = R2_TStar / T;

  double cp = -tau * tau * (R2_gamma_0_tau_tau(pi, tau) + R2_gamma_r_tau_tau(pi, tau)) * Rgas;
  double val1 = 1.0 + pi * R2_gamma_r_pi(pi, tau) - tau * pi * R2_gamma_r_pi_tau(pi, tau);
  double val2 = 1.0 - pi * pi * R2_gamma_r_pi_pi(pi, tau);

  return cp - (val1 * val1 / val2) * Rgas;
}

double R2_sound_speed(double p, double T)
{
  double pi = p / R2_pStar;
  double tau = R2_TStar / T;

  double gamma_r_pi = R2_gamma_r_pi(pi, tau);
  double val1 = 1.0 + 2.0 * pi * R2_gamma_r_pi(pi, tau) + pi * pi * gamma_r_pi * gamma_r_pi;
  double val2 = 1.0 - pi * pi * R2_gamma_r_pi_pi(pi, tau);
  double val3 = 1.0 + pi * R2_gamma_r_pi(pi, tau) - tau * pi * R2_gamma_r_pi_tau(pi, tau);
  double val4 = tau * tau * (R2_gamma_0_tau_tau(pi, tau) + R2_gamma_r_tau_tau(pi, tau));

  double w2 = val1 / (val2 + val3 * val3 / val4) * Rgas * T;

  return std::sqrt(w2);
}

// R2 Meta
double R2Meta_gamma_0(double pi, double tau)
{
  double gamma_0 = std::log(pi);
  for (int i = 0; i < 9; i++)
    gamma_0 += R2MetaCoef0[i][1] * std::pow(tau, R2MetaCoef0[i][0]);

  return gamma_0;
}

double R2Meta_gamma_r(double pi, double tau)
{
  double gamma_r = 0;
  for (int i = 0; i < 13; i++)
    gamma_r += R2MetaCoefr[i][2] * std::pow(pi, R2MetaCoefr[i][0]) * std::pow(tau - 0.5, R2MetaCoefr[i][1]);

  return gamma_r;
}

double R2Meta_gamma_0_pi(double pi, double /*tau*/)
{
  return 1.0 / pi;
}

double R2Meta_gamma_0_pi_pi(double pi, double tau)
{
  return -1.0 / (pi * pi);
}
double R2Meta_gamma_0_tau(double pi, double tau)
{
  double der = 0.0;
  for (int i = 0; i < 9; i++)
    der += R2MetaCoef0[i][1] * R2MetaCoef0[i][0] * std::pow(tau, R2MetaCoef0[i][0] - 1.0);

  return der;
}

double R2Meta_gamma_0_tau_tau(double pi, double tau)
{
  double der = 0.0;
  for (int i = 0; i < 9; i++)
    der += R2MetaCoef0[i][1] * R2MetaCoef0[i][0] * (R2MetaCoef0[i][0] - 1.0) * std::pow(tau, R2MetaCoef0[i][0] - 2.0);

  return der;
}

double R2Meta_gamma_0_pi_tau(double /*pi*/, double /*tau*/)
{
  return 0.0;
}

double R2Meta_gamma_r_pi(double pi, double tau)
{
  double der = 0.0;
  for (int i = 0; i < 13; i++)
    der += R2MetaCoefr[i][2] * R2MetaCoefr[i][0] * std::pow(pi, R2MetaCoefr[i][0] - 1.0)
            * std::pow(tau - 0.5, R2MetaCoefr[i][1]);

  return der;
}

double R2Meta_gamma_r_pi_pi(double pi, double tau)
{
  double der = 0.0;
  for (int i = 0; i < 13; i++)
    der += R2MetaCoefr[i][2] * R2MetaCoefr[i][0] * (R2MetaCoefr[i][0] - 1.0 )* std::pow(pi, R2MetaCoefr[i][0] - 2.0)
            * std::pow(tau - 0.5, R2MetaCoefr[i][1]);

  return der;
}

double R2Meta_gamma_r_tau(double pi, double tau)
{
  double der = 0.0;
  for (int i = 0; i < 13; i++)
    der += R2MetaCoefr[i][2] * std::pow(pi, R2MetaCoefr[i][0]) * R2MetaCoefr[i][1]
            * std::pow(tau - 0.5, R2MetaCoefr[i][1] - 1.0);

  return der;
}
double R2Meta_gamma_r_tau_tau(double pi, double tau)
{
  double der = 0.0;
  for (int i = 0; i < 13; i++)
    der += R2MetaCoefr[i][2] * std::pow(pi, R2MetaCoefr[i][0]) * R2MetaCoefr[i][1] * (R2MetaCoefr[i][1] - 1.0)
            * std::pow(tau - 0.5, R2MetaCoefr[i][1] - 2.0);

  return der;
}

double R2Meta_gamma_r_pi_tau(double pi, double tau)
{
  double der = 0.0;
  for (int i = 0; i < 13; i++)
    der += R2MetaCoefr[i][2] * R2MetaCoefr[i][0] * std::pow(pi, R2MetaCoefr[i][0] - 1.0)
            * R2MetaCoefr[i][1] * std::pow(tau - 0.5, R2MetaCoefr[i][1] - 1.0);

  return der;
}

double R2Meta_specific_volume(double p, double T)
{
  double pi = p / R2_pStar;
  double tau = R2_TStar / T;

  return pi * (R2Meta_gamma_0_pi(pi, tau) + R2Meta_gamma_r_pi(pi, tau)) * Rgas * T / p;
}

double R2Meta_specific_int_energy(double p, double T)
{
  double pi = p / R2_pStar;
  double tau = R2_TStar / T;

  double val = tau * (R2Meta_gamma_0_tau(pi, tau) + R2Meta_gamma_r_tau(pi, tau))
              - pi * (R2Meta_gamma_0_pi(pi, tau)  + R2Meta_gamma_r_pi(pi, tau));

  return val * Rgas * T;
}

double R2Meta_specific_entropy(double p, double T)
{
  double pi = p / R2_pStar;
  double tau = R2_TStar / T;

  double val = tau * (R2Meta_gamma_0_tau(pi, tau) + R2Meta_gamma_r_tau(pi, tau))
              - (R2Meta_gamma_0(pi, tau)  + R2Meta_gamma_r(pi, tau));

  return val * Rgas;
}

double R2Meta_specific_enthalpy(double p, double T)
{
  double pi = p / R2_pStar;
  double tau = R2_TStar / T;

  return tau * (R2Meta_gamma_0_tau(pi, tau) + R2Meta_gamma_r_tau(pi, tau)) * Rgas * T;
}

double R2Meta_cp(double p, double T)
{
  double pi = p / R2_pStar;
  double tau = R2_TStar / T;

  return -tau * tau * (R2Meta_gamma_0_tau_tau(pi, tau) + R2Meta_gamma_r_tau_tau(pi, tau)) * Rgas;
}

double R2Meta_cv(double p, double T)
{
  double pi = p / R2_pStar;
  double tau = R2_TStar / T;

  double cp = -tau * tau * (R2Meta_gamma_0_tau_tau(pi, tau) * R2Meta_gamma_r_tau_tau(pi, tau)) * Rgas;
  double val1 = 1.0 + pi * R2Meta_gamma_r_pi(pi, tau) - tau * pi * R2Meta_gamma_r_pi_tau(pi, tau);
  double val2 = 1.0 - pi * pi * R2Meta_gamma_r_pi_pi(pi, tau);

  return cp - val1 * val1 / val2;
}

double R2Meta_sound_speed(double p, double T)
{
  double pi = p / R2_pStar;
  double tau = R2_TStar / T;

  double gamma_r_pi = R2Meta_gamma_r_pi(pi, tau);
  double val1 = 1.0 + 2.0 * pi * R2Meta_gamma_r_pi(pi, tau) + pi * pi * gamma_r_pi * gamma_r_pi;
  double val2 = 1.0 - pi * pi * R2Meta_gamma_r_pi_pi(pi, tau);
  double val3 = 1.0 + pi * R2Meta_gamma_r_pi(pi, tau) - tau * pi * R2Meta_gamma_r_pi_tau(pi, tau);
  double val4 = tau * tau * (R2Meta_gamma_0_tau_tau(pi, tau) + R2Meta_gamma_r_tau_tau(pi, tau));

  double w2 = val1 / (val2 + val3 * val3 / val4) * Rgas * T;

  return std::sqrt(w2);
}

double B2bc_p_from_h(double h)
{
  double eta = h * 1.0e-3;

  return (B2bc_n[0] + B2bc_n[1] * eta + B2bc_n[2] * eta * eta) * 1.0e6;
}
double B2bc_h_from_p(double p)
{
  double pi = p * 1.0e-6;

  return (B2bc_n[3] + std::sqrt((pi - B2bc_n[4]) / B2bc_n[2])) * 1.0e3;
}

double R2a_T_from_p_h(double p, double h)
{
  double pi = p * 1.0e-6;
  double eta = h / 2.0e6;

  double theta = 0.0;
  for (int i = 0; i < 34; i++)
    theta += R2a_ph_Coef[i][2] * std::pow(pi, R2a_ph_Coef[i][0])
            * std::pow(eta - 2.1, R2a_ph_Coef[i][1]);

  return theta;  // theta = T
}

double R2b_T_from_p_h(double p, double h)
{
  double pi = p * 1.0e-6;
  double eta = h / 2.0e6;

  double theta = 0.0;
  for (int i = 0; i < 38; i++)
    theta += R2b_ph_Coef[i][2] * std::pow(pi - 2.0, R2b_ph_Coef[i][0])
            * std::pow(eta - 2.6, R2b_ph_Coef[i][1]);

  return theta;  // theta = T
}

double R2c_T_from_p_h(double p, double h)
{
  double pi = p * 1.0e-6;
  double eta = h / 2.0e6;

  double theta = 0.0;
  for (int i = 0; i < 23; i++)
    theta += R2c_ph_Coef[i][2] * std::pow(pi + 25.0, R2c_ph_Coef[i][0])
            * std::pow(eta - 1.8, R2c_ph_Coef[i][1]);

  return theta;  // theta = T
}

double R2a_T_from_p_s(double p, double s)
{
  double pi = p * 1.0e-6;
  double sigma = s / 2.0e3;

  double theta = 0.0;
  for (int i = 0; i < 46; i ++)
    theta += R2a_ps_Coef[i][2] * std::pow(pi, R2a_ps_Coef[i][0]) * std::pow(sigma - 2.0, R2a_ps_Coef[i][1]);

  return theta; // theta = T
}

double R2b_T_from_p_s(double p, double s)
{
  double pi = p * 1.0e-6;
  double sigma = s / 0.7853e3;

  double theta = 0.0;
  for (int i = 0; i < 44; i ++)
    theta += R2b_ps_Coef[i][2] * std::pow(pi, R2b_ps_Coef[i][0]) * std::pow(10.0 - sigma, R2b_ps_Coef[i][1]);

  return theta; // theta = T
}

double R2c_T_from_p_s(double p, double s)
{
  double pi = p * 1.0e-6;
  double sigma = s / 2.9251e3;

  double theta = 0.0;
  for (int i = 0; i < 30; i ++)
    theta += R2c_ps_Coef[i][2] * std::pow(pi, R2c_ps_Coef[i][0]) * std::pow(2.0 - sigma, R2c_ps_Coef[i][1]);

  return theta; // theta = T
}

double R3_phi(double delta, double tau)
{
  double phi = R3Coef[0][2] * std::log(delta);

  for (int i = 1; i < 40; i++)
    phi += R3Coef[i][2] * std::pow(delta, R3Coef[i][0]) * std::pow(tau, R3Coef[i][1]);

  return phi;
}

double R3_phi_delta(double delta, double tau)
{
  double der = R3Coef[0][2] / delta;

  for (int i = 1; i < 40; i++)
    der += R3Coef[i][2] * R3Coef[i][0] * std::pow(delta, R3Coef[i][0] - 1.0)
            * std::pow(tau, R3Coef[i][1]);

  return der;
}

double R3_phi_delta_delta(double delta, double tau)
{
  double der = -R3Coef[0][2] / (delta * delta);

  for (int i = 1; i < 40; i++)
    der += R3Coef[i][2] * R3Coef[i][0] * (R3Coef[i][0] - 1.0) * std::pow(delta, R3Coef[i][0] - 2.0)
            * std::pow(tau, R3Coef[i][1]);

  return der;
}

double R3_phi_tau(double delta, double tau)
{
  double der = 0.0;

  for (int i = 1; i < 40; i++)
    der += R3Coef[i][2] * std::pow(delta, R3Coef[i][0]) * R3Coef[i][1]
            * std::pow(tau, R3Coef[i][1] - 1.0);

  return der;
}

double R3_phi_tau_tau(double delta, double tau)
{
  double der = 0.0;

  for (int i = 1; i < 40; i++)
    der += R3Coef[i][2] * std::pow(delta, R3Coef[i][0]) * R3Coef[i][1] * (R3Coef[i][1] - 1.0)
            * std::pow(tau, R3Coef[i][1] - 2.0);

  return der;
}
double R3_phi_delta_tau(double delta, double tau)
{
  double der = 0.0;

  for (int i = 1; i < 40; i++)
    der += R3Coef[i][2] * R3Coef[i][0] * std::pow(delta, R3Coef[i][0] - 1.0)
            * R3Coef[i][1] * std::pow(tau, R3Coef[i][1] - 1.0);

  return der;
}

double R3_p(double rho, double T)
{
  double delta = rho / Rhocrit;
  double tau = Tcrit / T;

  return delta * R3_phi_delta(delta, tau) * rho * Rgas * T;
}

double R3_specific_int_energy(double rho, double T)
{
  double delta = rho / Rhocrit;
  double tau = Tcrit / T;

  return tau * R3_phi_tau(delta, tau) * Rgas * T;
}

double R3_specific_entropy(double rho, double T)
{
  double delta = rho / Rhocrit;
  double tau = Tcrit / T;

  return (tau * R3_phi_tau(delta, tau) - R3_phi(delta, tau)) * Rgas;
}

double R3_specific_enthalpy(double rho, double T)
{
  double delta = rho / Rhocrit;
  double tau = Tcrit / T;

  return (tau * R3_phi_tau(delta, tau) + delta * R3_phi_delta(delta, tau)) * Rgas * T;
}

double R3_cv(double rho, double T)
{
  double delta = rho / Rhocrit;
  double tau = Tcrit / T;

  return - tau * tau * R3_phi_tau_tau(delta, tau) * Rgas;
}

double R3_cp(double rho, double T)
{
  double delta = rho / Rhocrit;
  double tau = Tcrit / T;

  double val0 = - tau * tau * R3_phi_tau_tau(delta, tau);

  double val1 = delta * R3_phi_delta(delta, tau) - delta * tau * R3_phi_delta_tau(delta, tau);
  double val2 = 2.0 * delta * R3_phi_delta(delta, tau) + delta * delta * R3_phi_delta_delta(delta, tau);

  return (val0 + val1 * val1 / val2) * Rgas;
}

double R3_sound_speed(double rho, double T)
{
  double delta = rho / Rhocrit;
  double tau = Tcrit / T;

  double val1 = delta * R3_phi_delta(delta, tau) - delta * tau * R3_phi_delta_tau(delta, tau);
  double val2 = 2.0 * delta * R3_phi_delta(delta, tau) + delta * delta * R3_phi_delta_delta(delta, tau);

  double val = val2 - val1 * val1 / (tau * tau * R3_phi_tau_tau(delta, tau));

  return std::sqrt(val * Rgas * T);
}

double p_sat_from_T(double T)
{
  double theta = T + R4Coef[8] / (T - R4Coef[9]);
  double theta2 = theta * theta;
  double A = theta2 + R4Coef[0] * theta + R4Coef[1];
  double B = R4Coef[2] * theta2 + R4Coef[3] * theta + R4Coef[4];
  double C = R4Coef[5] * theta2 + R4Coef[6] * theta + R4Coef[7];

  double val = 2.0 * C / (-B + std::sqrt(B * B - 4.0 * A * C));

  return std::pow(val, 4) * 1.0e6;
}

double T_sat_from_p(double p)
{
  double beta = std::pow(p * 1.0e-6, 0.25);
  double beta2 = beta * beta;

  double E = beta2 + R4Coef[2] * beta + R4Coef[5];
  double F = R4Coef[0] * beta2 + R4Coef[3] * beta + R4Coef[6];
  double G = R4Coef[1] * beta2 + R4Coef[4] * beta + R4Coef[7];
  double D = 2.0 * G / (-F - std::sqrt(F * F - 4 * E * G));

  double val = (R4Coef[9] + D) * (R4Coef[9] + D) - 4.0 * (R4Coef[8] + R4Coef[9] * D);
  return 0.5 * (R4Coef[9] + D - std::sqrt(val));
}

double R5_gamma_0(double pi, double tau)
{
  double gamma = std::log(pi);

  for (int i = 0; i < 6; i++)
    gamma += R5Coef0[i][1] * std::pow(tau, R5Coef0[i][0]);

  return gamma;
}

double R5_gamma_r(double pi, double tau)
{
  double gamma_r = 0.0;
  for (int i = 0; i < 6; i++)
    gamma_r += R5Coefr[i][2] * std::pow(pi, R5Coefr[i][0]) * std::pow(tau, R5Coefr[i][1]);

  return gamma_r;
}

double R5_gamma_0_pi(double pi, double /*tau*/)
{
  return 1.0 / pi;
}

double R5_gamma_0_pi_pi(double pi, double /*tau*/)
{
  return -1.0 / (pi * pi);
}

double R5_gamma_0_tau(double /*pi*/, double tau)
{
  double der = 0.0;
  for (int i = 0; i < 6; i++)
    der += R5Coef0[i][1] * R5Coef0[i][0] * std::pow(tau, R5Coef0[i][0] - 1.0);

  return der;
}

double R5_gamma_0_tau_tau(double /*pi*/, double tau)
{
  double der = 0.0;
  for (int i = 0; i < 6; i++)
    der += R5Coef0[i][1] * R5Coef0[i][0] * (R5Coef0[i][0] - 1.0)
            * std::pow(tau, R5Coef0[i][0] - 2.0);

  return der;
}

double R5_gamma_0_pi_tau(double /*pi*/, double /*tau*/)
{
  return 0.0;
}

double R5_gamma_r_pi(double pi, double tau)
{
  double der = 0.0;
  for (int i = 0; i < 6; i++)
    der += R5Coefr[i][2] * R5Coefr[i][0] * std::pow(pi, R5Coefr[i][0] - 1.0)
            * std::pow(tau, R5Coefr[i][1]);

  return der;
}

double R5_gamma_r_pi_pi(double pi, double tau)
{
  double der = 0.0;
  for (int i = 0; i < 6; i++)
    der += R5Coefr[i][2] * R5Coefr[i][0] * (R5Coefr[i][0] - 1.0) * std::pow(pi, R5Coefr[i][0] - 2.0)
            * std::pow(tau, R5Coefr[i][1]);

  return der;
}

double R5_gamma_r_tau(double pi, double tau)
{
  double der = 0.0;
  for (int i = 0; i < 6; i++)
    der += R5Coefr[i][2] * std::pow(pi, R5Coefr[i][0]) * R5Coefr[i][1]
            * std::pow(tau, R5Coefr[i][1] - 1.0);

  return der;
}

double R5_gamma_r_tau_tau(double pi, double tau)
{
  double der = 0.0;
  for (int i = 0; i < 6; i++)
    der += R5Coefr[i][2] * std::pow(pi, R5Coefr[i][0]) * R5Coefr[i][1]
            * (R5Coefr[i][1] - 1.0) * std::pow(tau, R5Coefr[i][1] - 2.0);

  return der;
}

double R5_gamma_r_pi_tau(double pi, double tau)
{
  double der = 0.0;
  for (int i = 0; i < 6; i++)
    der += R5Coefr[i][2] * R5Coefr[i][0] * std::pow(pi, R5Coefr[i][0] - 1.0)
            * R5Coefr[i][1] * std::pow(tau, R5Coefr[i][1] - 1.0);

  return der;
}

double R5_specific_volume(double p, double T)
{
  double pi = p * 1.0e-6;
  double tau = 1000.0 / T;

  return pi * (R5_gamma_0_pi(pi, tau) + R5_gamma_r_pi(pi, tau)) * Rgas * T / p;
}

double R5_specific_int_energy(double p, double T)
{
  double pi = p * 1.0e-6;
  double tau = 1000.0 / T;

  double val = tau * (R5_gamma_0_tau(pi, tau) + R5_gamma_r_tau(pi, tau))
              - pi * (R5_gamma_0_pi(pi, tau)  + R5_gamma_r_pi(pi, tau));

  return val * Rgas * T;
}

double R5_specific_entropy(double p, double T)
{
  double pi = p * 1.0e-6;
  double tau = 1000.0 / T;

  double val = tau * (R5_gamma_0_tau(pi, tau) + R5_gamma_r_tau(pi, tau))
              - (R5_gamma_0(pi, tau)  + R5_gamma_r(pi, tau));

  return val * Rgas;
}

double R5_specific_enthalpy(double p, double T)
{
  double pi = p * 1.0e-6;
  double tau = 1000.0 / T;

  return tau * (R5_gamma_0_tau(pi, tau) + R5_gamma_r_tau(pi, tau)) * Rgas * T;
}

double R5_cp(double p, double T)
{
  double pi = p * 1.0e-6;
  double tau = 1000.0 / T;

  return -tau * tau * (R5_gamma_0_tau_tau(pi, tau) + R5_gamma_r_tau_tau(pi, tau)) * Rgas;
}

double R5_cv(double p, double T)
{
  double pi = p * 1.0e-6;
  double tau = 1000.0 / T;

  double cp = -tau * tau * (R5_gamma_0_tau_tau(pi, tau) * R5_gamma_r_tau_tau(pi, tau)) * Rgas;
  double val1 = 1.0 + pi * R5_gamma_r_pi(pi, tau) - tau * pi * R5_gamma_r_pi_tau(pi, tau);
  double val2 = 1.0 - pi * pi * R5_gamma_r_pi_pi(pi, tau);

  return cp - val1 * val1 / val2;
}

double R5_sound_speed(double p, double T)
{
  double pi = p * 1.0e-6;
  double tau = 1000.0 / T;

  double gamma_r_pi = R5_gamma_r_pi(pi, tau);
  double val1 = 1.0 + 2.0 * pi * R5_gamma_r_pi(pi, tau) + pi * pi * gamma_r_pi * gamma_r_pi;
  double val2 = 1.0 - pi * pi * R5_gamma_r_pi_pi(pi, tau);
  double val3 = 1.0 + pi * R5_gamma_r_pi(pi, tau) - tau * pi * R5_gamma_r_pi_tau(pi, tau);
  double val4 = tau * tau * (R5_gamma_0_tau_tau(pi, tau) + R5_gamma_r_tau_tau(pi, tau));

  double w2 = val1 / (val2 + val3 * val3 / val4) * Rgas * T;

  return std::sqrt(w2);
}
