#include "IAPWS1995_Rev.h"
#include <math.h>
#include <stdio.h>
#include <iostream>

//phi_0 and its derivatives
double
IAPWS1995Rev::func_phi_0(const double Temperature, const double rho)
{
  //non dimensional density and temperature
  double delta = rho / Rho_critical;
  double tau = T_critical / Temperature;

  //calculate value of phi_0  (Equation (5), page 3, Ref. [1])
  double phi_0 = log(delta) + n0[0] + n0[1] * tau + n0[2] * log(tau);

  for(int i = 4; i <= 8; i++)
    phi_0 += n0[i-1] * log(1.0 - exp(-gamma0[i-4] * tau));

  return phi_0;
}

double
IAPWS1995Rev::func_d_phi0_d_delta(const double Temperature, const double rho)
{
  // Talbe 4, page 12, Ref. [1]
  double delta = rho / Rho_critical;

  return 1.0 / delta;
}

double
IAPWS1995Rev::func_d_phi0_d_delta2(const double Temperature, const double rho)
{
  // Talbe 4, page 12, Ref. [1]
  double delta = rho / Rho_critical;

  return -1.0 / (delta * delta);
}

double
IAPWS1995Rev::func_d_phi0_d_tau(const double Temperature, const double rho)
{
  // Talbe 4, page 12, Ref. [1]
  double tau = T_critical / Temperature;

  double d_phi0_d_tau = n0[1] + n0[2] / tau;
  for(int i = 4; i <= 8; i++)
    d_phi0_d_tau += n0[i-1] * gamma0[i-4] * (1.0 / (1.0 - exp(-gamma0[i-4]*tau)) - 1.0);

  return d_phi0_d_tau;
}

double
IAPWS1995Rev::func_d_phi0_d_tau2(const double Temperature, const double rho)
{
  // Talbe 4, page 12, Ref. [1]
  double tau = T_critical / Temperature;

  double d_phi0_d_tau2 = -n0[2] / (tau * tau);
  for(int i = 4; i <= 8; i++)
  {
    double temp_val = exp(-gamma0[i-4] * tau);
    d_phi0_d_tau2 -= n0[i-1] * gamma0[i-4] * gamma0[i-4] * temp_val / ((1.0 - temp_val) * (1.0 - temp_val));
  }

  return d_phi0_d_tau2;
}

double
IAPWS1995Rev::func_d_phi0_d_tau_d_delta(const double Temperature, const double rho)
{
  // Talbe 4, page 12, Ref. [1]
  return 0.0;
}

void
IAPWS1995Rev::phi_0_and_its_derivatives(const double Temperature, const double rho,
              double & phi0,
              double & d_phi0_d_delta,
              double & d_phi0_d_delta2,
              double & d_phi0_d_tau,
              double & d_phi0_d_tau2,
              double & d_phi0_d_tau_d_delta)
{
  phi0 = func_phi_0(Temperature, rho);
  d_phi0_d_delta = func_d_phi0_d_delta(Temperature, rho);
  d_phi0_d_delta2 = func_d_phi0_d_delta2(Temperature, rho);
  d_phi0_d_tau = func_d_phi0_d_tau(Temperature, rho);
  d_phi0_d_tau2 = func_d_phi0_d_tau2(Temperature, rho);
  d_phi0_d_tau_d_delta = func_d_phi0_d_tau_d_delta(Temperature, rho);
}


//phi_r and its derivatives
double
IAPWS1995Rev::func_phi_r(const double Temperature, const double rho)
{
  //non dimensional density and temperature
  double delta = rho / Rho_critical;
  double tau = T_critical / Temperature;

  //calculate value of phi_r  (Equation (6), page 4, Ref. [1])
  double phi_r = 0.0;
  for(int i = 1; i <= 7; i++)
    phi_r += n1[i-1] * pow(delta, d1[i-1]) * pow(tau, t1[i-1]);

  for(int i = 8; i <= 51; i++)
    phi_r += n2[i-8] * pow(delta, d2[i-8]) * pow(tau, t2[i-8]) * exp(-pow(delta, c2[i-8]));

  for(int i = 52; i <= 54; i++)
  {
    double coeff = - alpha3[i-52] * (delta - epsil3[i-52]) * (delta - epsil3[i-52])
                   - beta3[i-52] * (tau - gamma3[i-52]) * (tau - gamma3[i-52]);

    phi_r += n3[i-52] * pow(delta, d3[i-52]) * pow(tau, t3[i-52]) * exp(coeff);
  }

  for(int i = 55; i <= 56; i++)
  {
    double theta = (1.0 - tau) + AA4[i-55] * pow((delta - 1.0) * (delta - 1.0), 0.5 / beta4[i-55]);
    double DELTA = theta * theta + BB4[i-55] * pow((delta - 1.0) * (delta - 1.0), a4[i-55]);
    double PSI = exp(-CC4[i-55] * (delta - 1.0) * (delta - 1.0) - DD4[i-55] * (tau - 1.0) * (tau - 1.0));
    phi_r += n4[i-55] * pow(DELTA, b4[i-55]) * delta * PSI;
  }

  return phi_r;
}

double
IAPWS1995Rev::func_d_phir_d_delta(const double Temperature, const double rho)
{
  double delta = rho / Rho_critical;
  double tau = T_critical / Temperature;

  double d_phir_d_delta = 0.0;

  // Table 5, page 13, Ref. [1]
  for(int i = 1; i <= 7; i++)
    d_phir_d_delta += n1[i-1] * d1[i-1] * pow(delta, d1[i-1] - 1.0) * pow(tau, t1[i-1]);

  for(int i = 8; i <= 51; i++)
    d_phir_d_delta += n2[i-8] * exp(-pow(delta, c2[i-8])) * pow(delta, d2[i-8] - 1.0) * pow(tau, t2[i-8]) * (d2[i-8] - c2[i-8] * pow(delta, c2[i-8]));

  for(int i = 52; i <= 54; i++)
  {
    double coeff = - alpha3[i-52] * (delta - epsil3[i-52]) * (delta - epsil3[i-52])
                   - beta3[i-52] * (tau - gamma3[i-52]) * (tau - gamma3[i-52]);

    d_phir_d_delta += n3[i-52] * pow(delta, d3[i-52]) * pow(tau, t3[i-52]) * exp(coeff) * (d3[i-52]/delta - 2.0 * alpha3[i-52] * (delta-epsil3[i-52]));
  }

  for(int i = 55; i <= 56; i++)
  {
    double theta = (1.0 - tau) + AA4[i-55] * pow((delta - 1.0) * (delta - 1.0), 0.5 / beta4[i-55]);
    double DELTA = theta * theta + BB4[i-55] * pow((delta - 1.0) * (delta - 1.0), a4[i-55]);
    double PSI = exp(-CC4[i-55] * (delta - 1.0) * (delta - 1.0) - DD4[i-55] * (tau - 1.0) * (tau - 1.0));
    double d_DELTA_d_delta = (delta - 1.0) *
          (
            AA4[i-55] * theta * 2.0 / beta4[i-55] * pow((delta-1.0) * (delta-1.0), 0.5/beta4[i-55] - 1.0)
            + 2.0 * BB4[i-55] *a4[i-55] * pow((delta-1.0)*(delta-1.0), a4[i-55] - 1.0)
          );

    double d_PSI_d_delta = -2.0 * CC4[i-55] * (delta - 1.0) * PSI;
    double d_DELTA_pow_bi_d_delta = b4[i-55] * pow(DELTA, b4[i-55] - 1.0) * d_DELTA_d_delta;

    d_phir_d_delta += n4[i-55] * (pow(DELTA, b4[i-55]) * (PSI + delta * d_PSI_d_delta) + d_DELTA_pow_bi_d_delta * delta * PSI);
  }

  return d_phir_d_delta;
}

double
IAPWS1995Rev::func_d_phir_d_delta2(const double Temperature, const double rho)
{
  double delta = rho / Rho_critical;
  double tau = T_critical / Temperature;

  double d_phir_d_delta2 = 0.0;

  // Table 5, page 13, Ref. [1]
  for(int i = 1; i <= 7; i++)
    d_phir_d_delta2 += n1[i-1] * d1[i-1] * (d1[i-1] - 1.0) * pow(delta, d1[i-1] - 2.0) * pow(tau, t1[i-1]);

  for(int i = 8; i <= 51; i++)
    d_phir_d_delta2 += n2[i-8] * exp(-pow(delta, c2[i-8])) * pow(delta, d2[i-8] - 2.0) * pow(tau, t2[i-8])
            * (
                (d2[i-8] - c2[i-8] * pow(delta, c2[i-8])) * (d2[i-8] - 1.0 - c2[i-8] * pow(delta, c2[i-8]))
                - c2[i-8] * c2[i-8] * pow(delta, c2[i-8])
              );

  for(int i = 52; i <= 54; i++)
  {
    double coeff = - alpha3[i-52] * (delta - epsil3[i-52]) * (delta - epsil3[i-52])
                   - beta3[i-52] * (tau - gamma3[i-52]) * (tau - gamma3[i-52]);

    d_phir_d_delta2 += n3[i-52] * pow(tau, t3[i-52]) * exp(coeff)
          * (
            - 2.0 * alpha3[i-52] * pow(delta, d3[i-52])
            + 4.0 * alpha3[i-52] * alpha3[i-52] * pow(delta, d3[i-52]) * (delta - epsil3[i-52]) * (delta - epsil3[i-52])
            - 4.0 * d3[i-52] * alpha3[i-52] * pow(delta, d3[i-52] - 1.0) * (delta - epsil3[i-52])
            + d3[i-52] * (d3[i-52] - 1.0) * pow(delta, d3[i-52] - 2.0)
          );
  }

  for(int i = 55; i <= 56; i++)
  {
    double theta = (1.0 - tau) + AA4[i-55] * pow((delta - 1.0) * (delta - 1.0), 0.5 / beta4[i-55]);
    double DELTA = theta * theta + BB4[i-55] * pow((delta - 1.0) * (delta - 1.0), a4[i-55]);
    double PSI = exp(-CC4[i-55] * (delta - 1.0) * (delta - 1.0) - DD4[i-55] * (tau - 1.0) * (tau - 1.0));

    double d_DELTA_d_delta = (delta - 1.0) *
          (
            AA4[i-55] * theta * 2.0 / beta4[i-55] * pow((delta-1.0)*(delta-1.0), 0.5/beta4[i-55]-1.0)
            + 2.0 * BB4[i-55] *a4[i-55] * pow((delta-1.0)*(delta-1.0), a4[i-55] - 1.0)
          );
    double d_DELTA2_d_delta2 = d_DELTA_d_delta / (delta - 1.0)
            + (delta - 1.0) * (delta - 1.0) *
              (
                4.0 * BB4[i-55] * a4[i-55] * (a4[i-55] - 1.0) * pow((delta-1.0)*(delta-1.0), a4[i-55] - 2.0)
                    + 2.0 * AA4[i-55] * AA4[i-55] / beta4[i-55] / beta4[i-55] * pow((delta-1.0) * (delta-1.0), 1.0/beta4[i-55] - 2.0)
                    + AA4[i-55] * theta * 4.0 / beta4[i-55] * (0.5/beta4[i-55] - 1.0) * pow((delta-1.0)*(delta-1.0), 0.5/beta4[i-55] - 2.0)
              );


    double d_PSI_d_delta = -2.0 * CC4[i-55] * (delta - 1.0) * PSI;
    double d_PSI2_d_delta2 = (2.0 * CC4[i-55] * (delta - 1.0) * (delta - 1.0) - 1.0) * 2.0 * CC4[i-55] * PSI;

    double d_DELTA_pow_bi_d_delta = b4[i-55] * pow(DELTA, b4[i-55] - 1.0) * d_DELTA_d_delta;
    double d_DELTA_pow_bi_d_delta2 = b4[i-55] * (
                    pow(DELTA, b4[i-55] - 1.0) * d_DELTA2_d_delta2
                    + (b4[i-55] - 1.0) * pow(DELTA, b4[i-55] - 2.0) * d_DELTA_d_delta * d_DELTA_d_delta
                   );

    d_phir_d_delta2 += n4[i-55] * (
              pow(DELTA, b4[i-55]) * (2.0 * d_PSI_d_delta + delta * d_PSI2_d_delta2)
                  + 2.0 * d_DELTA_pow_bi_d_delta * (PSI + delta * d_PSI_d_delta)
                  + d_DELTA_pow_bi_d_delta2 * delta * PSI
                );
  }

  return d_phir_d_delta2;
}

double
IAPWS1995Rev::func_d_phir_d_tau(const double Temperature, const double rho)
{
  double delta = rho / Rho_critical;
  double tau = T_critical / Temperature;

  double d_phir_d_tau = 0.0;

  // Table 5, page 13, Ref. [1]
  for(int i = 1; i <= 7; i++)
    d_phir_d_tau += n1[i-1] * t1[i-1] * pow(delta, d1[i-1]) * pow(tau, t1[i-1] - 1.0);

  for(int i = 8; i <= 51; i++)
    d_phir_d_tau += n2[i-8] * t2[i-8] * pow(delta, d2[i-8]) * pow(tau, t2[i-8] - 1.0) * exp(-pow(delta, c2[i-8]));

  for(int i = 52; i <= 54; i++)
  {
    double coeff = - alpha3[i-52] * (delta - epsil3[i-52]) * (delta - epsil3[i-52])
                   - beta3[i-52] * (tau - gamma3[i-52]) * (tau - gamma3[i-52]);

    d_phir_d_tau += n3[i-52] * pow(delta, d3[i-52]) * pow(tau, t3[i-52]) * exp(coeff) * (t3[i-52] / tau - 2.0 * beta3[i-52] * (tau - gamma3[i-52]));
  }

  for(int i = 55; i <= 56; i++) {
    double theta = (1.0 - tau) + AA4[i-55] * pow((delta - 1.0) * (delta - 1.0), 0.5 / beta4[i-55]);
    double DELTA = theta * theta + BB4[i-55] * pow((delta - 1.0) * (delta - 1.0), a4[i-55]);
    double PSI = exp(-CC4[i-55] * (delta - 1.0) * (delta - 1.0) - DD4[i-55] * (tau - 1.0) * (tau - 1.0));

    double d_PSI_d_tau = -2.0 * DD4[i-55] * (tau - 1.0) * PSI;
    double d_DELTA_pow_bi_d_tau = -2.0 * theta * b4[i-55] * pow(DELTA, b4[i-55] - 1.0);

    d_phir_d_tau += n4[i-55] * delta * (d_DELTA_pow_bi_d_tau * PSI + pow(DELTA, b4[i-55]) * d_PSI_d_tau);
  }

  return d_phir_d_tau;
}

double
IAPWS1995Rev::func_d_phir_d_tau2(const double Temperature, const double rho)
{
  double delta = rho / Rho_critical;
  double tau = T_critical / Temperature;

  double d_phir_d_tau2 = 0.0;

  // Table 5, page 13, Ref. [1]
  for(int i = 1; i <= 7; i++)
    d_phir_d_tau2 += n1[i-1] * t1[i-1] * (t1[i-1] - 1.0) * pow(delta, d1[i-1]) * pow(tau, t1[i-1] - 2.0);

  for(int i = 8; i <= 51; i++)
    d_phir_d_tau2 += n2[i-8] * t2[i-8] * (t2[i-8] - 1.0) * pow(delta, d2[i-8]) * pow(tau, t2[i-8] - 2.0) * exp(-pow(delta, c2[i-8]));

  for(int i = 52; i <= 54; i++)
  {
    double coeff = - alpha3[i-52] * (delta - epsil3[i-52]) * (delta - epsil3[i-52])
                   - beta3[i-52] * (tau - gamma3[i-52]) * (tau - gamma3[i-52]);

    d_phir_d_tau2 += n3[i-52] * pow(delta, d3[i-52]) * pow(tau, t3[i-52]) * exp(coeff) *
            (
              pow((t3[i-52] / tau - 2.0 * beta3[i-52] * (tau - gamma3[i-52])), 2.0)
              - t3[i-52] / tau / tau - 2.0 * beta3[i-52]
            );
  }

  for(int i = 55; i <= 56; i++)
  {
    double theta = (1.0 - tau) + AA4[i-55] * pow((delta - 1.0) * (delta - 1.0), 0.5 / beta4[i-55]);
    double DELTA = theta * theta + BB4[i-55] * pow((delta - 1.0) * (delta - 1.0), a4[i-55]);
    double PSI = exp(-CC4[i-55] * (delta - 1.0) * (delta - 1.0) - DD4[i-55] * (tau - 1.0) * (tau - 1.0));

    double d_PSI_d_tau = -2.0 * DD4[i-55] * (tau - 1.0) * PSI;
    double d_PSI_d_tau2 = (2.0 * DD4[i-55] * (tau - 1.0) * (tau - 1.0) - 1.0) * 2.0 * DD4[i-55] * PSI;
    double d_DELTA_pow_bi_d_tau = -2.0 * theta * b4[i-55] * pow(DELTA, b4[i-55] - 1.0);
    double d_DELTA_pow_bi_d_tau2 = 2.0 * b4[i-55] * pow(DELTA, b4[i-55] - 1.0) + 4.0 * theta * theta * b4[i-55] * (b4[i-55] - 1.0) * pow(DELTA, b4[i-55] - 2.0);

    d_phir_d_tau2 += n4[i-55] * delta * (d_DELTA_pow_bi_d_tau2 * PSI + 2.0 * d_DELTA_pow_bi_d_tau * d_PSI_d_tau + pow(DELTA, b4[i-55]) * d_PSI_d_tau2);
  }

  return d_phir_d_tau2;
}

double
IAPWS1995Rev::func_d_phir_d_tau_d_delta(const double Temperature, const double rho)
{
  double delta = rho / Rho_critical;
  double tau = T_critical / Temperature;

  double d_phir_d_tau_d_delta = 0.0;

  // Table 5, page 13, Ref. [1]
  for(int i = 1; i <= 7; i++)
    d_phir_d_tau_d_delta += n1[i-1] * d1[i-1] * t1[i-1] * pow(delta, d1[i-1] - 1.0) * pow(tau, t1[i-1] - 1.0);

  for(int i = 8; i <= 51; i++)
    d_phir_d_tau_d_delta += n2[i-8] * t2[i-8] * pow(delta, d2[i-8] - 1.0) * pow(tau, t2[i-8] - 1.0)
          * (d2[i-8] - c2[i-8] * pow(delta, c2[i-8])) * exp(-pow(delta, c2[i-8]));

  for(int i = 52; i <= 54; i++)
  {
    double coeff = - alpha3[i-52] * (delta - epsil3[i-52]) * (delta - epsil3[i-52])
                   - beta3[i-52] * (tau - gamma3[i-52]) * (tau - gamma3[i-52]);

    d_phir_d_tau_d_delta += n3[i-52] * pow(delta, d3[i-52]) * pow(tau, t3[i-52]) * exp(coeff)
          * ( d3[i-52] / delta - 2.0 * alpha3[i-52] * (delta - epsil3[i-52]) )
          * ( t3[i-52] / tau   - 2.0 *  beta3[i-52] * (tau   - gamma3[i-52]) );
  }

  for(int i = 55; i <= 56; i++)
  {
    double theta = (1.0 - tau) + AA4[i-55] * pow((delta - 1.0) * (delta - 1.0), 0.5 / beta4[i-55]);
    double DELTA = theta * theta + BB4[i-55] * pow((delta - 1.0) * (delta - 1.0), a4[i-55]);
    double PSI = exp(-CC4[i-55] * (delta - 1.0) * (delta - 1.0) - DD4[i-55] * (tau - 1.0) * (tau - 1.0));

    double d_DELTA_d_delta = (delta - 1.0) *
          (
            AA4[i-55] * theta * 2.0 / beta4[i-55] * pow((delta-1.0)*(delta-1.0), 0.5/beta4[i-55]-1.0)
            + 2.0 * BB4[i-55] *a4[i-55] * pow((delta-1.0)*(delta-1.0), a4[i-55] - 1.0)
          );
    double d_PSI_d_delta = -2.0 * CC4[i-55] * (delta - 1.0) * PSI;
    double d_PSI2_d_delta2 = (2.0 * CC4[i-55] * (delta - 1.0) * (delta - 1.0) - 1.0) * 2.0 * CC4[i-55] * PSI;
    double d_PSI_d_tau = -2.0 * DD4[i-55] * (tau - 1.0) * PSI;
    double d_PSI_d_tau2 = (2.0 * DD4[i-55] * (tau - 1.0) * (tau - 1.0) - 1.0) * 2.0 * DD4[i-55] * PSI;
    double d_PSI_d_delta_d_tau = 4.0 * CC4[i-55] * DD4[i-55] * (delta - 1.0) * (tau - 1.0) * PSI;

    double d_DELTA_pow_bi_d_delta = b4[i-55] * pow(DELTA, b4[i-55] - 1.0) * d_DELTA_d_delta;
    double d_DELTA_pow_bi_d_tau = -2.0 * theta * b4[i-55] * pow(DELTA, b4[i-55] - 1.0);
    double d_DELTA_pow_bi_d_delta_d_tau = - AA4[i-55] * b4[i-55] * 2.0 / beta4[i-55] * pow(DELTA, b4[i-55] - 1.0) * (delta - 1.0) *
                pow((delta - 1.0)*(delta - 1.0), 0.5/beta4[i-55] - 1.0)
                  - 2.0 * theta * b4[i-55] * (b4[i-55] - 1.0) * pow(DELTA, b4[i-55] - 2.0) * d_DELTA_d_delta;

    d_phir_d_tau_d_delta += n4[i-55] * (
              pow(DELTA, b4[i-55]) * (d_PSI_d_tau + delta * d_PSI_d_delta_d_tau)
                  + delta * d_DELTA_pow_bi_d_delta * d_PSI_d_tau
                  + d_DELTA_pow_bi_d_tau * (PSI + delta * d_PSI_d_delta)
                  + d_DELTA_pow_bi_d_delta_d_tau * delta * PSI
               );
  }

  return d_phir_d_tau_d_delta;
}

void
IAPWS1995Rev::phi_r_and_its_derivatives(const double Temperature, const double rho,
              double & phi_r,
              double & d_phir_d_delta,
              double & d_phir_d_delta2,
              double & d_phir_d_tau,
              double & d_phir_d_tau2,
              double & d_phir_d_tau_d_delta)
{
  phi_r = func_phi_r(Temperature, rho);
  d_phir_d_delta = func_d_phir_d_delta(Temperature, rho);
  d_phir_d_delta2 = func_d_phir_d_delta2(Temperature, rho);
  d_phir_d_tau = func_d_phir_d_tau(Temperature, rho);
  d_phir_d_tau2 = func_d_phir_d_tau2(Temperature, rho);
  d_phir_d_tau_d_delta = func_d_phir_d_tau_d_delta(Temperature, rho);
}

//Thermodynamics properties
double
IAPWS1995Rev::Pressure(const double Temperature, const double rho)
{
  double delta = rho / Rho_critical;
  double tau = T_critical / Temperature;
  double d_phir_d_delta = IAPWS1995Rev::func_d_phir_d_delta(Temperature, rho);

  return (1.0 + delta * d_phir_d_delta) * rho * IAPWS1995Rev::R * Temperature;
}


double
IAPWS1995Rev::Entropy(const double Temperature, const double rho)
{
  //double delta = rho / Rho_critical;
  double tau = T_critical / Temperature;

  double phi0 = IAPWS1995Rev::func_phi_0(Temperature, rho);
  double d_phi0_d_tau = IAPWS1995Rev::func_d_phi0_d_tau(Temperature, rho);

  double phi_r = IAPWS1995Rev::func_phi_r(Temperature, rho);
  double d_phir_d_tau = IAPWS1995Rev::func_d_phir_d_tau(Temperature, rho);

  return IAPWS1995Rev::R * (tau * (d_phi0_d_tau + d_phir_d_tau) - phi0 - phi_r);
}

double
IAPWS1995Rev::IntEnergy(const double Temperature, const double rho)
{
  //double delta = rho / Rho_critical;
  double tau = T_critical / Temperature;

  double d_phi0_d_tau = IAPWS1995Rev::func_d_phi0_d_tau(Temperature, rho);
  double d_phir_d_tau = IAPWS1995Rev::func_d_phir_d_tau(Temperature, rho);

  return IAPWS1995Rev::R * Temperature * tau * (d_phi0_d_tau + d_phir_d_tau);
}

double
IAPWS1995Rev::Enthalpy(const double Temperature, const double rho)
{
  double delta = rho / Rho_critical;
  double tau = T_critical / Temperature;

  double d_phi0_d_tau = IAPWS1995Rev::func_d_phi0_d_tau(Temperature, rho);
  double d_phir_d_tau = IAPWS1995Rev::func_d_phir_d_tau(Temperature, rho);
  double d_phir_d_delta = IAPWS1995Rev::func_d_phir_d_delta(Temperature, rho);

  double value = 1.0 + tau * (d_phi0_d_tau + d_phir_d_tau) + delta * d_phir_d_delta;

  return IAPWS1995Rev::R * Temperature * value;
}

double
IAPWS1995Rev::GibbsFreeEnergy(const double Temperature, const double rho)
{
  double delta = rho / Rho_critical;

  double phi_0 = IAPWS1995Rev::func_phi_0(Temperature, rho);
  double phi_r = IAPWS1995Rev::func_phi_r(Temperature, rho);
  double d_phir_d_delta = IAPWS1995Rev::func_d_phir_d_delta(Temperature, rho);

  return IAPWS1995Rev::R * Temperature * (1.0 + phi_0 + phi_r + delta * d_phir_d_delta);
}

double
IAPWS1995Rev::Cv(const double Temperature, const double rho)
{
  double tau = T_critical / Temperature;

  double d_phi0_d_tau2 = IAPWS1995Rev::func_d_phi0_d_tau2(Temperature, rho);
  double d_phir_d_tau2 = IAPWS1995Rev::func_d_phir_d_tau2(Temperature, rho);

  return -IAPWS1995Rev::R * tau * tau * (d_phi0_d_tau2 + d_phir_d_tau2);
}

double
IAPWS1995Rev::Cp(const double Temperature, const double rho)
{
  double delta = rho / Rho_critical;
  double tau = T_critical / Temperature;

  double d_phi0_d_tau2 = IAPWS1995Rev::func_d_phi0_d_tau2(Temperature, rho);
  double d_phir_d_tau2 = IAPWS1995Rev::func_d_phir_d_tau2(Temperature, rho);
  double d_phir_d_delta = IAPWS1995Rev::func_d_phir_d_delta(Temperature, rho);
  double d_phir_d_delta2 = IAPWS1995Rev::func_d_phir_d_delta2(Temperature, rho);
  double d_phir_d_tau_d_delta = IAPWS1995Rev::func_d_phir_d_tau_d_delta(Temperature, rho);

  return -IAPWS1995Rev::R * tau * tau * (d_phi0_d_tau2 + d_phir_d_tau2)
        + pow((1.0 + delta * d_phir_d_delta - delta * tau * d_phir_d_tau_d_delta), 2.0)
          / (1.0 + 2.0 * delta * d_phir_d_delta + delta * delta * d_phir_d_delta2);
}

double
IAPWS1995Rev::SoundSpeed(const double Temperature, const double rho)
{
  double delta = rho / Rho_critical;
  double tau = T_critical / Temperature;

  double d_phi0_d_tau2 = IAPWS1995Rev::func_d_phi0_d_tau2(Temperature, rho);
  double d_phir_d_tau2 = IAPWS1995Rev::func_d_phir_d_tau2(Temperature, rho);
  double d_phir_d_delta = IAPWS1995Rev::func_d_phir_d_delta(Temperature, rho);
  double d_phir_d_delta2 = IAPWS1995Rev::func_d_phir_d_delta2(Temperature, rho);
  double d_phir_d_tau_d_delta = IAPWS1995Rev::func_d_phir_d_tau_d_delta(Temperature, rho);

  double c_squared = 1.0 + 2.0 * delta * d_phir_d_delta + delta * delta * d_phir_d_delta2
          - pow((1.0 + delta * d_phir_d_delta - delta * tau * d_phir_d_tau_d_delta), 2.0)
            / tau / tau / (d_phi0_d_tau2 + d_phir_d_tau2);

  c_squared *= IAPWS1995Rev::R * Temperature;

  if(c_squared >= 0.0)
    return sqrt(c_squared);
  else
  {
    fprintf(stderr, "%s", "c_squared < 0!\n");
    exit(1);
  }
  return 0;
}

double
IAPWS1995Rev::JouleThomson_Coeff(const double Temperature, const double rho)
{
  //non dimensional density and temperature
  double delta = rho / Rho_critical;
  double tau = T_critical / Temperature;

  double d_phi0_d_tau2 = IAPWS1995Rev::func_d_phi0_d_tau2(Temperature, rho);
  double d_phir_d_tau2 = IAPWS1995Rev::func_d_phir_d_tau2(Temperature, rho);
  double d_phir_d_delta = IAPWS1995Rev::func_d_phir_d_delta(Temperature, rho);
  double d_phir_d_delta2 = IAPWS1995Rev::func_d_phir_d_delta2(Temperature, rho);
  double d_phir_d_tau_d_delta = IAPWS1995Rev::func_d_phir_d_tau_d_delta(Temperature, rho);

  double value = -(delta * d_phir_d_delta + delta * delta * d_phir_d_delta2 + delta * tau * d_phir_d_tau_d_delta) /
          (
            pow((1.0 + delta * d_phir_d_delta - delta * tau * d_phir_d_tau_d_delta), 2.0)
            - tau * tau * (d_phi0_d_tau2 + d_phir_d_tau2) * (1.0 + 2.0 * delta * d_phir_d_delta + delta * delta * d_phir_d_delta2)
          );

  return value / (IAPWS1995Rev::R * rho);
}

double
IAPWS1995Rev::IsoThrottling_Coeff(const double Temperature, const double rho)
{
  double delta = rho / Rho_critical;
  double tau = T_critical / Temperature;
  double d_phir_d_delta = IAPWS1995Rev::func_d_phir_d_delta(Temperature, rho);
  double d_phir_d_delta2 = IAPWS1995Rev::func_d_phir_d_delta2(Temperature, rho);
  double d_phir_d_tau_d_delta = IAPWS1995Rev::func_d_phir_d_tau_d_delta(Temperature, rho);

  double delta_T = 1.0 - (1.0 + delta * d_phir_d_delta - delta * tau * d_phir_d_tau_d_delta) / (1.0 + 2.0 * delta * d_phir_d_delta + delta * delta * d_phir_d_delta2);
  return delta_T / rho;
}

double
IAPWS1995Rev::Isentropic_TemperaturePressure_Coeff(const double Temperature, const double rho)
{
  double delta = rho / Rho_critical;
  double tau = T_critical / Temperature;

  double d_phi0_d_tau2 = IAPWS1995Rev::func_d_phi0_d_tau2(Temperature, rho);
  double d_phir_d_tau2 = IAPWS1995Rev::func_d_phir_d_tau2(Temperature, rho);
  double d_phir_d_delta = IAPWS1995Rev::func_d_phir_d_delta(Temperature, rho);
  double d_phir_d_delta2 = IAPWS1995Rev::func_d_phir_d_delta2(Temperature, rho);
  double d_phir_d_tau_d_delta = IAPWS1995Rev::func_d_phir_d_tau_d_delta(Temperature, rho);

  double beta_s = (1.0 + delta * d_phir_d_delta - delta * tau * d_phir_d_tau_d_delta) /
          (
            pow((1.0 + delta * d_phir_d_delta - delta * tau * d_phir_d_tau_d_delta), 2.0)
            - tau * tau * (d_phi0_d_tau2 + d_phir_d_tau2) * (1.0 + 2.0 * delta * d_phir_d_delta + delta * delta * d_phir_d_delta2)
          );

  return beta_s / (IAPWS1995Rev::R * rho);
}

double
IAPWS1995Rev::Second_Virial_Coeff(const double Temperature, const double rho)
{
  // TODO: Table 3, page 11. Ref. [1]
  return 0;
}

double
IAPWS1995Rev::Third_Virial_Coeff(const double Temperature, const double rho)
{
  // TODO: Table 3, page 11. Ref. [1]
  return 0;
}
