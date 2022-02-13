#include "IAPWS1995_Interface.h"
#include "IAPWS1995_Rev.h"
#include <iostream>
#include <iomanip>
#include <cmath>

double
IAPWS1995Rev::rho_l_from_pT(double p, double T)
{
  double Tsat = Tsat_from_p(p);
  double rho_l_sat = rho_l_sat_from_T(Tsat);

  double rho_min, rho_max, rho_guess;
  int it = 0;           int it_max = 1000;
  double rtol = 1e-9;   double error = 1.0;

  // setup the [min, max] bound
  if (T < Tsat)  // subcooled
  {
    rho_min = rho_l_sat;
    rho_max = 1100.0; // random max value
  }
  else  // superheated
  {
    rho_min = rho_g_sat_from_T(Tsat);
    rho_max = rho_l_sat;
  }

  // iterate to find correct value
  while (it < it_max)
  {
    rho_guess = 0.5 * (rho_min + rho_max);
    error = (Pressure(T, rho_guess) - p) / p;

    if (abs(error) < rtol)
      break;
    else
    {
      if (error > 0.0)    rho_max = rho_guess;
      else                rho_min = rho_guess;
    }
    it++;
  }

  if (it == it_max)
    std::cerr << "Max it number reached in IAPWS1995Rev::rho_l_from_pT\n";
  return rho_guess;
}

double
IAPWS1995Rev::e_l_from_pT(double p, double T)
{
  return 0;
}

double
IAPWS1995Rev::h_l_from_pT(double p, double T)
{
  return 0;
}

double
IAPWS1995Rev::s_l_from_pT(double p, double T)
{
  return 0;
}

double
IAPWS1995Rev::cv_l_from_pT(double p, double T)
{
  return 0;
}

double
IAPWS1995Rev::cp_l_from_pT(double p, double T)
{
  return 0;
}

double
IAPWS1995Rev::c_l_from_pT(double p, double T)
{
  return 0;
}

double
IAPWS1995Rev::k_l_from_pT(double p, double T)
{
  return 0;
}

double
IAPWS1995Rev::mu_l_from_pT(double p, double T)
{
  return 0;
}


double
IAPWS1995Rev::rho_g_from_pT(double p, double T)
{
  return 0;
}

double
IAPWS1995Rev::e_g_from_pT(double p, double T)
{
  return 0;
}

double
IAPWS1995Rev::h_g_from_pT(double p, double T)
{
  return 0;
}

double
IAPWS1995Rev::s_g_from_pT(double p, double T)
{
  return 0;
}

double
IAPWS1995Rev::cv_g_from_pT(double p, double T)
{
  return 0;
}

double
IAPWS1995Rev::cp_g_from_pT(double p, double T)
{
  return 0;
}

double
IAPWS1995Rev::c_g_from_pT(double p, double T)
{
  return 0;
}

double
IAPWS1995Rev::k_g_from_pT(double p, double T)
{
  return 0;
}

double
IAPWS1995Rev::mu_g_from_pT(double p, double T)
{
  return 0;
}
