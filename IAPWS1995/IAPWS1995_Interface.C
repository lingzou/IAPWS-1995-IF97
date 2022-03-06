#include <iostream>
#include <iomanip>
#include <cmath>

#include "IAPWS1995_Interface.h"
#include "IAPWS1995_Rev.h"
#include "ThermalConductivity.h"
#include "Viscosity.h"

#include "IF97.h" // for metastable vapor
#include "IF97_interface.h"

void
IAPWS1995Rev::rho_l_from_pT(double p, double T, double &rho, bool &search_failed)
{
  double Tsat = IAPWS1995Rev::Tsat_from_p(p);
  double rho_l_sat = IAPWS1995Rev::rho_l_sat_from_T(Tsat);

  double rho_min, rho_max, rho_guess;
  int it = 0;           int it_max = 1000;
  double rtol = 1e-9;   double error = 1.0;

  // setup the [min, max] bound
  if (T < Tsat)  // subcooled (stable)
  {
    rho_min = rho_from_pT(p, T) - 2.0; //IF97 call
    rho_max = rho_from_pT(p, T) + 2.0; //IF97 call
  }
  else  // superheated (metastable)
  {
    // from rho_max = rho_l_sat,
    // keep searching a [rho-1, rho] bracket
    rho_max = rho_l_sat + 0.5; // tolerance
    rho_min = rho_max - 1.0;
    for (int i = 0; i < int(rho_max); i++)
    {
      if (Pressure(T, rho_min) < p)
        break;
      else
      {
        rho_max = rho_min;
        rho_min = rho_max - 1.0;
      }

      if (i == int(rho_max) - 1) // we know [rho_min, rho_max] bracket does not capture the p
        it = 999; // so let's not waste time in the coming iteration to find rho
    }
  }

  // bi-section iterate to find correct value
  while (it < it_max)
  {
    rho_guess = 0.5 * (rho_min + rho_max);
    double p_guess = Pressure(T, rho_guess);
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
  {
    std::cerr << "Max it number reached in IAPWS1995Rev::rho_l_from_pT\n";
    rho = 0;
    search_failed = true;
  }
  else
  {
    rho = rho_guess;
    search_failed = false;
  }
}

double
IAPWS1995Rev::e_l_from_pT(double p, double T)
{
  double rho_l = 0; bool search_failed = true;
  rho_l_from_pT(p, T, rho_l, search_failed);
  return search_failed ? 0 : IntEnergy(T, rho_l);
}

double
IAPWS1995Rev::h_l_from_pT(double p, double T)
{
  double rho_l = 0; bool search_failed = true;
  rho_l_from_pT(p, T, rho_l, search_failed);
  return search_failed ? 0 : Enthalpy(T, rho_l);
}

double
IAPWS1995Rev::s_l_from_pT(double p, double T)
{
  double rho_l = 0; bool search_failed = true;
  rho_l_from_pT(p, T, rho_l, search_failed);
  return search_failed ? 0 : Entropy(T, rho_l);
}

double
IAPWS1995Rev::cv_l_from_pT(double p, double T)
{
  double rho_l = 0; bool search_failed = true;
  rho_l_from_pT(p, T, rho_l, search_failed);
  return search_failed ? 0 : Cv(T, rho_l);
}

double
IAPWS1995Rev::cp_l_from_pT(double p, double T)
{
  double rho_l = 0; bool search_failed = true;
  rho_l_from_pT(p, T, rho_l, search_failed);
  return search_failed ? 0 : Cp(T, rho_l);
}

double
IAPWS1995Rev::c_l_from_pT(double p, double T)
{
  double rho_l = 0; bool search_failed = true;
  rho_l_from_pT(p, T, rho_l, search_failed);
  return search_failed ? 0 : SoundSpeed(T, rho_l);
}

double
IAPWS1995Rev::k_l_from_pT(double p, double T)
{
  double rho_l = 0; bool search_failed = true;
  rho_l_from_pT(p, T, rho_l, search_failed);
  return search_failed ? 0 : thermal_conductivity_no_enhancement(T, rho_l);
}

double
IAPWS1995Rev::mu_l_from_pT(double p, double T)
{
  double rho_l = 0; bool search_failed = true;
  rho_l_from_pT(p, T, rho_l, search_failed);
  return search_failed ? 0 : viscosity(T, rho_l);
}

void
IAPWS1995Rev::liquid_properties_from_pT(double p, double T,
              double &rho, double &e, double &h, double &s, double &cv, double &cp,
              double &c, double &k, double &mu)
{
  bool search_failed = true;
  rho_l_from_pT(p, T, rho, search_failed);

  if (search_failed)  rho = 0.0;
  e   = search_failed ? 0 : IntEnergy(T, rho);
  h   = search_failed ? 0 : Enthalpy(T, rho);
  s   = search_failed ? 0 : Entropy(T, rho);
  cv  = search_failed ? 0 : Cv(T, rho);
  cp  = search_failed ? 0 : Cp(T, rho);
  c   = search_failed ? 0 : SoundSpeed(T, rho);
  k   = search_failed ? 0 : thermal_conductivity_no_enhancement(rho, T);
  mu  = search_failed ? 0 : viscosity(rho, T);
}


void
IAPWS1995Rev::rho_g_from_pT(double p, double T, double &rho, bool &search_failed)
{
  double Tsat = Tsat_from_p(p);
  double rho_l_sat = rho_l_sat_from_T(Tsat);
  double rho_g_sat = rho_g_sat_from_T(Tsat);

  double rho_min, rho_max, rho_guess;
  int it = 0;           int it_max = 1000;
  double rtol = 1e-9;   double error = 1.0;

  // For stable vapor region, using bisection iteration
  // For metastable vapor region, using IF97 formula
  if (T > Tsat)  // superheated (stable)
  {
    rho_min = 0.0; // an arbitrary small value for min value
    rho_max = rho_g_sat + 0.5; // tolerance

    while (it < it_max)
    {
      rho_guess = 0.5 * (rho_min + rho_max);
      double p_guess = Pressure(T, rho_guess);
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
    {
      std::cerr << "Max it number reached in IAPWS1995Rev::rho_g_from_pT\n";
      rho = 0;
      search_failed = true;
    }
    else
    {
      rho = rho_guess;
      search_failed = false;
    }
  }
  else  // subcooled (metastable)
  {
    double rho_g = 1.0 / R2Meta_specific_volume(p, T);
    if (rho_g < 0)
    {
      rho = 0;
      std::cerr << "rho_g < 0\n";
      search_failed = true;
    }
    else
    {
      rho = rho_g;
      search_failed = false;
    }
  }
}

double
IAPWS1995Rev::e_g_from_pT(double p, double T)
{
  double rho_g = 0; bool search_failed = true;
  rho_g_from_pT(p, T, rho_g, search_failed);
  return search_failed ? 0 : IntEnergy(T, rho_g);
}

double
IAPWS1995Rev::h_g_from_pT(double p, double T)
{
  double rho_g = 0; bool search_failed = true;
  rho_g_from_pT(p, T, rho_g, search_failed);
  return search_failed ? 0 : Enthalpy(T, rho_g);
}

double
IAPWS1995Rev::s_g_from_pT(double p, double T)
{
  double rho_g = 0; bool search_failed = true;
  rho_g_from_pT(p, T, rho_g, search_failed);
  return search_failed ? 0 : Entropy(T, rho_g);
}

double
IAPWS1995Rev::cv_g_from_pT(double p, double T)
{
  double rho_g = 0; bool search_failed = true;
  rho_g_from_pT(p, T, rho_g, search_failed);
  return search_failed ? 0 : Cv(T, rho_g);
}

double
IAPWS1995Rev::cp_g_from_pT(double p, double T)
{
  double rho_g = 0; bool search_failed = true;
  rho_g_from_pT(p, T, rho_g, search_failed);
  return search_failed ? 0 : Cp(T, rho_g);
}

double
IAPWS1995Rev::c_g_from_pT(double p, double T)
{
  double rho_g = 0; bool search_failed = true;
  rho_g_from_pT(p, T, rho_g, search_failed);
  return search_failed ? 0 : SoundSpeed(T, rho_g);
}

double
IAPWS1995Rev::k_g_from_pT(double p, double T)
{
  double rho_g = 0; bool search_failed = true;
  rho_g_from_pT(p, T, rho_g, search_failed);
  return search_failed ? 0 : thermal_conductivity_no_enhancement(rho_g, T);
}

double
IAPWS1995Rev::mu_g_from_pT(double p, double T)
{
  double rho_g = 0; bool search_failed = true;
  rho_g_from_pT(p, T, rho_g, search_failed);
  return search_failed ? 0 : viscosity(T, rho_g);
}

void
IAPWS1995Rev::vapor_properties_from_pT(double p, double T,
              double &rho, double &e, double &h, double &s, double &cv, double &cp,
              double &c, double &k, double &mu)
{
  bool search_failed = true;
  rho_g_from_pT(p, T, rho, search_failed);

  if (search_failed)  rho = 0.0;
  e   = search_failed ? 0 : IntEnergy(T, rho);
  h   = search_failed ? 0 : Enthalpy(T, rho);
  s   = search_failed ? 0 : Entropy(T, rho);
  cv  = search_failed ? 0 : Cv(T, rho);
  cp  = search_failed ? 0 : Cp(T, rho);
  c   = search_failed ? 0 : SoundSpeed(T, rho);
  k   = search_failed ? 0 : thermal_conductivity_no_enhancement(rho, T);
  mu  = search_failed ? 0 : viscosity(rho, T);
}
