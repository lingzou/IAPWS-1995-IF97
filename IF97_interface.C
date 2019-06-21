#include "IF97.h"
#include "IF97_helper.h"
#include "SurfaceTension.h"
#include "Viscosity.h"
#include "ThermalConductivity.h"

#include "IF97_interface.h"

double T_sat_from_p(double p)
{
  checkPSatValid(p);
  return R4_T_sat_from_p(p);
}

double p_sat_from_T(double T)
{
  checkTSatValid(T);
  return R4_p_sat_from_T(T);
}

double v_l_sat_from_T(double T)
{
  return 1.0 / rho_l_sat_from_T(T); /* This will check if T is valid */
}

double rho_l_sat_from_T(double T)
{
  checkTSatValid(T);

  if (T <= IF97_T_13)
    return 1.0 / R1_specific_volume(R4_p_sat_from_T(T), T);
  else
    return R3_rho_l_sat_from_T_ITER(T);
}

double e_l_sat_from_T(double T)
{
  checkTSatValid(T);

  if (T <= IF97_T_13)
    return R1_specific_int_energy(R4_p_sat_from_T(T), T);
  else
    return R3_specific_int_energy(R3_rho_l_sat_from_T_ITER(T), T);
}

double h_l_sat_from_T(double T)
{
  checkTSatValid(T);

  if (T <= IF97_T_13)
    return R1_specific_enthalpy(R4_p_sat_from_T(T), T);
  else
    return R3_specific_enthalpy(R3_rho_l_sat_from_T_ITER(T), T);
}

double s_l_sat_from_T(double T)
{
  checkTSatValid(T);

  if (T <= IF97_T_13)
    return R1_specific_entropy(R4_p_sat_from_T(T), T);
  else
    return R3_specific_entropy(R3_rho_l_sat_from_T_ITER(T), T);
}

double cv_l_sat_from_T(double T)
{
  checkTSatValid(T);

  if (T <= IF97_T_13)
    return R1_cv(R4_p_sat_from_T(T), T);
  else
    return R3_cv(R3_rho_l_sat_from_T_ITER(T), T);
}

double cp_l_sat_from_T(double T)
{
  checkTSatValid(T);

  if (T <= IF97_T_13)
    return R1_cp(R4_p_sat_from_T(T), T);
  else
    return R3_cp(R3_rho_l_sat_from_T_ITER(T), T);
}

double c_l_sat_from_T(double T)
{
  checkTSatValid(T);

  if (T <= IF97_T_13)
    return R1_sound_speed(R4_p_sat_from_T(T), T);
  else
    return R3_sound_speed(R3_rho_l_sat_from_T_ITER(T), T);
}

double k_l_sat_from_T(double T)
{
  checkTSatValid(T);

  if (T <= IF97_T_13)
    return thermal_conductivity_R1(R4_p_sat_from_T(T), T);
  else
    return thermal_conductivity_R3(R3_rho_l_sat_from_T_ITER(T), T);
}

double mu_l_sat_from_T(double T)
{
  double rho_l_sat = rho_l_sat_from_T(T); /* This will check if T is valid */
  return viscosity(rho_l_sat, T);
}

void liquid_sat_properties_from_T(double T, double & p, double & v, double & rho, double & e, double & h,
                                  double & s, double & cv, double & cp, double & c, double & k, double & mu)
{
  checkTSatValid(T);
  p = R4_p_sat_from_T(T);
  rho = rho_l_sat_from_T(T);
  v = 1.0 / rho;
  mu = viscosity(rho, T);

  if (T <= IF97_T_13)
  {
    e = R1_specific_int_energy(p, T);
    h = R1_specific_enthalpy(p, T);
    s = R1_specific_entropy(p, T);
    cv = R1_cv(p, T);
    cp = R1_cp(p, T);
    c = R1_sound_speed(p, T);
    k = thermal_conductivity_R1(p, T);
  }
  else
  {
    e = R3_specific_int_energy(rho, T);
    h = R3_specific_enthalpy(rho, T);
    s = R3_specific_entropy(rho, T);
    cv = R3_cv(rho, T);
    cp = R3_cp(rho, T);
    c = R3_sound_speed(rho, T);
    k = thermal_conductivity_R3(rho, T);
  }
}

double v_g_sat_from_T(double T)
{
  return 1.0 / rho_g_sat_from_T(T); /* This will check if T is valid */
}

double rho_g_sat_from_T(double T)
{
  checkTSatValid(T);

  if (T <= 623.15)
    return 1.0 / R2_specific_volume(R4_p_sat_from_T(T), T);
  else
    return R3_rho_g_sat_from_T_ITER(T);
}

double e_g_sat_from_T(double T)
{
  checkTSatValid(T);

  if (T <= IF97_T_13)
    return R2_specific_int_energy(R4_p_sat_from_T(T), T);
  else
    return R3_specific_int_energy(R3_rho_g_sat_from_T_ITER(T), T);
}

double h_g_sat_from_T(double T)
{
  checkTSatValid(T);

  if (T <= IF97_T_13)
    return R2_specific_enthalpy(R4_p_sat_from_T(T), T);
  else
    return R3_specific_enthalpy(R3_rho_g_sat_from_T_ITER(T), T);
}

double s_g_sat_from_T(double T)
{
  checkTSatValid(T);

  if (T <= IF97_T_13)
    return R2_specific_entropy(R4_p_sat_from_T(T), T);
  else
    return R3_specific_entropy(R3_rho_g_sat_from_T_ITER(T), T);
}

double cv_g_sat_from_T(double T)
{
  checkTSatValid(T);

  if (T <= IF97_T_13)
    return R2_cv(R4_p_sat_from_T(T), T);
  else
    return R3_cv(R3_rho_g_sat_from_T_ITER(T), T);
}

double cp_g_sat_from_T(double T)
{
  checkTSatValid(T);

  if (T <= IF97_T_13)
    return R2_cp(R4_p_sat_from_T(T), T);
  else
    return R3_cp(R3_rho_g_sat_from_T_ITER(T), T);
}

double c_g_sat_from_T(double T)
{
  checkTSatValid(T);

  if (T <= IF97_T_13)
    return R2_sound_speed(R4_p_sat_from_T(T), T);
  else
    return R3_sound_speed(R3_rho_g_sat_from_T_ITER(T), T);
}

double k_g_sat_from_T(double T)
{
  checkTSatValid(T);

  if (T <= IF97_T_13)
    return thermal_conductivity_R2(R4_p_sat_from_T(T), T);
  else
    return thermal_conductivity_R3(R3_rho_g_sat_from_T_ITER(T), T);
}

double mu_g_sat_from_T(double T)
{
  double rho_g_sat = rho_g_sat_from_T(T); /* This will check if T is valid */
  return viscosity(rho_g_sat, T);
}

void vapor_sat_properties_from_T(double T, double & p, double & v, double & rho, double & e, double & h,
                                  double & s, double & cv, double & cp, double & c, double & k, double & mu)
{
  checkTSatValid(T);
  p = R4_p_sat_from_T(T);
  rho = rho_g_sat_from_T(T);
  v = 1.0 / rho;
  mu = viscosity(rho, T);

  if (T <= IF97_T_13)
  {
    e = R2_specific_int_energy(p, T);
    h = R2_specific_enthalpy(p, T);
    s = R2_specific_entropy(p, T);
    cv = R2_cv(p, T);
    cp = R2_cp(p, T);
    c = R2_sound_speed(p, T);
    k = thermal_conductivity_R2(p, T);
  }
  else
  {
    e = R3_specific_int_energy(rho, T);
    h = R3_specific_enthalpy(rho, T);
    s = R3_specific_entropy(rho, T);
    cv = R3_cv(rho, T);
    cp = R3_cp(rho, T);
    c = R3_sound_speed(rho, T);
    k = thermal_conductivity_R3(rho, T);
  }
}

/***************************************************************
 * Saturation properties by p
 ***************************************************************/
double v_l_sat_from_p(double p)
{
  double T_sat = T_sat_from_p(p); // This will check if p is valid
  return v_l_sat_from_T(T_sat);
}

double rho_l_sat_from_p(double p)
{
  return rho_l_sat_from_T(T_sat_from_p(p));
}

double e_l_sat_from_p(double p)
{
  return e_l_sat_from_T(T_sat_from_p(p));
}

double h_l_sat_from_p(double p)
{
  return h_l_sat_from_T(T_sat_from_p(p));
}

double s_l_sat_from_p(double p)
{
  return s_l_sat_from_T(T_sat_from_p(p));
}

double cv_l_sat_from_p(double p)
{
  return cv_l_sat_from_T(T_sat_from_p(p));
}

double cp_l_sat_from_p(double p)
{
  return cp_l_sat_from_T(T_sat_from_p(p));
}

double c_l_sat_from_p(double p)
{
  return c_l_sat_from_T(T_sat_from_p(p));
}

double k_l_sat_from_p(double p)
{
  return k_l_sat_from_T(T_sat_from_p(p));
}

double mu_l_sat_from_p(double p)
{
  return mu_l_sat_from_T(T_sat_from_p(p));
}

void liquid_sat_properties_from_p(double p, double & T, double & v, double & rho, double & e, double & h,
                                  double & s, double & cv, double & cp, double & c, double & k, double & mu)
{
  T = T_sat_from_p(p); // This will check if p is valid
  rho = rho_l_sat_from_T(T);
  v = 1.0 / rho;
  mu = viscosity(rho, T);

  if (T <= IF97_T_13)
  {
    e = R1_specific_int_energy(p, T);
    h = R1_specific_enthalpy(p, T);
    s = R1_specific_entropy(p, T);
    cv = R1_cv(p, T);
    cp = R1_cp(p, T);
    c = R1_sound_speed(p, T);
    k = thermal_conductivity_R1(p, T);
  }
  else
  {
    e = R3_specific_int_energy(rho, T);
    h = R3_specific_enthalpy(rho, T);
    s = R3_specific_entropy(rho, T);
    cv = R3_cv(rho, T);
    cp = R3_cp(rho, T);
    c = R3_sound_speed(rho, T);
    k = thermal_conductivity_R3(rho, T);
  }
}

double v_g_sat_from_p(double p)
{
  return v_g_sat_from_T(T_sat_from_p(p));
}

double rho_g_sat_from_p(double p)
{
  return rho_g_sat_from_T(T_sat_from_p(p));
}

double e_g_sat_from_p(double p)
{
  return e_g_sat_from_T(T_sat_from_p(p));
}

double h_g_sat_from_p(double p)
{
  return h_g_sat_from_T(T_sat_from_p(p));
}

double s_g_sat_from_p(double p)
{
  return s_g_sat_from_T(T_sat_from_p(p));
}

double cv_g_sat_from_p(double p)
{
  return cv_g_sat_from_T(T_sat_from_p(p));
}

double cp_g_sat_from_p(double p)
{
  return cp_g_sat_from_T(T_sat_from_p(p));
}

double c_g_sat_from_p(double p)
{
  return c_g_sat_from_T(T_sat_from_p(p));
}

double k_g_sat_from_p(double p)
{
  return k_g_sat_from_T(T_sat_from_p(p));
}

double mu_g_sat_from_p(double p)
{
  return mu_g_sat_from_T(T_sat_from_p(p));
}

void vapor_sat_properties_from_p(double p, double & T, double & v, double & rho, double & e, double & h,
                                  double & s, double & cv, double & cp, double & c, double & k, double & mu)
{
  T = T_sat_from_p(p);
  rho = rho_g_sat_from_T(T);
  v = 1.0 / rho;
  mu = viscosity(rho, T);

  if (T <= IF97_T_13)
  {
    e = R2_specific_int_energy(p, T);
    h = R2_specific_enthalpy(p, T);
    s = R2_specific_entropy(p, T);
    cv = R2_cv(p, T);
    cp = R2_cp(p, T);
    c = R2_sound_speed(p, T);
    k = thermal_conductivity_R2(p, T);
  }
  else
  {
    e = R3_specific_int_energy(rho, T);
    h = R3_specific_enthalpy(rho, T);
    s = R3_specific_entropy(rho, T);
    cv = R3_cv(rho, T);
    cp = R3_cp(rho, T);
    c = R3_sound_speed(rho, T);
    k = thermal_conductivity_R3(rho, T);
  }
}

/***************************************************************
 * (p, T)-based properties
 ***************************************************************/
double v_from_pT(double p, double T)
{
  int region = locateRegion_from_pT(p, T);
  switch (region) {
    case 1:
      return R1_specific_volume(p, T);
    case 2:
      return R2_specific_volume(p, T);
    case 3:
      return 1.0 / R3_rho_from_p_T_ITER(p, T);
    case 5:
      return R5_specific_volume(p, T);
    default:
      fprintf(stderr, "%s", "Region not recognized!\n");
      exit(1);
      return 0.0;
  }
  //return (*v_func_from_pT_ptr[region - 1])(p, T);
}

double rho_from_pT(double p, double T)
{
  int region = locateRegion_from_pT(p, T);
  switch (region) {
    case 1:
      return 1.0 / R1_specific_volume(p, T);
    case 2:
      return 1.0 / R2_specific_volume(p, T);
    case 3:
      return R3_rho_from_p_T_ITER(p, T);
    case 5:
      return 1.0 / R5_specific_volume(p, T);
    default:
      fprintf(stderr, "%s", "Region not recognized!\n");
      exit(1);
      return 0.0;
  }
  //return (*rho_func_from_pT_ptr[region - 1])(p, T);
}


double e_from_pT(double p, double T)
{
  int region = locateRegion_from_pT(p, T);
  switch (region) {
    case 1:
      return R1_specific_int_energy(p, T);
    case 2:
      return R2_specific_int_energy(p, T);
    case 3:
      return R3_specific_int_energy(R3_rho_from_p_T_ITER(p, T), T);
    case 5:
      return R5_specific_int_energy(p, T);
    default:
      fprintf(stderr, "%s", "Region not recognized!\n");
      exit(1);
      return 0.0;
  }
}

double h_from_pT(double p, double T)
{
  int region = locateRegion_from_pT(p, T);
  switch (region) {
    case 1:
      return R1_specific_enthalpy(p, T);
    case 2:
      return R2_specific_enthalpy(p, T);
    case 3:
      return R3_specific_enthalpy(R3_rho_from_p_T_ITER(p, T), T);
    case 5:
      return R5_specific_enthalpy(p, T);
    default:
      fprintf(stderr, "%s", "Region not recognized!\n");
      exit(1);
      return 0.0;
  }
}

double s_from_pT(double p, double T)
{
  int region = locateRegion_from_pT(p, T);
  switch (region) {
    case 1:
      return R1_specific_entropy(p, T);
    case 2:
      return R2_specific_entropy(p, T);
    case 3:
      return R3_specific_entropy(R3_rho_from_p_T_ITER(p, T), T);
    case 5:
      return R5_specific_entropy(p, T);
    default:
      fprintf(stderr, "%s", "Region not recognized!\n");
      exit(1);
      return 0.0;
  }
}

double cv_from_pT(double p, double T)
{
  int region = locateRegion_from_pT(p, T);
  switch (region) {
    case 1:
      return R1_cv(p, T);
    case 2:
      return R2_cv(p, T);
    case 3:
      return R3_cv(R3_rho_from_p_T_ITER(p, T), T);
    case 5:
      return R5_cv(p, T);
    default:
      fprintf(stderr, "%s", "Region not recognized!\n");
      exit(1);
      return 0.0;
  }
}

double cp_from_pT(double p, double T)
{
  int region = locateRegion_from_pT(p, T);
  switch (region) {
    case 1:
      return R1_cp(p, T);
    case 2:
      return R2_cp(p, T);
    case 3:
      return R3_cp(R3_rho_from_p_T_ITER(p, T), T);
    case 5:
      return R5_cp(p, T);
    default:
      fprintf(stderr, "%s", "Region not recognized!\n");
      exit(1);
      return 0.0;
  }
}

double c_from_pT(double p, double T)
{
  int region = locateRegion_from_pT(p, T);
  switch (region) {
    case 1:
      return R1_sound_speed(p, T);
    case 2:
      return R2_sound_speed(p, T);
    case 3:
      return R3_sound_speed(R3_rho_from_p_T_ITER(p, T), T);
    case 5:
      return R5_sound_speed(p, T);
    default:
      fprintf(stderr, "%s", "Region not recognized!\n");
      exit(1);
      return 0.0;
  }
}

double k_from_pT(double p, double T)
{
  int region = locateRegion_from_pT(p, T);
  switch (region) {
    case 1:
      return thermal_conductivity_R1(p, T);
    case 2:
      return thermal_conductivity_R2(p, T);
    case 3:
      return thermal_conductivity_R3(R3_rho_from_p_T_ITER(p, T), T);
    case 5:
      return thermal_conductivity_R5(p, T);
    default:
      fprintf(stderr, "%s", "Region not recognized!\n");
      exit(1);
      return 0.0;
  }
}

double mu_from_pT(double p, double T)
{
  double rho = rho_from_pT(p, T);
  return viscosity(rho, T);
}
