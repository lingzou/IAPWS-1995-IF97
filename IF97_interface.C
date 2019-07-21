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

void liquid_sat_properties_from_T(double T, double * p, double * v, double * rho, double * e, double * h,
                                  double * s, double * cv, double * cp, double * c, double * k, double * mu)
{
  checkTSatValid(T);
  double p_val = R4_p_sat_from_T(T);
  double rho_val = rho_l_sat_from_T(T);
  if (p != NULL)   *p = p_val;
  if (v != NULL)   *v = 1.0 / rho_val;
  if (rho != NULL) *rho = rho_val;
  if (mu != NULL)  *mu = viscosity(rho_val, T);

  if (T <= IF97_T_13)
  {
    if (e != NULL)  *e = R1_specific_int_energy(p_val, T);
    if (h != NULL)  *h = R1_specific_enthalpy(p_val, T);
    if (s != NULL)  *s = R1_specific_entropy(p_val, T);
    if (cv != NULL) *cv = R1_cv(p_val, T);
    if (cp != NULL) *cp = R1_cp(p_val, T);
    if (c != NULL)  *c = R1_sound_speed(p_val, T);
    if (k != NULL)  *k = thermal_conductivity_R1(p_val, T);
  }
  else
  {
    if (e != NULL)  *e = R3_specific_int_energy(rho_val, T);
    if (h != NULL)  *h = R3_specific_enthalpy(rho_val, T);
    if (s != NULL)  *s = R3_specific_entropy(rho_val, T);
    if (cv != NULL) *cv = R3_cv(rho_val, T);
    if (cp != NULL) *cp = R3_cp(rho_val, T);
    if (c != NULL)  *c = R3_sound_speed(rho_val, T);
    if (k != NULL)  *k = thermal_conductivity_R3(rho_val, T);
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

void vapor_sat_properties_from_T(double T, double * p, double * v, double * rho, double * e, double * h,
                                  double * s, double * cv, double * cp, double * c, double * k, double * mu)
{
  checkTSatValid(T);
  double p_val = R4_p_sat_from_T(T);
  double rho_val = rho_g_sat_from_T(T);
  if (p != NULL)   *p = p_val;
  if (v != NULL)   *v = 1.0 / rho_val;
  if (rho != NULL) *rho = rho_val;
  if (mu != NULL)  *mu = viscosity(rho_val, T);

  if (T <= IF97_T_13)
  {
    if (e != NULL)  *e = R2_specific_int_energy(p_val, T);
    if (h != NULL)  *h = R2_specific_enthalpy(p_val, T);
    if (s != NULL)  *s = R2_specific_entropy(p_val, T);
    if (cv != NULL) *cv = R2_cv(p_val, T);
    if (cp != NULL) *cp = R2_cp(p_val, T);
    if (c != NULL)  *c = R2_sound_speed(p_val, T);
    if (k != NULL)  *k = thermal_conductivity_R2(p_val, T);
  }
  else
  {
    if (e != NULL)  *e = R3_specific_int_energy(rho_val, T);
    if (h != NULL)  *h = R3_specific_enthalpy(rho_val, T);
    if (s != NULL)  *s = R3_specific_entropy(rho_val, T);
    if (cv != NULL) *cv = R3_cv(rho_val, T);
    if (cp != NULL) *cp = R3_cp(rho_val, T);
    if (c != NULL)  *c = R3_sound_speed(rho_val, T);
    if (k != NULL)  *k = thermal_conductivity_R3(rho_val, T);
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

void liquid_sat_properties_from_p(double p, double * T, double * v, double * rho, double * e, double * h,
                                  double * s, double * cv, double * cp, double * c, double * k, double * mu)
{
  double T_val = T_sat_from_p(p); // This will check if p is valid
  double rho_val = rho_l_sat_from_T(T_val);
  if (T != NULL)   *T = T_val;
  if (v != NULL)   *v = 1.0 / rho_val;
  if (rho != NULL) *rho = rho_val;
  if (mu != NULL)  *mu = viscosity(rho_val, T_val);

  if (T_val <= IF97_T_13)
  {
    if (e != NULL)  *e = R1_specific_int_energy(p, T_val);
    if (h != NULL)  *h = R1_specific_enthalpy(p, T_val);
    if (s != NULL)  *s = R1_specific_entropy(p, T_val);
    if (cv != NULL) *cv = R1_cv(p, T_val);
    if (cp != NULL) *cp = R1_cp(p, T_val);
    if (c != NULL)  *c = R1_sound_speed(p, T_val);
    if (k != NULL)  *k = thermal_conductivity_R1(p, T_val);
  }
  else
  {
    if (e != NULL)  *e = R3_specific_int_energy(rho_val, T_val);
    if (h != NULL)  *h = R3_specific_enthalpy(rho_val, T_val);
    if (s != NULL)  *s = R3_specific_entropy(rho_val, T_val);
    if (cv != NULL) *cv = R3_cv(rho_val, T_val);
    if (cp != NULL) *cp = R3_cp(rho_val, T_val);
    if (c != NULL)  *c = R3_sound_speed(rho_val, T_val);
    if (k != NULL)  *k = thermal_conductivity_R3(rho_val, T_val);
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

void vapor_sat_properties_from_p(double p, double * T, double * v, double * rho, double * e, double * h,
                                  double * s, double * cv, double * cp, double * c, double * k, double * mu)
{
  double T_val = T_sat_from_p(p);
  double rho_val = rho_g_sat_from_T(T_val);
  if (T != NULL)   *T = T_val;
  if (v != NULL)   *v = 1.0 / rho_val;
  if (rho != NULL) *rho = rho_val;
  if (mu != NULL)  *mu = viscosity(rho_val, T_val);

  if (T_val <= IF97_T_13)
  {
    if (e != NULL)  *e = R2_specific_int_energy(p, T_val);
    if (h != NULL)  *h = R2_specific_enthalpy(p, T_val);
    if (s != NULL)  *s = R2_specific_entropy(p, T_val);
    if (cv != NULL) *cv = R2_cv(p, T_val);
    if (cp != NULL) *cp = R2_cp(p, T_val);
    if (c != NULL)  *c = R2_sound_speed(p, T_val);
    if (k != NULL)  *k = thermal_conductivity_R2(p, T_val);
  }
  else
  {
    if (e != NULL)  *e = R3_specific_int_energy(rho_val, T_val);
    if (h != NULL)  *h = R3_specific_enthalpy(rho_val, T_val);
    if (s != NULL)  *s = R3_specific_entropy(rho_val, T_val);
    if (cv != NULL) *cv = R3_cv(rho_val, T_val);
    if (cp != NULL) *cp = R3_cp(rho_val, T_val);
    if (c != NULL)  *c = R3_sound_speed(rho_val, T_val);
    if (k != NULL)  *k = thermal_conductivity_R3(rho_val, T_val);
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

void properties_from_pT(double p, double T, double * v, double * rho, double * e, double * h,
                         double * s, double * cv, double * cp, double * c, double * k, double * mu)
{
  int region = locateRegion_from_pT(p, T);
  switch (region) {
    case 1:
    {
      double v_val = R1_specific_volume(p, T);
      if (v != NULL)   *v = v_val;
      if (rho != NULL) *rho = 1.0 / v_val;
      if (e != NULL)   *e = R1_specific_int_energy(p, T);
      if (h != NULL)   *e = R1_specific_enthalpy(p, T);
      if (s != NULL)   *s = R1_specific_entropy(p, T);
      if (cv != NULL)  *cv = R1_cv(p, T);
      if (cp != NULL)  *cp = R1_cp(p, T);
      if (c != NULL)   *c = R1_sound_speed(p, T);
      if (k != NULL)   *k = thermal_conductivity_R1(p, T);
      if (mu != NULL)  *mu = mu_from_pT(p, T); // should have mu_from_rhoT
    }
    break;

    case 2:
    {
      double v_val = R2_specific_volume(p, T);
      if (v != NULL)   *v = v_val;
      if (rho != NULL) *rho = 1.0 / v_val;
      if (e != NULL)   *e = R2_specific_int_energy(p, T);
      if (h != NULL)   *e = R2_specific_enthalpy(p, T);
      if (s != NULL)   *s = R2_specific_entropy(p, T);
      if (cv != NULL)  *cv = R2_cv(p, T);
      if (cp != NULL)  *cp = R2_cp(p, T);
      if (c != NULL)   *c = R2_sound_speed(p, T);
      if (k != NULL)   *k = thermal_conductivity_R2(p, T);
      if (mu != NULL)  *mu = mu_from_pT(p, T); // should have mu_from_rhoT
    }
    break;

    case 3:
    {
      double rho_val = R3_rho_from_p_T_ITER(p, T);
      if (v != NULL)   *v = 1.0 / rho_val;
      if (rho != NULL) *rho = rho_val;
      if (e != NULL)   *e = R3_specific_int_energy(rho_val, T);
      if (h != NULL)   *e = R3_specific_enthalpy(rho_val, T);
      if (s != NULL)   *s = R3_specific_entropy(rho_val, T);
      if (cv != NULL)  *cv = R3_cv(rho_val, T);
      if (cp != NULL)  *cp = R3_cp(rho_val, T);
      if (c != NULL)   *c = R3_sound_speed(rho_val, T);
      if (k != NULL)   *k = thermal_conductivity_R3(rho_val, T);
      if (mu != NULL)  *mu = mu_from_pT(p, T); // should have mu_from_rhoT
    }
    break;

    case 5:
    {
      double v_val = R5_specific_volume(p, T);
      if (v != NULL)   *v = v_val;
      if (rho != NULL) *rho = 1.0 / v_val;
      if (e != NULL)   *e = R5_specific_int_energy(p, T);
      if (h != NULL)   *e = R5_specific_enthalpy(p, T);
      if (s != NULL)   *s = R5_specific_entropy(p, T);
      if (cv != NULL)  *cv = R5_cv(p, T);
      if (cp != NULL)  *cp = R5_cp(p, T);
      if (c != NULL)   *c = R5_sound_speed(p, T);
      if (k != NULL)   *k = thermal_conductivity_R5(p, T);
      if (mu != NULL)  *mu = mu_from_pT(p, T); // should have mu_from_rhoT
    }
    break;

    default:
      fprintf(stderr, "%s", "Region not recognized!\n");
      exit(1);
  }
}
/***************************************************************
 * (p, h)-based properties
 ***************************************************************/
int locateRegion_from_ph(double p, double h)
{
  // 21 -> 2a; 22 -> 2b; 23 -> 2c.
  if (p < IF97_SAT_P_MIN)
  {
    if (h < R2_specific_enthalpy(p, IF97_T_MIN))
    {
      fprintf(stderr, "p = %f; h = %f\n", p, h);
      fprintf(stderr, "Out of range: h < h(p, 273.15K)!\n");
      return -1;
    }
    else if (h < R2_specific_enthalpy(p, IF97_T_25))
      return 21;
    else if (h < R5_specific_enthalpy(p, IF97_T_MAX))
      return 5;
    else
    {
      fprintf(stderr, "p = %f; h = %f\n", p, h);
      fprintf(stderr, "Out of range: h > h(p, 2273.15K)!\n");
      return -2;
    }
  }
  else if (p < 4.0e6) // boundary between 2a and 2b (see page 21, Ref. [1])
  {
    if (h < R1_specific_enthalpy(p, IF97_T_MIN))
    {
      fprintf(stderr, "p = %f; h = %f\n", p, h);
      fprintf(stderr, "Out of range: h < h(p, 273.15K)!\n");
      return -1;
    }
    else if (h < h_l_sat_from_p(p))
      return 1;
    else if (h < h_g_sat_from_p(p))
      return 4;
    else if (h < R2_specific_enthalpy(p, IF97_T_25))
      return 21;
    else if (h < R5_specific_enthalpy(p, IF97_T_MAX))
      return 5;
    else
    {
      fprintf(stderr, "p = %f; h = %f\n", p, h);
      fprintf(stderr, "Out of range: h > h(p, 2273.15K)!\n");
      return -2;
    }
  }
  else if (p < 6.5467e6) // where 2b-2c boundary line intersects with saturation line, (see page 21, Ref. [1])
  {
    if (h < R1_specific_enthalpy(p, IF97_T_MIN))
    {
      fprintf(stderr, "p = %f; h = %f\n", p, h);
      fprintf(stderr, "Out of range: h < h(p, 273.15K)!\n");
      return -1;
    }
    else if (h < h_l_sat_from_p(p))
      return 1;
    else if (h < h_g_sat_from_p(p))
      return 4;
    else if (h < R2_specific_enthalpy(p, IF97_T_25))
      return 22;
    else if (h < R5_specific_enthalpy(p, IF97_T_MAX))
      return 5;
    else
    {
      fprintf(stderr, "p = %f; h = %f\n", p, h);
      fprintf(stderr, "Out of range: h > h(p, 2273.15K)!\n");
      return -2;
    }
  }
  else if (p < 1.65291643e7) // where 2-3 boundary line intersects with saturation line, (see page 6, Ref. [1])
  {
    if (h < R1_specific_enthalpy(p, IF97_T_MIN))
    {
      fprintf(stderr, "p = %f; h = %f\n", p, h);
      fprintf(stderr, "Out of range: h < h(p, 273.15K)!\n");
      return -1;
    }
    else if (h < h_l_sat_from_p(p))
      return 1;
    else if (h < h_g_sat_from_p(p))
      return 4;
    else if (h < B2bc_h_from_p(p))
      return 23;
    else if (h < R2_specific_enthalpy(p, IF97_T_25))
      return 22;
    else if (h < R5_specific_enthalpy(p, IF97_T_MAX))
      return 5;
    else
    {
      fprintf(stderr, "p = %f; h = %f\n", p, h);
      fprintf(stderr, "Out of range: h > h(p, 2273.15K)!\n");
      return -2;
    }
  }
  else if (p < P_CRIT)
  {
    if (h < R1_specific_enthalpy(p, IF97_T_MIN))
    {
      fprintf(stderr, "p = %f; h = %f\n", p, h);
      fprintf(stderr, "Out of range: h < h(p, 273.15K)!\n");
      return -1;
    }
    else if (h < R1_specific_enthalpy(p, IF97_T_13))
      return 1;
    else if (h < h_l_sat_from_p(p))
      return 3;
    else if (h < h_g_sat_from_p(p))
      return 4;
    else if (h < R2_specific_enthalpy(p, B23_T_from_p(p)))
      return 3;
    else if (h < B2bc_h_from_p(p))
      return 23;
    else if (h < R2_specific_enthalpy(p, IF97_T_25))
      return 22;
    else if (h < R5_specific_enthalpy(p, IF97_T_MAX))
      return 5;
    else
    {
      fprintf(stderr, "p = %f; h = %f\n", p, h);
      fprintf(stderr, "Out of range: h > h(p, 2273.15K)!\n");
      return -2;
    }
  }
  else if (p <= 50.0e6)
  {
    if (h < R1_specific_enthalpy(p, IF97_T_MIN))
    {
      fprintf(stderr, "p = %f; h = %f\n", p, h);
      fprintf(stderr, "Out of range: h < h(p, 273.15K)!\n");
      return -1;
    }
    else if (h < R1_specific_enthalpy(p, IF97_T_13))
      return 1;
    else if (h < R2_specific_enthalpy(p, B23_T_from_p(p)))
      return 3;
    else if (h < B2bc_h_from_p(p))
      return 23;
    else if (h < R2_specific_enthalpy(p, IF97_T_25))
      return 22;
    else if (h < R5_specific_enthalpy(p, IF97_T_MAX))
      return 5;
    else
    {
      fprintf(stderr, "p = %f; h = %f\n", p, h);
      fprintf(stderr, "Out of range: h > h(p, 2273.15K)!\n");
      return -2;
    }
  }
  else if (p <= IF97_P_MAX)
  {
    if (h < R1_specific_enthalpy(p, IF97_T_MIN))
    {
      fprintf(stderr, "p = %f; h = %f\n", p, h);
      fprintf(stderr, "Out of range: h < h(p, 273.15K)!\n");
      return -1;
    }
    else if (h < R1_specific_enthalpy(p, IF97_T_13))
      return 1;
    else if (h < R2_specific_enthalpy(p, B23_T_from_p(p)))
      return 3;
    else if (h < B2bc_h_from_p(p))
      return 23;
    else if (h < R2_specific_enthalpy(p, IF97_T_25))
      return 22;
    else
    {
      fprintf(stderr, "p = %f; h = %f\n", p, h);
      fprintf(stderr, "Out of range: h > h(p, 1073.15K)!\n");
      return -4;
    }
  }
  else
  {
    return -3;
  }
}

double v_from_ph(double p, double h)
{
  int region = locateRegion_from_ph(p, h);
  switch (region) {
    case 1:
      return R1_specific_volume(p, R1_T_from_p_h(p, h));
    case 21:
      return R2_specific_volume(p, R2a_T_from_p_h(p, h));
    case 22:
      return R2_specific_volume(p, R2b_T_from_p_h(p, h));
    case 23:
      return R2_specific_volume(p, R2c_T_from_p_h(p, h));
    case 3:
    {
      double rho = 0.0, T = 0.0, x = 0.0;
      R3_rho_T_x_from_p_h_ITER(p, h, rho, T, x); /* This function covers Region 3 and the subsetion of 4 within 3 */
      return 1.0 / rho;
    }
    case 4:
    {
      double v_l_sat = 1.0 / rho_l_sat_from_p(p);
      double v_g_sat = 1.0 / rho_g_sat_from_p(p);
      double h_l_sat = h_l_sat_from_p(p);
      double h_g_sat = h_g_sat_from_p(p);
      double x = (h - h_l_sat) / (h_g_sat - h_l_sat);
      return x * v_g_sat + (1.0 - x) * v_l_sat;
    }
    case 5:
      return R5_specific_volume(p, R5_T_from_p_h_ITER(p, h));
    default:
      fprintf(stderr, "%s", "Region not recognized!\n");
      exit(1);
      return 0.0;
  }
}

double rho_from_ph(double p, double h)
{
  int region = locateRegion_from_ph(p, h);
  switch (region) {
    case 1:
      return 1.0 / R1_specific_volume(p, R1_T_from_p_h(p, h));
    case 21:
      return 1.0 / R2_specific_volume(p, R2a_T_from_p_h(p, h));
    case 22:
      return 1.0 / R2_specific_volume(p, R2b_T_from_p_h(p, h));
    case 23:
      return 1.0 / R2_specific_volume(p, R2c_T_from_p_h(p, h));
    case 3:
    {
      double rho = 0.0, T = 0.0, x = 0.0;
      R3_rho_T_x_from_p_h_ITER(p, h, rho, T, x); /* This function covers Region 3 and the subsetion of 4 within 3 */
      return rho;
    }
    case 4:
    {
      double rho_l_sat = rho_l_sat_from_p(p);
      double rho_g_sat = rho_g_sat_from_p(p);
      double h_l_sat = h_l_sat_from_p(p);
      double h_g_sat = h_g_sat_from_p(p);
      double x = (h - h_l_sat) / (h_g_sat - h_l_sat);
      return 1.0 / ((1.0 - x) / rho_l_sat + x / rho_g_sat);
    }
    case 5:
      return 1.0 / R5_specific_volume(p, R5_T_from_p_h_ITER(p, h));
    default:
      fprintf(stderr, "%s", "Region not recognized!\n");
      exit(1);
      return 0.0;
  }
}

double e_from_ph(double p, double h)
{
  double v = v_from_ph(p, h);
  return h - p * v; // Definition, works for both single-phase and equilibrium two-phase region
}

double T_from_ph(double p, double h)
{
  int region = locateRegion_from_ph(p, h);
  switch (region) {
    case 1:
      return R1_T_from_p_h(p, h);
    case 21:
      return R2a_T_from_p_h(p, h);
    case 22:
      return R2b_T_from_p_h(p, h);
    case 23:
      return R2c_T_from_p_h(p, h);
    case 3:
    {
      double rho = 0.0, T = 0.0, x = 0.0;
      R3_rho_T_x_from_p_h_ITER(p, h, rho, T, x); /* This function covers Region 3 and the subsetion of 4 within 3 */
      return T;
    }
    case 4:
    {
      /*
      double h_l_sat = h_l_sat_from_p(p);
      double h_g_sat = h_g_sat_from_p(p);
      double x = (h - h_l_sat) / (h_g_sat - h_l_sat);*/
      return T_sat_from_p(p);
    }
    case 5:
      return R5_T_from_p_h_ITER(p, h);
    default:
      fprintf(stderr, "%s", "Region not recognized!\n");
      exit(1);
      return 0.0;
  }
}

double s_from_ph(double p, double h)
{
  int region = locateRegion_from_ph(p, h);
  switch (region) {
    case 1:
      return R1_specific_entropy(p, R1_T_from_p_h(p, h));
    case 21:
      return R2_specific_entropy(p, R2a_T_from_p_h(p, h));
    case 22:
      return R2_specific_entropy(p, R2b_T_from_p_h(p, h));
    case 23:
      return R2_specific_entropy(p, R2c_T_from_p_h(p, h));
    case 3:
    {
      double rho = 0.0, T = 0.0, x = 0.0;
      R3_rho_T_x_from_p_h_ITER(p, h, rho, T, x); /* This function covers Region 3 and the subsetion of 4 within 3 */
      return R3_specific_entropy(rho, T);
    }
    case 4:
    {
      double h_l_sat = h_l_sat_from_p(p);
      double h_g_sat = h_g_sat_from_p(p);
      double s_l_sat = s_l_sat_from_p(p);
      double s_g_sat = s_g_sat_from_p(p);
      double x = (h - h_l_sat) / (h_g_sat - h_l_sat);
      return x * s_g_sat + (1.0 - x) * s_l_sat; // Check correctness
    }
    case 5:
      return R5_specific_entropy(p, R5_T_from_p_h_ITER(p, h));
    default:
      fprintf(stderr, "%s", "Region not recognized!\n");
      exit(1);
      return 0.0;
  }
}

double cv_from_ph(double p, double h)
{
  int region = locateRegion_from_ph(p, h);
  switch (region) {
    case 1:
      return R1_cv(p, R1_T_from_p_h(p, h));
    case 21:
      return R2_cv(p, R2a_T_from_p_h(p, h));
    case 22:
      return R2_cv(p, R2b_T_from_p_h(p, h));
    case 23:
      return R2_cv(p, R2c_T_from_p_h(p, h));
    case 3:
    {
      double rho = 0.0, T = 0.0, x = 0.0;
      R3_rho_T_x_from_p_h_ITER(p, h, rho, T, x); /* This function covers Region 3 and the subsetion of 4 within 3 */
      return R3_cv(rho, T);
    }
    case 4:
    {
      return -1.0; // not defined
    }
    case 5:
      return R5_cv(p, R5_T_from_p_h_ITER(p, h));
    default:
      fprintf(stderr, "%s", "Region not recognized!\n");
      exit(1);
      return 0.0;
  }
}

double cp_from_ph(double p, double h)
{
  int region = locateRegion_from_ph(p, h);
  switch (region) {
    case 1:
      return R1_cp(p, R1_T_from_p_h(p, h));
    case 21:
      return R2_cp(p, R2a_T_from_p_h(p, h));
    case 22:
      return R2_cp(p, R2b_T_from_p_h(p, h));
    case 23:
      return R2_cp(p, R2c_T_from_p_h(p, h));
    case 3:
    {
      double rho = 0.0, T = 0.0, x = 0.0;
      R3_rho_T_x_from_p_h_ITER(p, h, rho, T, x); /* This function covers Region 3 and the subsetion of 4 within 3 */
      return R3_cp(rho, T);
    }
    case 4:
    {
      return -1.0; // not defined
    }
    case 5:
      return R5_cp(p, R5_T_from_p_h_ITER(p, h));
    default:
      fprintf(stderr, "%s", "Region not recognized!\n");
      exit(1);
      return 0.0;
  }
}

double c_from_ph(double p, double h)
{
  int region = locateRegion_from_ph(p, h);
  switch (region) {
    case 1:
      return R1_sound_speed(p, R1_T_from_p_h(p, h));
    case 21:
      return R2_sound_speed(p, R2a_T_from_p_h(p, h));
    case 22:
      return R2_sound_speed(p, R2b_T_from_p_h(p, h));
    case 23:
      return R2_sound_speed(p, R2c_T_from_p_h(p, h));
    case 3:
    {
      double rho = 0.0, T = 0.0, x = 0.0;
      R3_rho_T_x_from_p_h_ITER(p, h, rho, T, x); /* This function covers Region 3 and the subsetion of 4 within 3 */
      return R3_sound_speed(rho, T);
    }
    case 4:
    {
      return -1.0; // needs special model, not this package's interest
    }
    case 5:
      return R5_sound_speed(p, R5_T_from_p_h_ITER(p, h));
    default:
      fprintf(stderr, "%s", "Region not recognized!\n");
      exit(1);
      return 0.0;
  }
}

double k_from_ph(double p, double h)
{
  int region = locateRegion_from_ph(p, h);
  switch (region) {
    case 1:
      return thermal_conductivity_R1(p, R1_T_from_p_h(p, h));
    case 21:
      return thermal_conductivity_R2(p, R2a_T_from_p_h(p, h));
    case 22:
      return thermal_conductivity_R2(p, R2b_T_from_p_h(p, h));
    case 23:
      return thermal_conductivity_R2(p, R2c_T_from_p_h(p, h));
    case 3:
    {
      double rho = 0.0, T = 0.0, x = 0.0;
      R3_rho_T_x_from_p_h_ITER(p, h, rho, T, x); /* This function covers Region 3 and the subsetion of 4 within 3 */
      return thermal_conductivity_R3(rho, T);
    }
    case 4:
    {
      return -1.0; // not defined
    }
    case 5:
      return thermal_conductivity_R5(p, R5_T_from_p_h_ITER(p, h));
    default:
      fprintf(stderr, "%s", "Region not recognized!\n");
      exit(1);
      return 0.0;
  }
}

double mu_from_ph(double p, double h)
{
  int region = locateRegion_from_ph(p, h);
  switch (region) {
    case 1:
    case 21:
    case 22:
    case 23:
    case 3:
    case 5:
    {
      double rho = rho_from_ph(p, h);
      double T = T_from_ph(p, h);
      return viscosity(rho, T);
    }
    case 4:
    {
      return -1.0; // not defined
    }
    default:
      fprintf(stderr, "%s", "Region not recognized!\n");
      exit(1);
      return 0.0;
  }
}

void properties_from_ph(double p, double h, double * v, double * rho, double * e, double * T,
                         double * s, double * cv, double * cp, double * c, double * k, double * mu)
{
  int region = locateRegion_from_ph(p, h);

  switch (region) {
    case 1:
    {
      double T_val = R1_T_from_p_h(p, h);
      double v_val = R1_specific_volume(p, T_val);
      double rho_val = 1.0 / v_val;
      if (v != NULL)   *v = v_val;
      if (rho != NULL) *rho = rho_val;
      if (e != NULL)   *e = h - p * v_val;
      if (T != NULL)   *T = T_val;
      if (s != NULL)   *s = R1_specific_entropy(p, T_val);
      if (cv != NULL)  *cv = R1_cv(p, T_val);
      if (cp != NULL)  *cp = R1_cp(p, T_val);
      if (c != NULL)   *c = R1_sound_speed(p, T_val);
      if (k != NULL)   *k = thermal_conductivity_R1(p, T_val);
      if (mu != NULL)  *mu = viscosity(rho_val, T_val);
    }
    case 21:
    {
      double T_val = R2a_T_from_p_h(p, h);
      double v_val = R2_specific_volume(p, T_val);
      double rho_val = 1.0 / v_val;
      if (v != NULL)   *v = v_val;
      if (rho != NULL) *rho = rho_val;
      if (e != NULL)   *e = h - p * v_val;
      if (T != NULL)   *T = T_val;
      if (s != NULL)   *s = R2_specific_entropy(p, T_val);
      if (cv != NULL)  *cv = R2_cv(p, T_val);
      if (cp != NULL)  *cp = R2_cp(p, T_val);
      if (c != NULL)   *c = R2_sound_speed(p, T_val);
      if (k != NULL)   *k = thermal_conductivity_R2(p, T_val);
      if (mu != NULL)  *mu = viscosity(rho_val, T_val);
    }
    case 22:
    {
      double T_val = R2b_T_from_p_h(p, h);
      double v_val = R2_specific_volume(p, T_val);
      double rho_val = 1.0 / v_val;
      if (v != NULL)   *v = v_val;
      if (rho != NULL) *rho = rho_val;
      if (e != NULL)   *e = h - p * v_val;
      if (T != NULL)   *T = T_val;
      if (s != NULL)   *s = R2_specific_entropy(p, T_val);
      if (cv != NULL)  *cv = R2_cv(p, T_val);
      if (cp != NULL)  *cp = R2_cp(p, T_val);
      if (c != NULL)   *c = R2_sound_speed(p, T_val);
      if (k != NULL)   *k = thermal_conductivity_R2(p, T_val);
      if (mu != NULL)  *mu = viscosity(rho_val, T_val);
    }
    case 23:
    {
      double T_val = R2c_T_from_p_h(p, h);
      double v_val = R2_specific_volume(p, T_val);
      double rho_val = 1.0 / v_val;
      if (v != NULL)   *v = v_val;
      if (rho != NULL) *rho = rho_val;
      if (e != NULL)   *e = h - p * v_val;
      if (T != NULL)   *T = T_val;
      if (s != NULL)   *s = R2_specific_entropy(p, T_val);
      if (cv != NULL)  *cv = R2_cv(p, T_val);
      if (cp != NULL)  *cp = R2_cp(p, T_val);
      if (c != NULL)   *c = R2_sound_speed(p, T_val);
      if (k != NULL)   *k = thermal_conductivity_R2(p, T_val);
      if (mu != NULL)  *mu = viscosity(rho_val, T_val);
    }
    case 3:
    {
      double rho_val = 0.0, T_val = 0.0, x_val = 0.0;
      R3_rho_T_x_from_p_h_ITER(p, h, rho_val, T_val, x_val); /* This function covers Region 3 and the subsetion of 4 within 3 */
      if (v != NULL)   *v = 1.0 / rho_val;
      if (rho != NULL) *rho = rho_val;
      if (e != NULL)   *e = h - p / rho_val;
      if (T != NULL)   *T = T_val;
      if (s != NULL)   *s = R3_specific_entropy(rho_val, T_val);
      if (cv != NULL)  *cv = R3_cv(rho_val, T_val);
      if (cp != NULL)  *cp = R3_cp(rho_val, T_val);
      if (c != NULL)   *c = R3_sound_speed(rho_val, T_val);
      if (k != NULL)   *k = thermal_conductivity_R3(rho_val, T_val);
      if (mu != NULL)  *mu = viscosity(rho_val, T_val);
    }
    case 4:
    {
      double h_l_sat = h_l_sat_from_p(p);
      double h_g_sat = h_g_sat_from_p(p);
      double x = (h - h_l_sat) / (h_g_sat - h_l_sat);

      double rho_l_sat = rho_l_sat_from_p(p);
      double rho_g_sat = rho_g_sat_from_p(p);
      double v_l_sat = 1.0 / rho_l_sat;
      double v_g_sat = 1.0 / rho_g_sat;
      double v_val = x * v_g_sat + (1.0 - x) * v_l_sat;

      if (v != NULL)   *v = v_val;
      if (rho != NULL) *rho = 1.0 / v_val;
      if (e != NULL)   *e = h - p * v_val;
      if (T != NULL)   *T = T_sat_from_p(p);
      if (s != NULL)
      {
        double s_l_sat = s_l_sat_from_p(p);
        double s_g_sat = s_g_sat_from_p(p);
        *s = x * s_g_sat + (1.0 - x) * s_l_sat;
      }
      if (cv != NULL)  *cv = -1.0;
      if (cp != NULL)  *cp = -1.0;
      if (c != NULL)   *c = -1.0;
      if (k != NULL)   *k = -1.0;
      if (mu != NULL)  *mu = -1.0;
    }
    case 5:
    {
      double T_val = R5_T_from_p_h_ITER(p, h);
      double v_val = R5_specific_volume(p, T_val);
      double rho_val = 1.0 / v_val;
      if (v != NULL)   *v = v_val;
      if (rho != NULL) *rho = rho_val;
      if (e != NULL)   *e = h - p * v_val;
      if (T != NULL)   *T = T_val;
      if (s != NULL)   *s = R5_specific_entropy(p, T_val);
      if (cv != NULL)  *cv = R5_cv(p, T_val);
      if (cp != NULL)  *cp = R5_cp(p, T_val);
      if (c != NULL)   *c = R5_sound_speed(p, T_val);
      if (k != NULL)   *k = thermal_conductivity_R5(p, T_val);
      if (mu != NULL)  *mu = viscosity(rho_val, T_val);
    }
    default:
      fprintf(stderr, "%s", "Region not recognized!\n");
      exit(1);
  }
}
/***************************************************************
 * (p, s)-based properties
 ***************************************************************/
int locateRegion_from_ps(double p, double s)
{
  // 21 -> 2a; 22 -> 2b; 23 -> 2c.
  if (p < IF97_SAT_P_MIN)
  {
    if (s < R2_specific_entropy(p, IF97_T_MIN))
    {
      fprintf(stderr, "%s", "Out of range: s < s(p, 273.15K)!\n");
      return -1;
    }
    else if (s < R2_specific_entropy(p, IF97_T_25))
      return 21;
    else if (s < R5_specific_entropy(p, IF97_T_MAX))
      return 5;
    else
    {
      fprintf(stderr, "%s", "Out of range: s > s(p, 2273.15K)!\n");
      return -2;
    }
  }
  else if (p < 4.0e6) // boundary between 2a and 2b (see page 21, Ref. [1])
  {
    if (s < R1_specific_entropy(p, IF97_T_MIN))
    {
      fprintf(stderr, "%s", "Out of range: s < s(p, 273.15K)!\n");
      return -1;
    }
    else if (s < s_l_sat_from_p(p))
      return 1;
    else if (s < s_g_sat_from_p(p))
      return 4;
    else if (s < R2_specific_entropy(p, IF97_T_25))
      return 21;
    else if (s < R5_specific_entropy(p, IF97_T_MAX))
      return 5;
    else
    {
      fprintf(stderr, "%s", "Out of range: s > s(p, 2273.15K)!\n");
      return -2;
    }
  }
  else if (p < 6.5467e6) // where 2b-2c boundary line intersects with saturation line, (see page 21, Ref. [1])
  {
    if (s < R1_specific_entropy(p, IF97_T_MIN))
    {
      fprintf(stderr, "%s", "Out of range: s < s(p, 273.15K)!\n");
      return -1;
    }
    else if (s < s_l_sat_from_p(p))
      return 1;
    else if (s < s_g_sat_from_p(p))
      return 4;
    else if (s < R2_specific_entropy(p, IF97_T_25))
      return 22;
    else if (s < R5_specific_entropy(p, IF97_T_MAX))
      return 5;
    else
    {
      fprintf(stderr, "%s", "Out of range: s > s(p, 2273.15K)!\n");
      return -2;
    }
  }
  else if (p < 1.65291643e7) // where 2-3 boundary line intersects with saturation line, (see page 6, Ref. [1])
  {
    if (s < R1_specific_entropy(p, IF97_T_MIN))
    {
      fprintf(stderr, "%s", "Out of range: s < s(p, 273.15K)!\n");
      return -1;
    }
    else if (s < s_l_sat_from_p(p))
      return 1;
    else if (s < s_g_sat_from_p(p))
      return 4;
    else if (s < 5.85e3)
      return 23;
    else if (s < R2_specific_entropy(p, IF97_T_25))
      return 22;
    else if (s < R5_specific_entropy(p, IF97_T_MAX))
      return 5;
    else
    {
      fprintf(stderr, "%s", "Out of range: s > s(p, 2273.15K)!\n");
      return -2;
    }
  }
  else if (p < P_CRIT)
  {
    if (s < R1_specific_entropy(p, IF97_T_MIN))
    {
      fprintf(stderr, "%s", "Out of range: s < s(p, 273.15K)!\n");
      return -1;
    }
    else if (s < R1_specific_entropy(p, IF97_T_13))
      return 1;
    else if (s < s_l_sat_from_p(p))
      return 3;
    else if (s < s_g_sat_from_p(p))
      return 4;
    else if (s < R2_specific_entropy(p, B23_T_from_p(p)))
      return 3;
    else if (s < 5.85e3)
      return 23;
    else if (s < R2_specific_entropy(p, IF97_T_25))
      return 22;
    else if (s < R5_specific_entropy(p, IF97_T_MAX))
      return 5;
    else
    {
      fprintf(stderr, "%s", "Out of range: s > s(p, 2273.15K)!\n");
      return -2;
    }
  }
  else if (p < 50.0e6)
  {
    if (s < R1_specific_entropy(p, IF97_T_MIN))
    {
      fprintf(stderr, "%s", "Out of range: s < s(p, 273.15K)!\n");
      return -1;
    }
    else if (s < R1_specific_entropy(p, IF97_T_13))
      return 1;
    else if (s < R2_specific_entropy(p, B23_T_from_p(p)))
      return 3;
    else if (s < 5.85e3)
      return 23;
    else if (s < R2_specific_entropy(p, IF97_T_25))
      return 22;
    else if (s < R5_specific_entropy(p, IF97_T_MAX))
      return 5;
    else
    {
      fprintf(stderr, "%s", "Out of range: s > s(p, 2273.15K)!\n");
      return -2;
    }
  }
  else if (p <= IF97_P_MAX)
  {
    if (s < R1_specific_entropy(p, IF97_T_MIN))
    {
      fprintf(stderr, "%s", "Out of range: s < s(p, 273.15K)!\n");
      return -1;
    }
    else if (s < R1_specific_entropy(p, IF97_T_13))
      return 1;
    else if (s < R2_specific_entropy(p, B23_T_from_p(p)))
      return 3;
    else if (s < 5.85e3)
      return 23;
    else if (s < R2_specific_entropy(p, IF97_T_25))
      return 22;
    else
    {
      fprintf(stderr, "%s", "Out of range: s > s(p, 1073.15K)!\n");
      return -4;
    }
  }
  else
  {
    return -3;
  }
}

double v_from_ps(double p, double s)
{
  int region = locateRegion_from_ps(p, s);
  switch (region) {
    case 1:
      return R1_specific_volume(p, R1_T_from_p_s(p, s));
    case 21:
      return R2_specific_volume(p, R2a_T_from_p_s(p, s));
    case 22:
      return R2_specific_volume(p, R2b_T_from_p_s(p, s));
    case 23:
      return R2_specific_volume(p, R2c_T_from_p_s(p, s));
    case 3:
    {
      double rho = 0.0, T = 0.0, x = 0.0;
      R3_rho_T_x_from_p_s_ITER(p, s, rho, T, x); /* This function covers Region 3 and the subsetion of 4 within 3 */
      return 1.0 / rho;
    }
    case 4:
    {
      double v_l_sat = 1.0 / rho_l_sat_from_p(p);
      double v_g_sat = 1.0 / rho_g_sat_from_p(p);
      double s_l_sat = s_l_sat_from_p(p);
      double s_g_sat = s_g_sat_from_p(p);
      double x = (s - s_l_sat) / (s_g_sat - s_l_sat); // TODO: correct?
      return x * v_g_sat + (1.0 - x) * v_l_sat;
    }
    case 5:
      return R5_specific_volume(p, R5_T_from_p_s_ITER(p, s));
    default:
      fprintf(stderr, "%s", "Region not recognized!\n");
      exit(1);
      return 0.0;
  }
}

double rho_from_ps(double p, double s)
{
  int region = locateRegion_from_ps(p, s);
  switch (region) {
    case 1:
      return 1.0 / R1_specific_volume(p, R1_T_from_p_s(p, s));
    case 21:
      return 1.0 / R2_specific_volume(p, R2a_T_from_p_s(p, s));
    case 22:
      return 1.0 / R2_specific_volume(p, R2b_T_from_p_s(p, s));
    case 23:
      return 1.0 / R2_specific_volume(p, R2c_T_from_p_s(p, s));
    case 3:
    {
      double rho = 0.0, T = 0.0, x = 0.0;
      R3_rho_T_x_from_p_s_ITER(p, s, rho, T, x); /* This function covers Region 3 and the subsetion of 4 within 3 */
      return rho;
    }
    case 4:
    {
      double rho_l_sat = rho_l_sat_from_p(p);
      double rho_g_sat = rho_g_sat_from_p(p);
      double s_l_sat = s_l_sat_from_p(p);
      double s_g_sat = s_g_sat_from_p(p);
      double x = (s - s_l_sat) / (s_g_sat - s_l_sat);
      return 1.0 / ((1.0 - x) / rho_l_sat + x / rho_g_sat);
    }
    case 5:
      return 1.0 / R5_specific_volume(p, R5_T_from_p_s_ITER(p, s));
    default:
      fprintf(stderr, "%s", "Region not recognized!\n");
      exit(1);
      return 0.0;
  }
}

double e_from_ps(double p, double s)
{
  int region = locateRegion_from_ps(p, s);
  switch (region) {
    case 1:
      return 1.0 / R1_specific_int_energy(p, R1_T_from_p_s(p, s));
    case 21:
      return 1.0 / R2_specific_int_energy(p, R2a_T_from_p_s(p, s));
    case 22:
      return 1.0 / R2_specific_int_energy(p, R2b_T_from_p_s(p, s));
    case 23:
      return 1.0 / R2_specific_int_energy(p, R2c_T_from_p_s(p, s));
    case 3:
    {
      double rho = 0.0, T = 0.0, x = 0.0;
      R3_rho_T_x_from_p_s_ITER(p, s, rho, T, x); /* This function covers Region 3 and the subsetion of 4 within 3 */
      return R3_specific_int_energy(rho, T);
    }
    case 4:
    {
      double e_l_sat = e_l_sat_from_p(p);
      double e_g_sat = e_g_sat_from_p(p);
      double s_l_sat = s_l_sat_from_p(p);
      double s_g_sat = s_g_sat_from_p(p);
      double x = (s - s_l_sat) / (s_g_sat - s_l_sat);
      return e_g_sat * x + e_l_sat * (1.0 - x);
    }
    case 5:
      return 1.0 / R5_specific_int_energy(p, R5_T_from_p_s_ITER(p, s));
    default:
      fprintf(stderr, "%s", "Region not recognized!\n");
      exit(1);
      return 0.0;
  }
}

double T_from_ps(double p, double s)
{
  int region = locateRegion_from_ps(p, s);
  switch (region) {
    case 1:
      return R1_T_from_p_s(p, s);
    case 21:
      return R2a_T_from_p_s(p, s);
    case 22:
      return R2b_T_from_p_s(p, s);
    case 23:
      return R2c_T_from_p_s(p, s);
    case 3:
    {
      double rho = 0.0, T = 0.0, x = 0.0;
      R3_rho_T_x_from_p_s_ITER(p, s, rho, T, x); /* This function covers Region 3 and the subsetion of 4 within 3 */
      return T;
    }
    case 4:
    {
      /*
      double h_l_sat = h_l_sat_from_p(p);
      double h_g_sat = h_g_sat_from_p(p);
      double x = (h - h_l_sat) / (h_g_sat - h_l_sat);*/
      return T_sat_from_p(p);
    }
    case 5:
      return R5_T_from_p_s_ITER(p, s);
    default:
      fprintf(stderr, "%s", "Region not recognized!\n");
      exit(1);
      return 0.0;
  }
}

double h_from_ps(double p, double s)
{
  int region = locateRegion_from_ps(p, s);
  switch (region) {
    case 1:
      return R1_specific_enthalpy(p, R1_T_from_p_s(p, s));
    case 21:
      return R2_specific_enthalpy(p, R2a_T_from_p_s(p, s));
    case 22:
      return R2_specific_enthalpy(p, R2b_T_from_p_s(p, s));
    case 23:
      return R2_specific_enthalpy(p, R2c_T_from_p_s(p, s));
    case 3:
    {
      double rho = 0.0, T = 0.0, x = 0.0;
      R3_rho_T_x_from_p_s_ITER(p, s, rho, T, x); /* This function covers Region 3 and the subsetion of 4 within 3 */
      return R3_specific_enthalpy(rho, T);
    }
    case 4:
    {
      double h_l_sat = h_l_sat_from_p(p);
      double h_g_sat = h_g_sat_from_p(p);
      double s_l_sat = s_l_sat_from_p(p);
      double s_g_sat = s_g_sat_from_p(p);
      double x = (s - s_l_sat) / (s_g_sat - s_l_sat);
      return x * h_g_sat + (1.0 - x) * h_l_sat; // Check correctness
    }
    case 5:
      return R5_specific_enthalpy(p, R5_T_from_p_s_ITER(p, s));
    default:
      fprintf(stderr, "%s", "Region not recognized!\n");
      exit(1);
      return 0.0;
  }
}

double cv_from_ps(double p, double s)
{
  int region = locateRegion_from_ps(p, s);
  switch (region) {
    case 1:
      return R1_cv(p, R1_T_from_p_s(p, s));
    case 21:
      return R2_cv(p, R2a_T_from_p_s(p, s));
    case 22:
      return R2_cv(p, R2b_T_from_p_s(p, s));
    case 23:
      return R2_cv(p, R2c_T_from_p_s(p, s));
    case 3:
    {
      double rho = 0.0, T = 0.0, x = 0.0;
      R3_rho_T_x_from_p_s_ITER(p, s, rho, T, x); /* This function covers Region 3 and the subsetion of 4 within 3 */
      return R3_cv(rho, T);
    }
    case 4:
    {
      return -1.0; // not defined
    }
    case 5:
      return R5_cv(p, R5_T_from_p_s_ITER(p, s));
    default:
      fprintf(stderr, "%s", "Region not recognized!\n");
      exit(1);
      return 0.0;
  }
}

double cp_from_ps(double p, double s)
{
  int region = locateRegion_from_ps(p, s);
  switch (region) {
    case 1:
      return R1_cp(p, R1_T_from_p_s(p, s));
    case 21:
      return R2_cp(p, R2a_T_from_p_s(p, s));
    case 22:
      return R2_cp(p, R2b_T_from_p_s(p, s));
    case 23:
      return R2_cp(p, R2c_T_from_p_s(p, s));
    case 3:
    {
      double rho = 0.0, T = 0.0, x = 0.0;
      R3_rho_T_x_from_p_s_ITER(p, s, rho, T, x); /* This function covers Region 3 and the subsetion of 4 within 3 */
      return R3_cp(rho, T);
    }
    case 4:
    {
      return -1.0; // not defined
    }
    case 5:
      return R5_cp(p, R5_T_from_p_s_ITER(p, s));
    default:
      fprintf(stderr, "%s", "Region not recognized!\n");
      exit(1);
      return 0.0;
  }
}

double c_from_ps(double p, double s)
{
  int region = locateRegion_from_ps(p, s);
  switch (region) {
    case 1:
      return R1_sound_speed(p, R1_T_from_p_s(p, s));
    case 21:
      return R2_sound_speed(p, R2a_T_from_p_s(p, s));
    case 22:
      return R2_sound_speed(p, R2b_T_from_p_s(p, s));
    case 23:
      return R2_sound_speed(p, R2c_T_from_p_s(p, s));
    case 3:
    {
      double rho = 0.0, T = 0.0, x = 0.0;
      R3_rho_T_x_from_p_s_ITER(p, s, rho, T, x); /* This function covers Region 3 and the subsetion of 4 within 3 */
      return R3_sound_speed(rho, T);
    }
    case 4:
    {
      return -1.0; // needs special model, not this package's interest
    }
    case 5:
      return R5_sound_speed(p, R5_T_from_p_s_ITER(p, s));
    default:
      fprintf(stderr, "%s", "Region not recognized!\n");
      exit(1);
      return 0.0;
  }
}

double k_from_ps(double p, double s)
{
  int region = locateRegion_from_ps(p, s);
  switch (region) {
    case 1:
      return thermal_conductivity_R1(p, R1_T_from_p_s(p, s));
    case 21:
      return thermal_conductivity_R2(p, R2a_T_from_p_s(p, s));
    case 22:
      return thermal_conductivity_R2(p, R2b_T_from_p_s(p, s));
    case 23:
      return thermal_conductivity_R2(p, R2c_T_from_p_s(p, s));
    case 3:
    {
      double rho = 0.0, T = 0.0, x = 0.0;
      R3_rho_T_x_from_p_s_ITER(p, s, rho, T, x); /* This function covers Region 3 and the subsetion of 4 within 3 */
      return thermal_conductivity_R3(rho, T);
    }
    case 4:
    {
      return -1.0; // not defined
    }
    case 5:
      return thermal_conductivity_R5(p, R5_T_from_p_s_ITER(p, s));
    default:
      fprintf(stderr, "%s", "Region not recognized!\n");
      exit(1);
      return 0.0;
  }
}
double mu_from_ps(double p, double s)
{
  int region = locateRegion_from_ps(p, s);
  switch (region) {
    case 1:
    case 21:
    case 22:
    case 23:
    case 3:
    case 5:
    {
      double rho = rho_from_ps(p, s);
      double T = T_from_ps(p, s);
      return viscosity(rho, T);
    }
    case 4:
    {
      return -1.0; // not defined
    }
    default:
      fprintf(stderr, "%s", "Region not recognized!\n");
      exit(1);
      return 0.0;
  }
}

/***************************************************************
 * (h, v)-based properties
 ***************************************************************/
double p_from_hv(double h, double v)
{
  double p_min = IF97_SAT_P_MIN;
  double p_max = p_max_from_h(h);

  double p_find, v_find, p_error = 1.0;
  unsigned int it = 0;
  while ((p_error > 1.0e-9) && it < 1000)
  {
    p_find = 0.5 * (p_min + p_max);
    v_find = v_from_ph(p_find, h);

    if (v_find > v)   p_min = p_find;
    else              p_max = p_find;

    p_error = fabs((p_max - p_min) / p_find);
    it ++;
    /*if (it > 900)
      printf("it > 900\n");*/
  }

  return p_find;
}
