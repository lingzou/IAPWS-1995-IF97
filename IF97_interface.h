#ifndef IF97_INTERFACE_H
#define IF97_INTERFACE_H

extern "C"
{
/***************************************************************
 * Saturation line
 ***************************************************************/
double T_sat_from_p(double p);
double p_sat_from_T(double T);

double v_l_sat_from_T(double T);
double rho_l_sat_from_T(double T);
double e_l_sat_from_T(double T);
double h_l_sat_from_T(double T);
double s_l_sat_from_T(double T);
double cv_l_sat_from_T(double T);
double cp_l_sat_from_T(double T);
double c_l_sat_from_T(double T);
double k_l_sat_from_T(double T);
double mu_l_sat_from_T(double T);
void liquid_sat_properties_from_T(double T, double & p, double & v, double & rho, double & e, double & h,
                                  double & s, double & cv, double & cp, double & c, double & k, double & mu);

double v_g_sat_from_T(double T);
double rho_g_sat_from_T(double T);
double e_g_sat_from_T(double T);
double h_g_sat_from_T(double T);
double s_g_sat_from_T(double T);
double cv_g_sat_from_T(double T);
double cp_g_sat_from_T(double T);
double c_g_sat_from_T(double T);
double k_g_sat_from_T(double T);
double mu_g_sat_from_T(double T);
void vapor_sat_properties_from_T(double T, double & p, double & v, double & rho, double & e, double & h,
                                  double & s, double & cv, double & cp, double & c, double & k, double & mu);

}
#endif /*IF97_INTERFACE_H*/
