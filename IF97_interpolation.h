#ifndef IF97_INTERPOLATION_H
#define IF97_INTERPOLATION_H

#include <vector>
#include <map>

#include "IF97.h"
#include "IF97_helper.h"

struct IF97_Interpolation
{
  IF97_Interpolation();
  ~IF97_Interpolation();

  enum INTPL_PROPERTY {
    P = 0,
    RHO_L = 1,  V_L = 2,  E_L = 3,  S_L = 4,  H_L = 5,  CV_L = 6,  CP_L = 7,  C_L = 8,
    RHO_G = 9,  V_G = 10, E_G = 11, S_G = 12, H_G = 13, CV_G = 14, CP_G = 15, C_G = 16
  };

  void loadSatLineData();
  double INTPL_rho_l_sat_from_T(double T);
  double INTPL_rho_g_sat_from_T(double T);
  double INTPL_property_from_T(double T, INTPL_PROPERTY property);

private:
  std::vector<double> T_sat;
  std::vector<double> p_sat;

  std::vector<double> rho_l_sat;
  std::vector<double> v_l_sat;
  std::vector<double> e_l_sat;
  std::vector<double> s_l_sat;
  std::vector<double> h_l_sat;
  std::vector<double> cv_l_sat;
  std::vector<double> cp_l_sat;
  std::vector<double> c_l_sat;

  std::vector<double> rho_g_sat;
  std::vector<double> v_g_sat;
  std::vector<double> e_g_sat;
  std::vector<double> s_g_sat;
  std::vector<double> h_g_sat;
  std::vector<double> cv_g_sat;
  std::vector<double> cp_g_sat;
  std::vector<double> c_g_sat;

  std::map<INTPL_PROPERTY, std::vector<double> *> enum_to_vector_map;
};

#endif /*IF97_INTERPOLATION_H*/
