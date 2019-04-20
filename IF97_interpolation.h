#ifndef IF97_INTERPOLATION_H
#define IF97_INTERPOLATION_H

#include <vector>

#include "IF97.h"
#include "IF97_helper.h"

struct IF97_Interpolation
{
  IF97_Interpolation();
  ~IF97_Interpolation();

  void loadSatLineData();
  double INTPL_rho_l_sat_from_T(double T);
  double INTPL_rho_g_sat_from_T(double T);

private:
  std::vector<double> T_sat;
  std::vector<double> rho_l_sat;
  std::vector<double> rho_g_sat;
};

#endif /*IF97_INTERPOLATION_H*/
