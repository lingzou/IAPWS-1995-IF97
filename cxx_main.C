#include <iostream>
#include <iomanip>
#include <cstring>

#include "IF97.h"
#include "IF97_helper.h"
//#include "IF97_interpolation.h"
#include "unit_tests.h"
#include "SurfaceTension.h"
#include "Viscosity.h"
#include "ThermalConductivity.h"
#include "IF97_interface.h"

int main(int argc, char *argv[])
{
  std::cout << "Hello, IF97" << std::endl;

  if ((argc == 2) && (strcmp(argv[1], "-test_all") == 0))
  {
    Unit_Test_All();
    return 0;
  }
  else if ((argc == 3) && (strcmp(argv[1], "-unit_test") == 0))
  {
    std::string str(argv[2]);
    Unit_Test(str);
    return 0;
  }

  std::cout << std::scientific << std::setprecision(12) << "psat = " << p_sat_from_T(623.15) << std::endl;
  /*
  std::cout << std::scientific << std::setprecision(8);
  for (int i = 0; i < 100; i++)
  {
    double r = ((double) rand()) / (double) RAND_MAX;
    double T = 623.15 + (T_CRIT - 623.15) * r;

    std::cout << T << "; " << rho_l_sat_from_T(T) << std::endl;
  }

  std::cout << "p = 6.54670 MPa, T_sat = " << T_sat_from_p(6.54670e6)
            << std::endl
            << "h_sat = " << h_g_sat_from_p(6.54670e6) << std::endl
            << "h_2bc = " << B2bc_h_from_p(6.54670e6)  << std::endl;*/

  std::cout << "h_l_sat(T=T_crit) = " << h_l_sat_from_T(T_CRIT) << std::endl;
  std::cout << "h_g_sat(T=T_crit) = " << h_g_sat_from_T(T_CRIT) << std::endl;
  std::cout << "T_sat_from_p(p=p_crit) = " << T_sat_from_p(P_CRIT) << std::endl;
  std::cout << "h_l_sat(p=p_crit) = " << h_l_sat_from_p(P_CRIT) << std::endl;
  std::cout << "h_g_sat(p=p_crit) = " << h_g_sat_from_p(P_CRIT) << std::endl;

  std::cout << "R3_crit = " << R3_specific_enthalpy(RHO_CRIT, T_CRIT) << std::endl;

  std::cout << "h(p=1.0bar, T=273.15) = " << h_from_pT(1.0e5, 273.15) << std::endl;
  //std::cout << "T(p=1.0bar, h=0.0) = " << T_from_ph(1.0e5, 0.0) << std::endl;
  std::cout << "s(p=1.0bar, T=273.15) = " << s_from_pT(1.0e5, 273.15) << std::endl;

  std::cout << "h(p=50MPa, T=2253.15) = " << h_from_pT(50.e6, 2253.15) << std::endl;
  std::cout << "s(p=50MPa, T=2253.15) = " << s_from_pT(50.e6, 2253.15) << std::endl;

  std::cout << "T_sat(p=611.213Pa) = " << T_sat_from_p(611.213) << std::endl;
  std::cout << "T_sat(p=611.213Pa) = " << T_sat_from_p(6.112126774443e2) << std::endl;
  std::cout << "p_sat(p=273.15K) = " << p_sat_from_T(273.15) << std::endl;

  std::cout << "1: " << IF97_H_PMAX_TMIN << std::endl;
  std::cout << "2: " << IF97_H_PMAX_T25 << std::endl;
  std::cout << "3: " << IF97_H_PMID_T25 << std::endl;
  std::cout << "4: " << IF97_H_PMID_TMAX << std::endl;
/*
  genR3_sat_line();

  IF97_Interpolation IF97InterP;
  //IF97InterP.init();
  std::cout << std::scientific << std::setprecision(8);
  double T = 250.0;
  //std::cout << "T = " << T << IF97InterP.INTPL_rho_l_sat_from_T(T) << std::endl;
  //std::cout << "T = " << T << IF97InterP.INTPL_rho_g_sat_from_T(T) << std::endl;
  std::cout << IF97InterP.INTPL_property_from_T(T, IF97_Interpolation::RHO_L) << std::endl;
  std::cout << IF97InterP.INTPL_property_from_T(T, IF97_Interpolation::RHO_G) << std::endl;

  T = 700.0;
  //std::cout << "T = " << T << IF97InterP.INTPL_rho_l_sat_from_T(T) << std::endl;
  //std::cout << "T = " << T << IF97InterP.INTPL_rho_g_sat_from_T(T) << std::endl;
  std::cout << IF97InterP.INTPL_property_from_T(T, IF97_Interpolation::RHO_L) << std::endl;
  std::cout << IF97InterP.INTPL_property_from_T(T, IF97_Interpolation::RHO_G) << std::endl;

  T = 273.5;
  //std::cout << "T = " << T << IF97InterP.INTPL_rho_l_sat_from_T(T) << std::endl;
  //std::cout << "T = " << T << IF97InterP.INTPL_rho_g_sat_from_T(T) << std::endl;
  std::cout << IF97InterP.INTPL_property_from_T(T, IF97_Interpolation::RHO_L) << std::endl;
  std::cout << IF97InterP.INTPL_property_from_T(T, IF97_Interpolation::RHO_G) << std::endl;

  T = 300.0;
  //std::cout << "T = " << T << IF97InterP.INTPL_rho_l_sat_from_T(T) << std::endl;
  //std::cout << "T = " << T << IF97InterP.INTPL_rho_g_sat_from_T(T) << std::endl;
  std::cout << IF97InterP.INTPL_property_from_T(T, IF97_Interpolation::RHO_L) << std::endl;
  std::cout << IF97InterP.INTPL_property_from_T(T, IF97_Interpolation::RHO_G) << std::endl;

  T = 410.0;
  //std::cout << "T = " << T << IF97InterP.INTPL_rho_l_sat_from_T(T) << std::endl;
  //std::cout << "T = " << T << IF97InterP.INTPL_rho_g_sat_from_T(T) << std::endl;
  std::cout << IF97InterP.INTPL_property_from_T(T, IF97_Interpolation::RHO_L) << std::endl;
  std::cout << IF97InterP.INTPL_property_from_T(T, IF97_Interpolation::RHO_G) << std::endl;

  T = 555.0;
  //std::cout << "T = " << T << IF97InterP.INTPL_rho_l_sat_from_T(T) << std::endl;
  //std::cout << "T = " << T << IF97InterP.INTPL_rho_g_sat_from_T(T) << std::endl;
  std::cout << IF97InterP.INTPL_property_from_T(T, IF97_Interpolation::RHO_L) << std::endl;
  std::cout << IF97InterP.INTPL_property_from_T(T, IF97_Interpolation::RHO_G) << std::endl;

  T = 623.0;
  //std::cout << "T = " << T << IF97InterP.INTPL_rho_l_sat_from_T(T) << std::endl;
  //std::cout << "T = " << T << IF97InterP.INTPL_rho_g_sat_from_T(T) << std::endl;
  std::cout << IF97InterP.INTPL_property_from_T(T, IF97_Interpolation::RHO_L) << std::endl;
  std::cout << IF97InterP.INTPL_property_from_T(T, IF97_Interpolation::RHO_G) << std::endl;

  T = 6.47095550e+02;
  //std::cout << "T = " << T << IF97InterP.INTPL_rho_l_sat_from_T(T) << std::endl;
  //std::cout << "T = " << T << IF97InterP.INTPL_rho_g_sat_from_T(T) << std::endl;
  std::cout << IF97InterP.INTPL_property_from_T(T, IF97_Interpolation::RHO_L) << std::endl;
  std::cout << IF97InterP.INTPL_property_from_T(T, IF97_Interpolation::RHO_G) << std::endl;


  double T_array[4] = {298.15, 298.15, 298.15, 873.15};
  double rho_array[4] = {0.0, 998.0, 1200.0, 0.0};
  for (int i = 0; i < 4; i++)
    std::cout << thermal_conductivity_no_enhancement(rho_array[i], T_array[i]) << std::endl;
*/
  return 0;
}
