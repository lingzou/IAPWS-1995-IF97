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
#include "IAPWS1995_Rev.h"
#include "IAPWS1995_Interface.h"

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

  //std::cout << std::scientific << std::setprecision(12) << "psat = " << IAPWS1995Rev::psat_from_T(623.15) << std::endl;
  //std::cout << std::scientific << std::setprecision(12) << "Tsat = " << IAPWS1995Rev::Tsat_from_p(1.652933992355e7) << std::endl;
  //std::cout << std::scientific << std::setprecision(12) << "rho_l(1.652933992355e7, 620.15) = " << IAPWS1995Rev::rho_l_from_pT(1.652933992355e7, 620.15) << std::endl;
  //std::cout << std::scientific << std::setprecision(12) << "rho_l(1.652933992355e7, 625.15) = " << IAPWS1995Rev::rho_l_from_pT(1.652933992355e7, 625.15) << std::endl;


  std::cout << "p = " << IAPWS1995Rev::Pressure(620.00, 3.46015969e+02) << std::endl;
  std::cout << "p = " << IAPWS1995Rev::Pressure(620.00, 6.77938862e+02) << std::endl;
  std::cout << "p = " << IAPWS1995Rev::Pressure(636.00, 154.8611) << std::endl;
  std::cout << "p = " << IAPWS1995Rev::Pressure(636.00, 150.0885) << std::endl;  

  std::cout << "Tsat(19.3MPa) = " << IAPWS1995Rev::Tsat_from_p(19.3e6) << std::endl;

  double rho_l = 0; bool search_failed = true;
  /*
  std::cout << "T_sat(1e7) = " << IAPWS1995Rev::Tsat_from_p(1e7) << std::endl;
  IAPWS1995Rev::rho_l_from_pT(1e7, 529, rho_l, search_failed);
  std::cout << std::scientific << std::setprecision(12) << "rho_l(1e7, 529) = " << rho_l << std::endl;
  IAPWS1995Rev::rho_l_from_pT(1e7, 553, rho_l, search_failed);
  std::cout << std::scientific << std::setprecision(12) << "rho_l(1e7, 553) = " << rho_l << std::endl;
  IAPWS1995Rev::rho_l_from_pT(1e7, 585, rho_l, search_failed);
  std::cout << std::scientific << std::setprecision(12) << "rho_l(1e7, 585) = " << rho_l << std::endl;
  IAPWS1995Rev::rho_l_from_pT(1e7, 598, rho_l, search_failed);
  std::cout << std::scientific << std::setprecision(12) << "rho_l(1e7, 598) = " << rho_l << std::endl;
*/
  IAPWS1995Rev::rho_l_from_pT(16.7e6, 624, rho_l, search_failed);
  std::cout << std::scientific << std::setprecision(12) << "rho_l(16.7e6, 624) = " << (search_failed ? 0 : rho_l) << std::endl;
  std::cout << std::scientific << std::setprecision(12) << "e_l(16.7e6, 624) = " << IAPWS1995Rev::e_l_from_pT(16.7e6, 624) << std::endl;

  std::cout << std::scientific << std::setprecision(12) << "rho_g_meta(4.7e6, 455) = " << 1.0 / R2Meta_specific_volume(4.7e6, 455) << std::endl;
/*
  std::cout << std::scientific << std::setprecision(12) << "rho_g(1e7, 561) = " << IAPWS1995Rev::rho_g_from_pT(1e7, 561) << std::endl;
  std::cout << std::scientific << std::setprecision(12) << "rho_g(1e7, 584) = " << IAPWS1995Rev::rho_g_from_pT(1e7, 584) << std::endl;
  std::cout << std::scientific << std::setprecision(12) << "rho_g(1e7, 606) = " << IAPWS1995Rev::rho_g_from_pT(1e7, 606) << std::endl;
*/

  //std::cout << "p = " << IAPWS1995Rev::Pressure(585, 3.438828230716e2) << std::endl;
/*
  for (int i = 0; i < 400; i++)
  {
    double rho = 7.563471314695e2 - i*2;
    double T = 585;
    std::cout << i << ": " << IAPWS1995Rev::Pressure(T, rho)
              << "; rho = " << rho << std::endl;
  }*/
  /*
  double T_min = 273.16;
  double T_max = 647.096;
  int N = 100;
  double dT = (T_max - T_min) / N;
  for (int i = 0; i < N + 1; i++)
  {
    double T = T_min + dT * i;
    std::cout << p_sat_from_T(T) << ", "
              << rho_l_sat_from_T(T) << ", "
              << rho_g_sat_from_T(T) << ", "
              << e_l_sat_from_T(T) << ", "
              << e_g_sat_from_T(T) << ", "
              << h_l_sat_from_T(T) << ", "
              << h_g_sat_from_T(T) << ", "
              << s_l_sat_from_T(T) << ", "
              << s_g_sat_from_T(T) << ", "
              << IAPWS1995Rev::psat_from_T(T) << ", "
              << IAPWS1995Rev::rho_l_sat_from_T(T) << ", "
              << IAPWS1995Rev::rho_g_sat_from_T(T) << ", "
              << IAPWS1995Rev::e_l_sat_from_T(T) << ", "
              << IAPWS1995Rev::e_g_sat_from_T(T) << ", "
              << IAPWS1995Rev::h_l_sat_from_T(T) << ", "
              << IAPWS1995Rev::h_g_sat_from_T(T) << ", "
              << IAPWS1995Rev::s_l_sat_from_T(T) << ", "
              << IAPWS1995Rev::s_g_sat_from_T(T) << std::endl;
  }*/
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

  std::cout << "h(p=1Pa, T=2273.15) = " << h_from_pT(1.0, 2273.15) << std::endl;
  std::cout << "h(p=50MPa, T=2273.15) = " << h_from_pT(50.e6, 2273.15) << std::endl;
  std::cout << "s(p=50MPa, T=2273.15) = " << s_from_pT(50.e6, 2273.15) << std::endl;

  std::cout << "T_sat(p=611.213Pa) = " << T_sat_from_p(611.213) << std::endl;
  std::cout << "T_sat(p=611.213Pa) = " << T_sat_from_p(6.112126774443e2) << std::endl;
  std::cout << "p_sat(p=273.15K) = " << p_sat_from_T(273.15) << std::endl;

  std::cout << "1: " << IF97_H_PMAX_TMIN << std::endl;
  std::cout << "2: " << IF97_H_PMAX_T25 << std::endl;
  std::cout << "3: " << IF97_H_PMID_T25 << std::endl;
  std::cout << "4: " << IF97_H_PMID_TMAX << std::endl;

  std::cout << "p_max_from_h(50e3) = " << p_max_from_h(50e3) << std::endl;
  std::cout << "p_max_from_h(200e3) = " << p_max_from_h(200e3) << std::endl;
  std::cout << "p_max_from_h(800e3) = " << p_max_from_h(800e3) << std::endl;
  std::cout << "p_max_from_h(3600e3) = " << p_max_from_h(3600e3) << std::endl;
  std::cout << "p_max_from_h(3800e3) = " << p_max_from_h(3800e3) << std::endl;
  std::cout << "p_max_from_h(4000e3) = " << p_max_from_h(4000e3) << std::endl;

  //p_from_hv(0.0, 1.0e-3);

  std::cout << "IF97_V_PMID_TMAX = " << IF97_V_PMID_TMAX << std::endl;
  std::cout << "IF97_H_PMID_TMAX = " << IF97_H_PMID_TMAX << std::endl;

  std::cout << "h_max_from_v(1.0e-3) = " << h_max_from_v(1.0e-3) << std::endl;
  std::cout << "h_max_from_v(1.5e-3) = " << h_max_from_v(1.5e-3) << std::endl;
  std::cout << "h_max_from_v(5.0e-3) = " << h_max_from_v(5.0e-3) << std::endl;
  std::cout << "h_max_from_v(1.0e-2) = " << h_max_from_v(1.0e-2) << std::endl;

  std::cout << "IF97_VMIN_GLOBAL = " << IF97_VMIN_GLOBAL << std::endl;
  std::cout << "IF97_HMIN_GLOBAL = " << IF97_HMIN_GLOBAL << std::endl;
  std::cout << "IF97_V_SAT_LIQ_TMIN = " << IF97_V_SAT_LIQ_TMIN << std::endl;
  std::cout << "IF97_V_SAT_VAP_TMIN = " << IF97_V_SAT_VAP_TMIN << std::endl;

  std::cout << "IF97_H_SAT_LIQ_TMIN = " << IF97_H_SAT_LIQ_TMIN << std::endl;
  std::cout << "IF97_H_SAT_VAP_TMIN = " << IF97_H_SAT_VAP_TMIN << std::endl;

  std::cout << "+++++++++++++ 9.8e-4\n";
  double val = h_min_from_v(9.8e-4);
  std::cout << val << std::endl;
  std::cout << "+++++++++++++ 9.9e-4\n";
  val = h_min_from_v(9.9e-4);
  std::cout << val << std::endl;
  std::cout << "+++++++++++++ 1.0e-3\n";
  val = h_min_from_v(1.0e-3);
  std::cout << val << std::endl;
  std::cout << "+++++++++++++ 1.02e-3\n";
  val = h_min_from_v(1.02e-3);
  std::cout << val << std::endl;
/*
  double p = 20e6;
  double T_sat = T_sat_from_p(p);
  double rho_l_sat = rho_l_sat_from_p(p);
  double rho_g_sat = rho_g_sat_from_p(p);
  double drho = (rho_l_sat - rho_g_sat) / 100;

  double TT = T_sat + 2.0;
  std::cout << "TT = " << TT << std::endl;
  for (unsigned i = 0; i < 100; i++)
  {
    //double rho = rho_l_sat - 1.0;
    double rho = rho_l_sat - i * drho;
    //double TT = T_sat + 0.002 * i;
    double pp = R3_p(rho, TT);
    std::cout << i << "; " << pp << "; " << rho << std::endl;
    //std::cout << i << "; " << pp << std::endl;
  }
  double rho_g_superheated = R3_rho_from_p_T_ITER(p, TT);
  std::cout << "rho_g_superheated = " << rho_g_superheated << std::endl;
*/
  std::cout << "p = " << 20e6
            << "T = " << 274
            << "rho = " << 1.0 / R1_specific_volume(20e6, 274)
            << "rho = " << rho_from_pT(20e6, 274) << std::endl;
/*
  std::cout << "p_from_hv(200e3, 9.8e-4) = " << p_from_hv(200e3, 9.8e-4) << std::endl;
  std::cout << "p_from_hv(200e3, 0.1) = " << p_from_hv(200e3, 0.1) << std::endl;
  std::cout << "p_from_hv(400e3, 1.02e-3) = " << p_from_hv(400e3, 1.02e-3) << std::endl;
  std::cout << "p_from_hv(1000e3, 0.002) = " << p_from_hv(1000e3, 0.002) << std::endl;*/
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
