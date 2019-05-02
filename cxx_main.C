#include <iostream>
#include <iomanip>
#include <cstring>

#include "IF97.h"
#include "IF97_helper.h"
#include "IF97_interpolation.h"
#include "unit_tests.h"

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

  //genR4_sat_line();
/*
  std::cout << R2_cp(1.22818387e3, 283.15) << std::endl;
  std::cout << R2_cv(1.22818387e3, 283.15) << std::endl;
  std::cout << std::scientific << std::setprecision(8) << B23_p_from_T(863.15) << std::endl;
  std::cout << std::scientific << std::setprecision(8) << B23_T_from_p(100.e6) << std::endl;
*/


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

/*
  double p = 50.0e6;
  T = 650.0;
  double p23 = B23_p_from_T(T);
  double rho_min = 1.0 / R2_specific_volume(p23, T);
  double rho_max = 1.0 / R1_specific_volume(100.0e6, T13);
  double dRho = (rho_max - rho_min) / 100.0;
  double tau = Tcrit / T;
  for (int i = 0; i < 101; i++)
  {
    double rho = rho_min + dRho * i;
    double delta = rho / Rhocrit;

    std::cout << rho << "; " << R3_p(rho, T) << std::endl;
  }*/

  //double p = 50.0e6;
  //double T = 650.0;
  //std::cout << "p = " << p << "; T = " << T << "; rho = " << R3_rho_from_p_T_ITER(p, T) << std::endl;

  double T_array[8] = {630.0, 637.0, 647.0, 626.16, 629.0, 697.0, 636.0, 854.0};
  double p_array[8] = {33.5e6, 20.e6, 21.0e6, 17.9e6, 17.37e6, 31.0e6, 97.0e6, 97.7e6};
  for (int i = 0; i < 8; i++)
  {
    double rho = R3_rho_from_p_T_ITER(p_array[i], T_array[i]);
    std::cout
      << "    p = " << std::scientific << std::setprecision(8) << p_array[i]
      << ";   T = " << std::scientific << std::setprecision(8) << T_array[i]
      << "; rho = " << std::scientific << std::setprecision(8) << rho
      << ";   h = " << std::scientific << std::setprecision(8) << R3_specific_enthalpy(rho, T_array[i])
      << std::endl;
  }

  std::cout << std::endl;
  double TT, xx;
  R3_T_x_from_p_h_ITER(3.35e7, 1.64344943e6, TT, xx);
  std::cout << TT << std::endl;
  R3_T_x_from_p_h_ITER(2.0e7, 1.79130276e6, TT, xx);
  std::cout << TT << std::endl;
  R3_T_x_from_p_h_ITER(2.1e7, 2.48929822e+06, TT, xx);
  std::cout << TT << std::endl;
  R3_T_x_from_p_h_ITER(1.79e7, 1.68821864e+06, TT, xx);
  std::cout << TT << std::endl;
  R3_T_x_from_p_h_ITER(1.737e7, 2.56582491e+06, TT, xx);
  std::cout << TT << std::endl;
  R3_T_x_from_p_h_ITER(3.1e7, 2.55197817e+06, TT, xx);
  std::cout << TT << std::endl;
  R3_T_x_from_p_h_ITER(9.7e7, 1.61414193e+06, TT, xx);
  std::cout << TT << std::endl;
  R3_T_x_from_p_h_ITER(9.77e7, 2.77581149e+06, TT, xx);
  std::cout << TT << std::endl;
  std::cout << R3_p(1.19044937e2, 629.0) << std::endl;
  std::cout << R3_p(1.0/8.400177202e-3, 629.0) << std::endl;



/*
  double p = Pcrit;
  double rho = Rhocrit;
  double T = Tcrit;

  std::cout << "P_crit_cal = " << R3_p(rho, T) << std::endl;
*/
/*
  double T13 = 623.15;
  double p_sat13 = p_sat_from_T(T13);
  //std::cout << "T13 = " << T13 << ";" << "p_sat = " << p_sat_from_T(T13) << std::endl;

  double rho_1 = 1.0 / R1_specific_volume(p_sat13, T13);
  double rho_2 = 1.0 / R2_specific_volume(p_sat13, T13);

  std::cout << rho_1 << std::endl;
  std::cout << rho_2 << std::endl;

  std::cout << "p/(RT rho)_1 = " << std::scientific << std::setprecision(8)
    << p_sat13 / (Rgas * T13 * rho_1) << std::endl;
  std::cout << "delta * phi_delta = " << std::scientific << std::setprecision(8)
    << rho_1 / Rhocrit * R3_phi_delta(rho_1 / Rhocrit, Tcrit / T13) << std::endl;

  std::cout << "p/(RT rho)_2 = " << std::scientific << std::setprecision(8)
    << p_sat13 / (Rgas * T13 * rho_2) << std::endl;
  std::cout << "delta * phi_delta = " << std::scientific << std::setprecision(8)
    << rho_2 / Rhocrit * R3_phi_delta(rho_2 / Rhocrit, Tcrit / T13) << std::endl;
*/

/*
  double T13 = 623.15;
  //std::cout << "T13 = " << T13 << ";" << "p_sat = " << p_sat_from_T(T13) << std::endl;

  double p = 22.063e6;
  double Ts = T_sat_from_p(p);
  std::cout << "Ts = " << Ts << std::endl;
  double T23 = B23_T_from_p(p);
  std::cout << "T23 = " << T23 << std::endl;
  double rho_13 = 1.0 / R1_specific_volume(p, T13);
  double rho_23 = 1.0 / R2_specific_volume(p, T23);
  std::cout << "rho_13 = " << rho_13 << std::endl;
  std::cout << "rho_23 = " << rho_23 << std::endl;

  rho_13 = 322.0 + 10.0;
  rho_23 = 322.0 - 10.0;

  for (int i = 0; i < 101; i++)
  {
    double rho = rho_13 + (rho_23 - rho_13) / 100 * i;
    double delta = rho / Rhocrit;
    double tau = Tcrit / Ts;

    double val = p / Rgas / Ts / rho - delta * R3_phi_delta(delta, tau);

    //std::cout << "i = " << i << "; rho = " << rho << " val = " << val << std::endl;
    std::cout << rho << ";" << val << std::endl;
  }*/
/*
  //double rho = 100.0;
  double T_min = 273.15;
  double T_max = 623.15;
  double p_max = 100.0e6;
  int N = 100;

  for (double T = T_min; T <= T_max+0.001; T+=5.0)
  {
    //std::cout << T << ";" << R3_p(rho, T) << std::endl;)
    double ps = p_sat_from_T(T);
    double dp = (p_max - ps) / N;
    for (int i = 0; i < N+1; i++)
    {
      double p = ps + dp * i;
      std::cout << std::scientific << std::setprecision(8) << std::setw(20) << T
      << std::scientific << std::setprecision(8) << std::setw(20) << p
      << std::scientific << std::setprecision(8) << std::setw(20) << 1.0 / R1_specific_volume(p, T) << std::endl;
    }
  }*/

  return 0;
}
