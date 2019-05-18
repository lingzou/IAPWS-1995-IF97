#include <iostream>
#include <iomanip>
#include <cstring>

#include "IF97.h"
#include "IF97_helper.h"
#include "IF97_interpolation.h"
#include "unit_tests.h"
#include "SurfaceTension.h"
#include "Viscosity.h"

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
/*
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
      << ";   s = " << std::scientific << std::setprecision(8) << R3_specific_entropy(rho, T_array[i])
      << std::endl;
  }
*/
/*
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
*/
/*
  std::cout << std::endl;
  double TT, xx;
  R3_T_x_from_p_s_ITER(3.35e7, 3.69009147e+03, TT, xx);
  std::cout << TT << std::endl;
  R3_T_x_from_p_s_ITER(2.0e7, 3.95927258e+03, TT, xx);
  std::cout << TT << std::endl;
  R3_T_x_from_p_s_ITER(2.1e7, 5.04168152e+03, TT, xx);
  std::cout << TT << std::endl;
  R3_T_x_from_p_s_ITER(1.79e7, 3.80223220e+03, TT, xx);
  std::cout << TT << std::endl;
  R3_T_x_from_p_s_ITER(1.737e7, 5.20290950e+03, TT, xx);
  std::cout << TT << std::endl;
  R3_T_x_from_p_s_ITER(3.1e7, 5.05425538e+03, TT, xx);
  std::cout << TT << std::endl;
  R3_T_x_from_p_s_ITER(9.7e7, 3.49976082e+03, TT, xx);
  std::cout << TT << std::endl;
  R3_T_x_from_p_s_ITER(9.77e7, 5.06161454e+03, TT, xx);
  std::cout << TT << std::endl;
*/


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
/*
  double T_array[8] = {1100.0, 1200.0, 1400.0, 1600.0, 1800.0, 2000.0, 2200.0, 2250.0};
  double p_array[8] = {1.0e3, 1.0e6, 5.0e6, 10.0e6, 20.0e6, 30.0e6, 40.0e6, 50.0e6};
  //double h[64], s[64];
  for (int i = 0; i < 8; i++)
  {
    //for (int j = 0; j < 8; j++)
    //{
      //h[i*8+j] = R5_specific_enthalpy(p_array[i], T_array[j]);
      //s[i*8+j] = R5_specific_entropy(p_array[i], T_array[j]);
      std::cout
        << "    p = " << std::scientific << std::setprecision(8) << p_array[i]
        << ";   T = " << std::scientific << std::setprecision(8) << T_array[i]
        << "; rho = " << std::scientific << std::setprecision(8) << 1.0 / R5_specific_volume(p_array[i], T_array[i])
        << ";   h = " << std::scientific << std::setprecision(8) << R5_specific_enthalpy(p_array[i], T_array[i])
        << ";   s = " << std::scientific << std::setprecision(8) << R5_specific_entropy(p_array[i], T_array[i])
        << std::endl;
    //}
  }*/
/*
  for (int i = 0; i < 8; i++)
  {
    for (int j = 0; j < 8; j++)
    {
      std::cout << std::scientific << std::setprecision(8) << R5_T_from_p_h_ITER(p_array[i], h[i*8+j]) << std::endl;
      std::cout << std::scientific << std::setprecision(8) << R5_T_from_p_s_ITER(p_array[i], s[i*8+j]) << std::endl;
    }
  }*/
/*
  std::cout << std::scientific << std::setprecision(8) << R5_T_from_p_h_ITER(1.0e3, 7.30921691e6) << std::endl;
  std::cout << std::scientific << std::setprecision(8) << R5_T_from_p_s_ITER(1.0e3, 1.36477789e4) << std::endl;

  std::cout << std::scientific << std::setprecision(8) << R5_T_from_p_h_ITER(2.0e7, 7.15701330e6) << std::endl;
  std::cout << std::scientific << std::setprecision(8) << R5_T_from_p_s_ITER(2.0e7, 9.00461822e3) << std::endl;
  */
/*
  double T_array[9] = {0.01, 75.0, 100.0, 145.0, 195.0, 225.0, 300.0, 350.0, 370.0}; // in [C]
  for (int i = 0; i < 9; i++)
  {
    std::cout << T_array[i] << "; ";
    std::cout << std::scientific << std::setprecision(3) << SURF_TENSION::surf_tension(T_array[i] + 273.15) << std::endl;
  }
*/
/*
  double T_array[11] = {298.15, 298.15, 373.15, 433.15, 433.15, 873.15, 873.15, 873.15, 1173.15, 1173.15, 1173.15};
  double rho_array[11] = {998, 1200, 1000, 1, 1000, 1, 100, 600, 1, 100, 400};
  for (int i = 0; i < 11; i++)
    std::cout << std::scientific << std::setprecision(9) << viscosity(rho_array[i], T_array[i]) << std::endl;

  std::cout << std::endl << std::endl;
  std::cout << viscosity(122.0, 647.35, true) << std::endl;
*/

  double rho_array[6] = {122, 222, 272, 322, 372, 422};
  double xi_array[6] = {0.309247, 1.571405, 5.266522, 16.590209, 5.603768, 1.876244};
  for (int i = 0; i < 6; i++)
    std::cout << viscosity(rho_array[i], 647.35, true) << std::endl << std::endl;
    //std::cout << mu2_bar(rho_array[i]/322.0, 647.35/647.096, xi_array[i]) << std::endl << std::endl;

  return 0;
}
