#include <iostream>
#include <iomanip>

#include "IF97.h"

int main()
{
  std::cout << "Hello, IF97" << std::endl;

  double T = 0.62315e3;
  double p = 0.165291643e8;

  std::cout << "p     = " << std::setprecision(9) << p << std::endl;
  std::cout << "T     = " << std::setprecision(9) << T << std::endl;
  std::cout << "p_cal = " << std::setprecision(9) << B23_p_from_T(T) << std::endl;
  std::cout << "T_cal = " << std::setprecision(9) << B23_T_from_p(p) << std::endl;

  std::cout << "v = " << std::scientific << std::setprecision(8) << std::setw(16)
      << R1_specific_volume(3.0e6, 300) << " "
      << R1_specific_volume(80.0e6, 300) << " "
      << R1_specific_volume(3.0e6, 500) << std::endl;
  std::cout << "h = " << std::scientific << std::setprecision(8) << std::setw(16)
      << R1_specific_enthalpy(3.0e6, 300) << " "
      << R1_specific_enthalpy(80.0e6, 300) << " "
      << R1_specific_enthalpy(3.0e6, 500) << std::endl;
  std::cout << "u = " << std::scientific << std::setprecision(8) << std::setw(16)
      << R1_specific_int_energy(3.0e6, 300) << " "
      << R1_specific_int_energy(80.0e6, 300) << " "
      << R1_specific_int_energy(3.0e6, 500) << std::endl;

  std::cout << "s = " << std::scientific << std::setprecision(8) << std::setw(16)
      << R1_specific_entropy(3.0e6, 300) << " "
      << R1_specific_entropy(80.0e6, 300) << " "
      << R1_specific_entropy(3.0e6, 500) << std::endl;

  std::cout << "cp = " << std::scientific << std::setprecision(8) << std::setw(16)
      << R1_cp(3.0e6, 300) << " "
      << R1_cp(80.0e6, 300) << " "
      << R1_cp(3.0e6, 500) << std::endl;

  std::cout << "w = " << std::scientific << std::setprecision(8) << std::setw(16)
      << R1_sound_speed(3.0e6, 300) << " "
      << R1_sound_speed(80.0e6, 300) << " "
      << R1_sound_speed(3.0e6, 500) << std::endl;
  std::cout << std::endl;

  // R1 backward T(p,h)
  std::cout << "T = " << std::scientific << std::setprecision(8) << std::setw(16)
      << R1_T_from_p_h(3.0e6, 500e3) << std::endl;
  std::cout << "T = " << std::scientific << std::setprecision(8) << std::setw(16)
      << R1_T_from_p_h(80.0e6, 500e3) << std::endl;
  std::cout << "T = " << std::scientific << std::setprecision(8) << std::setw(16)
      << R1_T_from_p_h(80.0e6, 1500e3) << std::endl;
  std::cout << std::endl;

  // R1 backward T(p,s)
  std::cout << "T = " << std::scientific << std::setprecision(8) << std::setw(16)
      << R1_T_from_p_s(3.0e6, 0.5e3) << std::endl;
  std::cout << "T = " << std::scientific << std::setprecision(8) << std::setw(16)
      << R1_T_from_p_s(80.0e6, 0.5e3) << std::endl;
  std::cout << "T = " << std::scientific << std::setprecision(8) << std::setw(16)
      << R1_T_from_p_s(80.0e6, 3e3) << std::endl;
  std::cout << std::endl;

  // R2
  std::cout << "v = " << std::scientific << std::setprecision(8) << std::setw(16)
      << R2_specific_volume(0.0035e6, 300) << " "
      << R2_specific_volume(0.0035e6, 700) << " "
      << R2_specific_volume(30.0e6, 700) << std::endl;
  std::cout << "h = " << std::scientific << std::setprecision(8) << std::setw(16)
      << R2_specific_enthalpy(0.0035e6, 300) << " "
      << R2_specific_enthalpy(0.0035e6, 700) << " "
      << R2_specific_enthalpy(30.0e6, 700) << std::endl;

  std::cout << "u = " << std::scientific << std::setprecision(8) << std::setw(16)
      << R2_specific_int_energy(0.0035e6, 300) << " "
      << R2_specific_int_energy(0.0035e6, 700) << " "
      << R2_specific_int_energy(30.0e6, 700) << std::endl;

  std::cout << "s = " << std::scientific << std::setprecision(8) << std::setw(16)
      << R2_specific_entropy(0.0035e6, 300) << " "
      << R2_specific_entropy(0.0035e6, 700) << " "
      << R2_specific_entropy(30.0e6, 700) << std::endl;

  std::cout << "cp = " << std::scientific << std::setprecision(8) << std::setw(16)
      << R2_cp(0.0035e6, 300) << " "
      << R2_cp(0.0035e6, 700) << " "
      << R2_cp(30.0e6, 700) << std::endl;

  std::cout << "w = " << std::scientific << std::setprecision(8) << std::setw(16)
      << R2_sound_speed(0.0035e6, 300) << " "
      << R2_sound_speed(0.0035e6, 700) << " "
      << R2_sound_speed(30.0e6, 700) << std::endl;
  std::cout << std::endl;

  // R2Meta
  std::cout << "v = " << std::scientific << std::setprecision(8) << std::setw(16)
      << R2Meta_specific_volume(1.0e6, 450) << " "
      << R2Meta_specific_volume(1.0e6, 440) << " "
      << R2Meta_specific_volume(1.5e6, 450) << std::endl;
  std::cout << "h = " << std::scientific << std::setprecision(8) << std::setw(16)
      << R2Meta_specific_enthalpy(1.0e6, 450) << " "
      << R2Meta_specific_enthalpy(1.0e6, 440) << " "
      << R2Meta_specific_enthalpy(1.5e6, 450) << std::endl;

  std::cout << "u = " << std::scientific << std::setprecision(8) << std::setw(16)
      << R2Meta_specific_int_energy(1.0e6, 450) << " "
      << R2Meta_specific_int_energy(1.0e6, 440) << " "
      << R2Meta_specific_int_energy(1.5e6, 450) << std::endl;

  std::cout << "s = " << std::scientific << std::setprecision(8) << std::setw(16)
      << R2Meta_specific_entropy(1.0e6, 450) << " "
      << R2Meta_specific_entropy(1.0e6, 440) << " "
      << R2Meta_specific_entropy(1.5e6, 450) << std::endl;

  std::cout << "cp = " << std::scientific << std::setprecision(8) << std::setw(16)
      << R2Meta_cp(1.0e6, 450) << " "
      << R2Meta_cp(1.0e6, 440) << " "
      << R2Meta_cp(1.5e6, 450) << std::endl;

  std::cout << "w = " << std::scientific << std::setprecision(8) << std::setw(16)
      << R2Meta_sound_speed(1.0e6, 450) << " "
      << R2Meta_sound_speed(1.0e6, 440) << " "
      << R2Meta_sound_speed(1.5e6, 450) << std::endl;
  std::cout << std::endl;

  std::cout << "p     = " << std::setprecision(9) << 100e6 << std::endl;
  std::cout << "h     = " << std::setprecision(9) << 0.3516004323e7 << std::endl;
  std::cout << "p_cal = " << std::setprecision(9) << B2bc_p_from_h(0.3516004323e7) << std::endl;
  std::cout << "h_cal = " << std::setprecision(9) << B2bc_h_from_p(100e6) << std::endl;

  std::cout << std::endl;
  // R2 backward (p, h)
  std::cout << "T     = " << std::setprecision(8) << R2a_T_from_p_h(0.001e6, 3.0e6) << std::endl;
  std::cout << "T     = " << std::setprecision(8) << R2a_T_from_p_h(3.0e6, 3.0e6) << std::endl;
  std::cout << "T     = " << std::setprecision(8) << R2a_T_from_p_h(3.0e6, 4.0e6) << std::endl << std::endl;

  std::cout << "T     = " << std::setprecision(8) << R2b_T_from_p_h(5.0e6, 3.5e6) << std::endl;
  std::cout << "T     = " << std::setprecision(8) << R2b_T_from_p_h(5.0e6, 4.0e6) << std::endl;
  std::cout << "T     = " << std::setprecision(8) << R2b_T_from_p_h(25.0e6, 3.5e6) << std::endl << std::endl;

  std::cout << "T     = " << std::setprecision(8) << R2c_T_from_p_h(40.0e6, 2.7e6) << std::endl;
  std::cout << "T     = " << std::setprecision(8) << R2c_T_from_p_h(60.0e6, 2.7e6) << std::endl;
  std::cout << "T     = " << std::setprecision(8) << R2c_T_from_p_h(60.0e6, 3.2e6) << std::endl << std::endl;

  // R2 backward (p, s)
  std::cout << "T     = " << std::setprecision(8) << R2a_T_from_p_s(0.1e6, 7.5e3) << std::endl;
  std::cout << "T     = " << std::setprecision(8) << R2a_T_from_p_s(0.1e6, 8.0e3) << std::endl;
  std::cout << "T     = " << std::setprecision(8) << R2a_T_from_p_s(2.5e6, 8.0e3) << std::endl << std::endl;

  std::cout << "T     = " << std::setprecision(8) << R2b_T_from_p_s(8.0e6, 6.0e3) << std::endl;
  std::cout << "T     = " << std::setprecision(8) << R2b_T_from_p_s(8.0e6, 7.5e3) << std::endl;
  std::cout << "T     = " << std::setprecision(8) << R2b_T_from_p_s(90.0e6, 6.0e3) << std::endl << std::endl;

  std::cout << "T     = " << std::setprecision(8) << R2c_T_from_p_s(20.0e6, 5.75e3) << std::endl;
  std::cout << "T     = " << std::setprecision(8) << R2c_T_from_p_s(80.0e6, 5.25e3) << std::endl;
  std::cout << "T     = " << std::setprecision(8) << R2c_T_from_p_s(80.0e6, 5.75e3) << std::endl << std::endl;

  // R3
  std::cout << "v = " << std::scientific << std::setprecision(8) << std::setw(16)
      << R3_p(500, 650) << " "
      << R3_p(200, 650) << " "
      << R3_p(500, 750) << std::endl;
  std::cout << "h = " << std::scientific << std::setprecision(8) << std::setw(16)
      << R3_specific_enthalpy(500, 650) << " "
      << R3_specific_enthalpy(200, 650) << " "
      << R3_specific_enthalpy(500, 750) << std::endl;

  std::cout << "u = " << std::scientific << std::setprecision(8) << std::setw(16)
      << R3_specific_int_energy(500, 650) << " "
      << R3_specific_int_energy(200, 650) << " "
      << R3_specific_int_energy(500, 750) << std::endl;

  std::cout << "s = " << std::scientific << std::setprecision(8) << std::setw(16)
      << R3_specific_entropy(500, 650) << " "
      << R3_specific_entropy(200, 650) << " "
      << R3_specific_entropy(500, 750) << std::endl;

  std::cout << "cp = " << std::scientific << std::setprecision(8) << std::setw(16)
      << R3_cp(500, 650) << " "
      << R3_cp(200, 650) << " "
      << R3_cp(500, 750) << std::endl;

  std::cout << "w = " << std::scientific << std::setprecision(8) << std::setw(16)
      << R3_sound_speed(500, 650) << " "
      << R3_sound_speed(200, 650) << " "
      << R3_sound_speed(500, 750) << std::endl;
  std::cout << std::endl;

  std::cout << "p_sat     = " << std::setprecision(8) << p_sat_from_T(300.0) << std::endl;
  std::cout << "p_sat     = " << std::setprecision(8) << p_sat_from_T(500.0) << std::endl;
  std::cout << "p_sat     = " << std::setprecision(8) << p_sat_from_T(600.0) << std::endl << std::endl;

  std::cout << "T_sat     = " << std::setprecision(8) << T_sat_from_p(0.1e6) << std::endl;
  std::cout << "T_sat     = " << std::setprecision(8) << T_sat_from_p(1.0e6) << std::endl;
  std::cout << "T_sat     = " << std::setprecision(8) << T_sat_from_p(10.0e6) << std::endl << std::endl;

  // R5
  // R5
  std::cout << "v = " << std::scientific << std::setprecision(8) << std::setw(16)
      << R5_specific_volume(0.5e6, 1500) << " "
      << R5_specific_volume(30.0e6, 1500) << " "
      << R5_specific_volume(30.0e6, 2000) << std::endl;
  std::cout << "h = " << std::scientific << std::setprecision(8) << std::setw(16)
      << R5_specific_enthalpy(0.5e6, 1500) << " "
      << R5_specific_enthalpy(30.0e6, 1500) << " "
      << R5_specific_enthalpy(30.0e6, 2000) << std::endl;

  std::cout << "u = " << std::scientific << std::setprecision(8) << std::setw(16)
      << R5_specific_int_energy(0.5e6, 1500) << " "
      << R5_specific_int_energy(30.0e6, 1500) << " "
      << R5_specific_int_energy(30.0e6, 2000) << std::endl;

  std::cout << "s = " << std::scientific << std::setprecision(8) << std::setw(16)
      << R5_specific_entropy(0.5e6, 1500) << " "
      << R5_specific_entropy(30.0e6, 1500) << " "
      << R5_specific_entropy(30.0e6, 2000) << std::endl;

  std::cout << "cp = " << std::scientific << std::setprecision(8) << std::setw(16)
      << R5_cp(0.5e6, 1500) << " "
      << R5_cp(30.0e6, 1500) << " "
      << R5_cp(30.0e6, 2000) << std::endl;

  std::cout << "w = " << std::scientific << std::setprecision(8) << std::setw(16)
      << R5_sound_speed(0.5e6, 1500) << " "
      << R5_sound_speed(30.0e6, 1500) << " "
      << R5_sound_speed(30.0e6, 2000) << std::endl;
  std::cout << std::endl;
  
  return 0;
}
