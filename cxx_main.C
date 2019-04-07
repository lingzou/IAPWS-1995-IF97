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

  return 0;
}
