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
  return 0;
}
