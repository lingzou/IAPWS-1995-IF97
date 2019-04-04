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
  return 0;
}
