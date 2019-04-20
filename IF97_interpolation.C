#include <iostream>
#include <iomanip>
#include <cstdio>

#include "IF97_interpolation.h"

IF97_Interpolation::IF97_Interpolation()
{
  loadSatLineData();
}

IF97_Interpolation::~IF97_Interpolation()
{
}

void IF97_Interpolation::loadSatLineData()
{
  double T, rho_l_sat_val, rho_g_sat_val;

  FILE * ptr_File;
  ptr_File = fopen("tables/rho_sat_from_T.dat", "r");

  if (ptr_File == NULL)
    std::cerr << "File not found" << std::endl;

  //int i = 0;
  while (fscanf(ptr_File, "%lf %lf %lf", &T, &rho_l_sat_val, &rho_g_sat_val) != EOF)
  {
    T_sat.push_back(T);
    rho_l_sat.push_back(rho_l_sat_val);
    rho_g_sat.push_back(rho_g_sat_val);
  }
  std::cout << "T_sat has " << T_sat.size() << " entries.\n";

  fclose(ptr_File);
}

double IF97_Interpolation::INTPL_rho_l_sat_from_T(double T)
{
  std::vector<double>::iterator low;
  unsigned int idx = std::lower_bound(T_sat.begin(), T_sat.end(), T) - T_sat.begin();
  if (idx == 0 || idx == T_sat.size())
  {
    std::cerr << "T = " << T << " is out of bound." << std::endl;
    //exit(1);
    return -999.9;
  }
  else
  {
    std::cout << "idx = " << idx << std::endl;
    std::cout << "T: " << T_sat[idx-1] << "; " << T << "; " << T_sat[idx] << std::endl;
    double ratio = (T - T_sat[idx-1]) / (T_sat[idx] - T_sat[idx-1]);
    double rho = rho_l_sat[idx-1] + (rho_l_sat[idx] - rho_l_sat[idx-1]) * ratio;
    std::cout << "rho_l: " << rho_l_sat[idx-1] << "; " << rho << "; " << rho_l_sat[idx] << std::endl;
    //return rho_l_sat[idx-1] + (rho_l_sat[idx] - rho_l_sat[idx-1]) * ratio;
    return rho;
  }
}

double IF97_Interpolation::INTPL_rho_g_sat_from_T(double T)
{
  std::vector<double>::iterator low;
  unsigned int idx = std::lower_bound(T_sat.begin(), T_sat.end(), T) - T_sat.begin();
  if (idx == 0 || idx == T_sat.size())
  {
    std::cerr << "T = " << T << " is out of bound." << std::endl;
    //exit(1);
    return -999.9;
  }
  else
  {
    std::cout << "idx = " << idx << std::endl;
    std::cout << "T: " << T_sat[idx-1] << "; " << T << "; " << T_sat[idx] << std::endl;
    double ratio = (T - T_sat[idx-1]) / (T_sat[idx] - T_sat[idx-1]);
    double rho = rho_g_sat[idx-1] + (rho_g_sat[idx] - rho_g_sat[idx-1]) * ratio;
    std::cout << "rho_g: " << rho_g_sat[idx-1] << "; " << rho << "; " << rho_g_sat[idx] << std::endl;
    //return rho_l_sat[idx-1] + (rho_l_sat[idx] - rho_l_sat[idx-1]) * ratio;
    return rho;
  }
}
