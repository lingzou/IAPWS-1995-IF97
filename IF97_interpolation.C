#include <iostream>
#include <iomanip>
#include <cstdio>

#include "IF97_interpolation.h"

IF97_Interpolation::IF97_Interpolation()
{
  loadSatLineData();

  enum_to_vector_map[P] = &p_sat;

  enum_to_vector_map[RHO_L] = &rho_l_sat;
  enum_to_vector_map[V_L] = &v_l_sat;
  enum_to_vector_map[E_L] = &e_l_sat;
  enum_to_vector_map[S_L] = &s_l_sat;
  enum_to_vector_map[H_L] = &h_l_sat;
  enum_to_vector_map[CV_L] = &cv_l_sat;
  enum_to_vector_map[CP_L] = &cp_l_sat;
  enum_to_vector_map[C_L] = &c_l_sat;

  enum_to_vector_map[RHO_G] = &rho_g_sat;
  enum_to_vector_map[V_G] = &v_g_sat;
  enum_to_vector_map[E_G] = &e_g_sat;
  enum_to_vector_map[S_G] = &s_g_sat;
  enum_to_vector_map[H_G] = &h_g_sat;
  enum_to_vector_map[CV_G] = &cv_g_sat;
  enum_to_vector_map[CP_G] = &cp_g_sat;
  enum_to_vector_map[C_G] = &c_g_sat;
}

IF97_Interpolation::~IF97_Interpolation()
{
}

void IF97_Interpolation::loadSatLineData()
{
  double T, p;
  double rho_l, v_l, e_l, s_l, h_l, cv_l, cp_l, c_l;
  double rho_g, v_g, e_g, s_g, h_g, cv_g, cp_g, c_g;

  FILE * ptr_File;
  ptr_File = fopen("tables/sat_line.dat", "r");

  if (ptr_File == NULL)
    std::cerr << "File not found" << std::endl;

  char buffer[1000];
  fgets(buffer, 1000, ptr_File); // skip the 1st line
  while (fscanf(ptr_File, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
    &T, &p,
    &rho_l, &v_l, &e_l, &s_l, &h_l, &cv_l, &cp_l, &c_l,
    &rho_g, &v_g, &e_g, &s_g, &h_g, &cv_g, &cp_g, &c_g) != EOF)
  {
    T_sat.push_back(T);
    p_sat.push_back(p);

    rho_l_sat.push_back(rho_l);
    v_l_sat.push_back(v_l);
    e_l_sat.push_back(e_l);
    s_l_sat.push_back(s_l);
    h_l_sat.push_back(h_l);
    cv_l_sat.push_back(cv_l);
    cp_l_sat.push_back(cp_l);
    c_l_sat.push_back(c_l);

    rho_g_sat.push_back(rho_g);
    v_g_sat.push_back(v_g);
    e_g_sat.push_back(e_g);
    s_g_sat.push_back(s_g);
    h_g_sat.push_back(h_g);
    cv_g_sat.push_back(cv_g);
    cp_g_sat.push_back(cp_g);
    c_g_sat.push_back(c_g);
  }

  fclose(ptr_File);
}

double IF97_Interpolation::INTPL_property_from_T(double T, INTPL_PROPERTY property)
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
    double ratio = (T - T_sat[idx-1]) / (T_sat[idx] - T_sat[idx-1]);
    std::vector<double> * property_vec = enum_to_vector_map[property];
    return (*property_vec)[idx-1] + ((*property_vec)[idx] - (*property_vec)[idx-1]) * ratio;
  }
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
