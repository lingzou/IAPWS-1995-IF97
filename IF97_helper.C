#include <iostream>
#include <iomanip>

#include "IF97_helper.h"

int findRegion(double p, double T)
{
  if (T < Tmin)
  {
    std::cerr << "Out of range: T < Tmin" << std::endl;
    return -1;
  }

  if (T > Tmax2)
  {
    std::cerr << "Out of range: T > Tmax2" << std::endl;
    return -2;
  }

  if (p > Pmax)
  {
    std::cerr << "Out of range: p < Pmax" << std::endl;
    return -3;
  }

  if ((p > 50.0e6) && (T > Tmax))
  {
    std::cerr << "Out of range: (p > 50.0e6) && (T > Tmax)" << std::endl;
    return -4;
  }

  if (T <= T13)
  {
    double ps = p_sat_from_T(T);
    return (p > ps) ? 1 : 2;
  }
  else
  {
    double p_bc = B23_p_from_T(T);
    return (p > p_bc) ? 3 : 2;
  }
}

void genR3_sat_line()
{
  FILE * ptr_sat_line_File;
  ptr_sat_line_File = fopen("UnitTest/R3_sat_line.dat", "w");
  for (int i = 0; i < 54; i++)
  {
    double bracket_size = (i < 49) ? 10.0 : 5.0;
    if (i == 53) bracket_size = 2.0;
    double T = R3_T_list[i];
    double ps = p_sat_from_T(T);
    double rho_l_min = R3_rho_l_sat_guess[i] - bracket_size;
    double rho_l_max = R3_rho_l_sat_guess[i] + bracket_size;
    double rho_g_min = R3_rho_g_sat_guess[i] - bracket_size;
    double rho_g_max = R3_rho_g_sat_guess[i] + bracket_size;

    double delta_l_min = rho_l_min / Rhocrit;
    double delta_l_max = rho_l_max / Rhocrit;
    double delta_g_min = rho_g_min / Rhocrit;
    double delta_g_max = rho_g_max / Rhocrit;
    double tau = Tcrit / T;

    double val_l_min = ps / Rgas / T / rho_l_min - delta_l_min * R3_phi_delta(delta_l_min, tau);
    double val_l_max = ps / Rgas / T / rho_l_max - delta_l_max * R3_phi_delta(delta_l_max, tau);

    double val_g_min = ps / Rgas / T / rho_g_min - delta_g_min * R3_phi_delta(delta_g_min, tau);
    double val_g_max = ps / Rgas / T / rho_g_max - delta_g_max * R3_phi_delta(delta_g_max, tau);
/*
    std::cout << T << " " << val_l_min << " " << val_l_max << " " << val_g_min << " " << val_g_max;
    if (val_l_min * val_l_max < 0.0) std::cout << " OK ";
    else std::cout << " NO ";
    if (val_g_min * val_g_max < 0.0) std::cout << " OK ";
    else std::cout << " NO ";
    std::cout << std::endl;
*/
    //double rho_l_sat = 0.5 * (rho_l_min + rho_l_max);
    double rho_error = 1.0;
    double min = rho_l_min;
    double max = rho_l_max;
    double rho_l_find = 0.0;
    while (rho_error > 1.0e-9)
    {
      rho_l_find = 0.5 * (min + max);
      double delta = rho_l_find / Rhocrit;
      double val = ps / Rgas / T / rho_l_find - delta * R3_phi_delta(delta, tau);

      if (val > 0.0) min = rho_l_find;
      else           max = rho_l_find;

      rho_error = max - min;
    }

    //std::cout << std::setprecision(8) << std::setw(20) << T << " "
    //  << std::scientific << std::setprecision(8) << std::setw(20) << rho_l_find;


    rho_error = 1.0;
    min = rho_g_min;
    max = rho_g_max;
    double rho_g_find = 0.0;
    while (rho_error > 1.0e-9)
    {
      rho_g_find = 0.5 * (min + max);
      double delta = rho_g_find / Rhocrit;
      double val = ps / Rgas / T / rho_g_find - delta * R3_phi_delta(delta, tau);

      if (val > 0.0) min = rho_g_find;
      else           max = rho_g_find;

      rho_error = max - min;
    }
    //std::cout << " " << std::scientific << std::setprecision(8) << std::setw(20) << rho_g_find << std::endl;
    fprintf (ptr_sat_line_File, "%20.8e%20.8e%20.8e\n", T, rho_l_find, rho_g_find);

    FILE * ptr_File;
    std::string file_name = "output/" + std::to_string(i) + ".dat";
    ptr_File = fopen(file_name.c_str(), "w");
    for (int j = 0; j < 51; j++)
    {
      double rho = rho_g_min + (rho_g_max - rho_g_min) / 50 * j;
      double delta = rho / Rhocrit;
      double val = ps / Rgas / T / rho - delta * R3_phi_delta(delta, tau);
      fprintf (ptr_File, "%20.8e%20.8e\n", rho, val);
    }
    for (int j = 0; j < 51; j++)
    {
      double rho = rho_g_max + (rho_l_min - rho_g_max) / 50 * j;
      double delta = rho / Rhocrit;
      double val = ps / Rgas / T / rho - delta * R3_phi_delta(delta, tau);
      fprintf (ptr_File, "%20.8e%20.8e\n", rho, val);
    }
    for (int j = 0; j < 51; j++)
    {
      double rho = rho_l_min + (rho_l_max - rho_l_min) / 50 * j;
      double delta = rho / Rhocrit;
      double val = ps / Rgas / T / rho - delta * R3_phi_delta(delta, tau);
      fprintf (ptr_File, "%20.8e%20.8e\n", rho, val);
    }

    fprintf (ptr_File, "%20.8e%20.8e\n", rho_l_find,
      ps / Rgas / T / rho_l_find - (rho_l_find / Rhocrit) * R3_phi_delta(rho_l_find / Rhocrit, tau));

    fprintf (ptr_File, "%20.8e%20.8e\n", rho_g_find,
      ps / Rgas / T / rho_g_find - (rho_g_find / Rhocrit) * R3_phi_delta(rho_g_find / Rhocrit, tau));

    fclose(ptr_File);
  }
  fprintf (ptr_sat_line_File, "%20.8e%20.8e%20.8e", Tcrit, Rhocrit, Rhocrit);
  fclose(ptr_sat_line_File);
}

void genR4_sat_line()
{
  for (int i = 0; i < 351; i++)
  {
    double T = 273.15 + i;
    double ps = p_sat_from_T(T);

    std::cout << std::scientific << std::setprecision(8) << std::setw(20) << T
      << std::scientific << std::setprecision(8) << std::setw(20) << 1.0 / R1_specific_volume(ps, T)
      << std::scientific << std::setprecision(8) << std::setw(20) << 1.0 / R2_specific_volume(ps, T) << std::endl;
  }
}

double R3_rho_from_p_T_ITER(double p, double T)
{
  double p23 = B23_p_from_T(T);
  double rho_min = 1.0 / R2_specific_volume(p23, T);
  double rho_max = 1.0 / R1_specific_volume(100.0e6, T13);

  double tau = Tcrit / T;
  double rho_find = 0.0;
  double rho_error = 1.0;
  while (rho_error > 1.0e-9)
  {
    rho_find = 0.5 * (rho_min + rho_max);
    double delta = rho_find / Rhocrit;
    double p_guess = R3_p(rho_find, T);

    if (p_guess > p) rho_max = rho_find;
    else             rho_min = rho_find;

    rho_error = rho_max - rho_min;
  }

  return rho_find;
}
