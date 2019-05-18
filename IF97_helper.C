#include <iostream>
#include <iomanip>

#include "IF97_helper.h"

IF97_Interpolation global_IF97_INTPL;

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
  fprintf (ptr_sat_line_File, "%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s\n",
            "T [K]", "p_Sat [Pa]",
            "rho_l_sat [kg/m^3]",
            "v_l_sat [m^3/kg]",
            "e_l_sat [J/kg]",
            "s_l_sat [J/kg-K]",
            "h_l_sat [J/kg]",
            "cv_l_sat [J/kg-K]",
            "cp_l_sat [J/kg-K]",
            "c_l_sat [m/s]",
            "rho_g_sat [kg/m^3]",
            "v_g_sat [m^3/kg]",
            "e_g_sat [J/kg]",
            "s_g_sat [J/kg-K]",
            "h_g_sat [J/kg]",
            "cv_g_sat [J/kg-K]",
            "cp_g_sat [J/kg-K]",
            "c_g_sat [m/s]"
          );
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

    // Iterate to find rho_l_sat
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

    // Iterate to find rho_g_sat
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
    /*
    double LEFT = ps / Rgas / T * (1.0 / rho_g_find - 1.0 / rho_l_find);
    double RIGHT = R3_phi(rho_l_find / Rhocrit, tau) - R3_phi(rho_g_find / Rhocrit, tau);
    std::cout << "T = " << std::scientific << std::setprecision(8) << std::setw(20) << T
      << "LEFT = " << std::scientific << std::setprecision(8) << std::setw(20) << LEFT
      << "; RIGHT = " << std::scientific << std::setprecision(8) << std::setw(20) << RIGHT
      << "; DIFF = " << std::scientific << std::setprecision(8) << std::setw(20) << (LEFT - RIGHT)/LEFT
      << std::endl;*/

    fprintf (ptr_sat_line_File, "%20.8e%20.8e%20.8e%20.8e%20.8e%20.8e%20.8e%20.8e%20.8e%20.8e%20.8e%20.8e%20.8e%20.8e%20.8e%20.8e%20.8e%20.8e\n",
              T, ps,
              rho_l_find, 1.0/rho_l_find, R3_specific_int_energy(rho_l_find, T), R3_specific_entropy(rho_l_find, T),
              R3_specific_enthalpy(rho_l_find, T), R3_cv(rho_l_find, T), R3_cp(rho_l_find, T), R3_sound_speed(rho_l_find, T),
              rho_g_find, 1.0/rho_g_find, R3_specific_int_energy(rho_g_find, T), R3_specific_entropy(rho_g_find, T),
              R3_specific_enthalpy(rho_g_find, T), R3_cv(rho_g_find, T), R3_cp(rho_g_find, T), R3_sound_speed(rho_g_find, T)
            );

    /*
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

    fclose(ptr_File);*/
  }
  fprintf (ptr_sat_line_File, "%20.8e%20.8e%20.8e%20.8e%20.8e%20.8e%20.8e%20.8e%20.8e%20.8e%20.8e%20.8e%20.8e%20.8e%20.8e%20.8e%20.8e%20.8e\n",
    Tcrit, Pcrit,
    Rhocrit, 1.0/Rhocrit, R3_specific_int_energy(Rhocrit, Tcrit), R3_specific_entropy(Rhocrit, Tcrit),
    R3_specific_enthalpy(Rhocrit, Tcrit), R3_cv(Rhocrit, Tcrit), R3_cp(Rhocrit, Tcrit), R3_sound_speed(Rhocrit, Tcrit),
    Rhocrit, 1.0/Rhocrit, R3_specific_int_energy(Rhocrit, Tcrit), R3_specific_entropy(Rhocrit, Tcrit),
    R3_specific_enthalpy(Rhocrit, Tcrit), R3_cv(Rhocrit, Tcrit), R3_cp(Rhocrit, Tcrit), R3_sound_speed(Rhocrit, Tcrit)
  );
  fclose(ptr_sat_line_File);
}

void genR4_sat_line()
{
  FILE * ptr_sat_line_File;
  ptr_sat_line_File = fopen("UnitTest/R4_sat_line.dat", "w");
  fprintf (ptr_sat_line_File, "%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s\n",
            "T [K]", "p_Sat [Pa]",
            "rho_l_sat [kg/m^3]",
            "v_l_sat [m^3/kg]",
            "e_l_sat [J/kg]",
            "s_l_sat [J/kg-K]",
            "h_l_sat [J/kg]",
            "cv_l_sat [J/kg-K]",
            "cp_l_sat [J/kg-K]",
            "c_l_sat [m/s]",
            "rho_g_sat [kg/m^3]",
            "v_g_sat [m^3/kg]",
            "e_g_sat [J/kg]",
            "s_g_sat [J/kg-K]",
            "h_g_sat [J/kg]",
            "cv_g_sat [J/kg-K]",
            "cp_g_sat [J/kg-K]",
            "c_g_sat [m/s]"
          );
  for (int i = 0; i < 351; i++)
  {
    double T = 273.15 + i;
    double ps = p_sat_from_T(T);

    double v_l_sat = R1_specific_volume(ps, T);
    double rho_l_sat = 1.0 / v_l_sat;
    double v_g_sat = R2_specific_volume(ps, T);
    double rho_g_sat = 1.0 / v_g_sat;

    fprintf (ptr_sat_line_File, "%20.8e%20.8e%20.8e%20.8e%20.8e%20.8e%20.8e%20.8e%20.8e%20.8e%20.8e%20.8e%20.8e%20.8e%20.8e%20.8e%20.8e%20.8e\n",
              T, ps,
              rho_l_sat, v_l_sat, R1_specific_int_energy(ps, T), R1_specific_entropy(ps, T),
              R1_specific_enthalpy(ps, T), R1_cv(ps, T), R1_cp(ps, T), R1_sound_speed(ps, T),
              rho_g_sat, v_g_sat, R2_specific_int_energy(ps, T), R2_specific_entropy(ps, T),
              R2_specific_enthalpy(ps, T), R2_cv(ps, T), R2_cp(ps, T), R2_sound_speed(ps, T)
            );
  }
  fclose(ptr_sat_line_File);
}

double R3_rho_from_p_T_ITER(double p, double T)
{
  //IF97_Interpolation global_IF97_INTPL;

  double p23 = B23_p_from_T(T);
  double rho_min = 0.0, rho_max = 0.0;
  if (p < Pcrit)
  {
    double Ts = T_sat_from_p(p);
    if (T > Ts) // vapor phase
    {
      rho_min = 1.0 / R2_specific_volume(p23, T);
      rho_max = global_IF97_INTPL.INTPL_property_from_T(Ts, IF97_Interpolation::RHO_G);
    }
    else
    {
      rho_min = global_IF97_INTPL.INTPL_property_from_T(Ts, IF97_Interpolation::RHO_L);
      rho_max = 1.0 / R1_specific_volume(p, T13);
    }
  }
  else
  {
    rho_min = 1.0 / R2_specific_volume(p23, T);
    rho_max = 1.0 / R1_specific_volume(p, T13);
  }

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

    rho_error = std::abs((rho_max - rho_min) / rho_find);
  }

  return rho_find;
}

void R3_T_x_from_p_h_ITER(double p, double h, double &T, double &x)
{
  //IF97_Interpolation global_IF97_INTPL;
  double T_min = T13;
  double T_max = B23_T_from_p(p);

  if (p < Pcrit)
  {
    double Ts = T_sat_from_p(p);
    double h_l_sat = global_IF97_INTPL.INTPL_property_from_T(Ts, IF97_Interpolation::H_L);
    double h_g_sat = global_IF97_INTPL.INTPL_property_from_T(Ts, IF97_Interpolation::H_G);
/*
    std::cout << "p = " << p << std::endl;
    std::cout << "Ts = " << Ts << std::endl;
    std::cout << "h_l_sat = " << h_l_sat << std::endl;
    std::cout << "h_g_sat = " << h_g_sat << std::endl;*/

    if ((h >= h_l_sat) && (h <= h_g_sat))
    {
      T = Ts;
      x = (h - h_l_sat) / (h_g_sat - h_l_sat);
    }
    else
    {
      double h_min, h_max;
      if (h < h_l_sat)
      {
        h_min = R1_specific_enthalpy(p, T_min);
        h_max = h_l_sat;
        T_max = Ts;
        x = 0.0;
      }
      else
      {
        T_min = Ts;
        h_min = h_g_sat;
        h_max = R2_specific_enthalpy(p, T_max);
        x = 1.0;
      }

      double T_error = 1.0;
      while (T_error > 1.0e-9)
      {
        T = 0.5 * (T_min + T_max);
        double rho = R3_rho_from_p_T_ITER(p, T);
        double h_find = R3_specific_enthalpy(rho, T);

        if (h_find > h)   T_max = T;
        else              T_min = T;

        T_error = std::abs((T_max - T_min) / T);
      }
    }
  }
  else
  {
    double T_error = 1.0;
    double h_min = R1_specific_enthalpy(p, T_min);
    double h_max = R2_specific_enthalpy(p, T_max);
    while (T_error > 1.0e-9)
    {
      T = 0.5 * (T_min + T_max);
      double rho = R3_rho_from_p_T_ITER(p, T);
      double h_find = R3_specific_enthalpy(rho, T);

      if (h_find > h)   T_max = T;
      else              T_min = T;

      T_error = std::abs((T_max - T_min) / T);
    }

    x = -1.0;
  }
}

void R3_T_x_from_p_s_ITER(double p, double s, double &T, double &x)
{
  double T_min = T13;
  double T_max = B23_T_from_p(p);

  if (p < Pcrit)
  {
    double Ts = T_sat_from_p(p);
    double s_l_sat = global_IF97_INTPL.INTPL_property_from_T(Ts, IF97_Interpolation::S_L);
    double s_g_sat = global_IF97_INTPL.INTPL_property_from_T(Ts, IF97_Interpolation::S_G);

    if ((s >= s_l_sat) && (s <= s_g_sat))
    {
      T = Ts;
      x = (s - s_l_sat) / (s_g_sat - s_l_sat);
    }
    else
    {
      double s_min, s_max;
      if (s < s_l_sat)
      {
        s_min = R1_specific_entropy(p, T_min);
        s_max = s_l_sat;
        T_max = Ts;
        x = 0.0;
      }
      else
      {
        T_min = Ts;
        s_min = s_g_sat;
        s_max = R2_specific_entropy(p, T_max);
        x = 1.0;
      }

      double T_error = 1.0;
      while (T_error > 1.0e-9)
      {
        T = 0.5 * (T_min + T_max);
        double rho = R3_rho_from_p_T_ITER(p, T);
        double s_find = R3_specific_entropy(rho, T);

        if (s_find > s)   T_max = T;
        else              T_min = T;

        T_error = std::abs((T_max - T_min) / T);
      }
    }
  }
  else
  {
    double T_error = 1.0;
    double s_min = R1_specific_entropy(p, T_min);
    double s_max = R2_specific_entropy(p, T_max);
    while (T_error > 1.0e-9)
    {
      T = 0.5 * (T_min + T_max);
      double rho = R3_rho_from_p_T_ITER(p, T);
      double s_find = R3_specific_entropy(rho, T);

      if (s_find > s)   T_max = T;
      else              T_min = T;

      T_error = std::abs((T_max - T_min) / T);
    }

    x = -1.0;
  }
}

double R3_dp_ddelta(double delta, double tau)
{
  return Rhocrit * Rgas * Tcrit / tau * (2.0 * delta * R3_phi_delta(delta, tau) + delta * delta * R3_phi_delta_delta(delta, tau));
}

double R5_T_from_p_h_ITER(double p, double h)
{
  double T_min = Tmax;
  double T_max = Tmax2;
  double T_find, T_error = 1.0;

  while (T_error > 1.0e-9)
  {
    T_find = 0.5 * (T_min + T_max);
    double h_find = R5_specific_enthalpy(p, T_find);

    if (h_find > h)   T_max = T_find;
    else              T_min = T_find;

    T_error = std::abs((T_max - T_min) / T_find);
  }

  return T_find;
}

double R5_T_from_p_s_ITER(double p, double s)
{
  double T_min = Tmax;
  double T_max = Tmax2;
  double T_find, T_error = 1.0;

  while (T_error > 1.0e-9)
  {
    T_find = 0.5 * (T_min + T_max);
    double s_find = R5_specific_entropy(p, T_find);

    if (s_find > s)   T_max = T_find;
    else              T_min = T_find;

    T_error = std::abs((T_max - T_min) / T_find);
  }

  return T_find;
}
