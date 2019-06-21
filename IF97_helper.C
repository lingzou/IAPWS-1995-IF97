#include <math.h>
#include <stdio.h>

#include "IF97_helper.h"


int locateRegion_from_pT(double p, double T)
{
  if (T < IF97_T_MIN)
  {
    fprintf(stderr, "%s", "Out of range: T < IF97_T_MIN!\n");
    return -1;
  }

  if (T > IF97_T_MAX)
  {
    fprintf(stderr, "%s", "Out of range: T > IF97_T_MAX!\n");
    return -2;
  }

  if (p > IF97_P_MAX)
  {
    fprintf(stderr, "%s", "Out of range: p > IF97_P_MAX!\n");
    return -3;
  }

  if ((p > 50.0e6) && (T > IF97_T_25))
  {
    fprintf(stderr, "%s", "Out of range: (p > 50.0e6) && (T > IF97_T_25)!\n");
    return -4;
  }

  if (T <= IF97_T_13)
  {
    double ps = R4_p_sat_from_T(T);
    return (p > ps) ? 1 : 2;
  }
  else
  {
    double p_bc = B23_p_from_T(T);
    return (p > p_bc) ? 3 : 2;
  }
}

unsigned int
find_T_lower_bound(double T)
{
  if (T < 646.15)
    return (unsigned int)(T - 623.15);
  else if (T < 647.05)
    return (unsigned int)((T - 646.15) / 0.05) + 23;
  else if (T < 647.09)
    return (unsigned int)((T - 647.05) / 0.01) + 41;
  else if (T < 647.095)
    return (unsigned int)((T - 647.09) / 0.001) + 45;
  else if (T < 647.0955)
    return 50;
  else
    return 51;
}
/*
double rho_l_sat_from_T(double T)
{
  if ((T < IF97_T_MIN) || (T > T_CRIT))
  {
    fprintf(stderr, "%s", "Temperature is out of bound!\n");
    exit(1);
  }
  else if (T <= 623.15)
    return 1.0 / R1_specific_volume(R4_p_sat_from_T(T), T);
  else
    return R3_rho_l_sat_from_T_ITER(T);
}

double rho_g_sat_from_T(double T)
{
  if ((T < IF97_T_MIN) || (T > T_CRIT))
  {
    fprintf(stderr, "%s", "Temperature is out of bound!\n");
    exit(1);
  }
  else if (T <= 623.15)
    return 1.0 / R2_specific_volume(R4_p_sat_from_T(T), T);
  else
    return R3_rho_g_sat_from_T_ITER(T);
}*/

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
  for (int i = 0; i < 51; i++)
  {
    double bracket_size = (i < 50) ? 8.0 : 2.0;
    double T = R3_T_list[i];
    double ps = R4_p_sat_from_T(T);
    double rho_l_min = R3_rho_l_sat_guess[i] - bracket_size;
    double rho_l_max = R3_rho_l_sat_guess[i] + bracket_size;
    double rho_g_min = R3_rho_g_sat_guess[i] - bracket_size;
    double rho_g_max = R3_rho_g_sat_guess[i] + bracket_size;

    double delta_l_min = rho_l_min / RHO_CRIT;
    double delta_l_max = rho_l_max / RHO_CRIT;
    double delta_g_min = rho_g_min / RHO_CRIT;
    double delta_g_max = rho_g_max / RHO_CRIT;
    double tau = T_CRIT / T;

    double val_l_min = ps / R_GAS / T / rho_l_min - delta_l_min * R3_phi_delta(delta_l_min, tau);
    double val_l_max = ps / R_GAS / T / rho_l_max - delta_l_max * R3_phi_delta(delta_l_max, tau);

    double val_g_min = ps / R_GAS / T / rho_g_min - delta_g_min * R3_phi_delta(delta_g_min, tau);
    double val_g_max = ps / R_GAS / T / rho_g_max - delta_g_max * R3_phi_delta(delta_g_max, tau);

    // Iterate to find rho_l_sat
    double rho_error = 1.0;
    double min = rho_l_min;
    double max = rho_l_max;
    double rho_l_find = 0.0;
    while (rho_error > 1.0e-9)
    {
      rho_l_find = 0.5 * (min + max);
      double delta = rho_l_find / RHO_CRIT;
      double val = ps / R_GAS / T / rho_l_find - delta * R3_phi_delta(delta, tau);

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
      double delta = rho_g_find / RHO_CRIT;
      double val = ps / R_GAS / T / rho_g_find - delta * R3_phi_delta(delta, tau);

      if (val > 0.0) min = rho_g_find;
      else           max = rho_g_find;

      rho_error = max - min;
    }
    /*
    double LEFT = ps / R_GAS / T * (1.0 / rho_g_find - 1.0 / rho_l_find);
    double RIGHT = R3_phi(rho_l_find / RHO_CRIT, tau) - R3_phi(rho_g_find / RHO_CRIT, tau);
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
      double delta = rho / RHO_CRIT;
      double val = ps / R_GAS / T / rho - delta * R3_phi_delta(delta, tau);
      fprintf (ptr_File, "%20.8e%20.8e\n", rho, val);
    }
    for (int j = 0; j < 51; j++)
    {
      double rho = rho_g_max + (rho_l_min - rho_g_max) / 50 * j;
      double delta = rho / RHO_CRIT;
      double val = ps / R_GAS / T / rho - delta * R3_phi_delta(delta, tau);
      fprintf (ptr_File, "%20.8e%20.8e\n", rho, val);
    }
    for (int j = 0; j < 51; j++)
    {
      double rho = rho_l_min + (rho_l_max - rho_l_min) / 50 * j;
      double delta = rho / RHO_CRIT;
      double val = ps / R_GAS / T / rho - delta * R3_phi_delta(delta, tau);
      fprintf (ptr_File, "%20.8e%20.8e\n", rho, val);
    }

    fprintf (ptr_File, "%20.8e%20.8e\n", rho_l_find,
      ps / R_GAS / T / rho_l_find - (rho_l_find / RHO_CRIT) * R3_phi_delta(rho_l_find / RHO_CRIT, tau));

    fprintf (ptr_File, "%20.8e%20.8e\n", rho_g_find,
      ps / R_GAS / T / rho_g_find - (rho_g_find / RHO_CRIT) * R3_phi_delta(rho_g_find / RHO_CRIT, tau));

    fclose(ptr_File);*/
  }
  fprintf (ptr_sat_line_File, "%20.8e%20.8e%20.8e%20.8e%20.8e%20.8e%20.8e%20.8e%20.8e%20.8e%20.8e%20.8e%20.8e%20.8e%20.8e%20.8e%20.8e%20.8e\n",
    T_CRIT, P_CRIT,
    RHO_CRIT, 1.0/RHO_CRIT, R3_specific_int_energy(RHO_CRIT, T_CRIT), R3_specific_entropy(RHO_CRIT, T_CRIT),
    R3_specific_enthalpy(RHO_CRIT, T_CRIT), R3_cv(RHO_CRIT, T_CRIT), R3_cp(RHO_CRIT, T_CRIT), R3_sound_speed(RHO_CRIT, T_CRIT),
    RHO_CRIT, 1.0/RHO_CRIT, R3_specific_int_energy(RHO_CRIT, T_CRIT), R3_specific_entropy(RHO_CRIT, T_CRIT),
    R3_specific_enthalpy(RHO_CRIT, T_CRIT), R3_cv(RHO_CRIT, T_CRIT), R3_cp(RHO_CRIT, T_CRIT), R3_sound_speed(RHO_CRIT, T_CRIT)
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
    double ps = R4_p_sat_from_T(T);

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

double R3_rho_l_sat_from_T_ITER(double T)
{
  if (T > 647.0955)
  {
    double low = R3_rho_l_sat_exact[51];
    double high = R3_rho_l_sat_exact[52];
    return (T - 647.0955) / 0.0005 * (high - low) + low;
  }
  else
  {
    unsigned int i = find_T_lower_bound(T);
    double max = R3_rho_l_sat_exact[i] + 0.1;
    double min = R3_rho_l_sat_exact[i+1] - 0.1;
    double rho_error = 1.0;

    double rho_l_find = 0.0;
    double tau = T_CRIT / T;
    double ps = R4_p_sat_from_T(T);
    while (rho_error > 1.0e-9)
    {
      rho_l_find = 0.5 * (min + max);
      double delta = rho_l_find / RHO_CRIT;
      double val = ps / R_GAS / T / rho_l_find - delta * R3_phi_delta(delta, tau);

      if (val > 0.0) min = rho_l_find;
      else           max = rho_l_find;

      rho_error = max - min;
    }

    return rho_l_find;
  }
}

double R3_rho_g_sat_from_T_ITER(double T)
{
  if (T > 647.0955)
  {
    double low = R3_rho_g_sat_exact[51];
    double high = R3_rho_g_sat_exact[52];
    return (T - 647.0955) / 0.0005 * (high - low) + low;
  }
  else
  {
    unsigned int i = find_T_lower_bound(T);
    double min = R3_rho_g_sat_exact[i] - 0.1;
    double max = R3_rho_g_sat_exact[i+1] + 0.1;
    double rho_error = 1.0;

    double rho_g_find = 0.0;
    double tau = T_CRIT / T;
    double ps = R4_p_sat_from_T(T);
    while (rho_error > 1.0e-9)
    {
      rho_g_find = 0.5 * (min + max);
      double delta = rho_g_find / RHO_CRIT;
      double val = ps / R_GAS / T / rho_g_find - delta * R3_phi_delta(delta, tau);

      if (val > 0.0) min = rho_g_find;
      else           max = rho_g_find;

      rho_error = max - min;
    }

    return rho_g_find;
  }
}

double R3_rho_from_p_T_ITER(double p, double T)
{
  double p23 = B23_p_from_T(T);
  double rho_min = 0.0, rho_max = 0.0;
  if (p < P_CRIT)
  {
    double Ts = R4_T_sat_from_p(p);
    if (T > Ts) // vapor phase
    {
      rho_min = 1.0 / R2_specific_volume(p23, T);
      rho_max = R3_rho_g_sat_from_T_ITER(Ts);
    }
    else
    {
      rho_min = R3_rho_l_sat_from_T_ITER(Ts);
      rho_max = 1.0 / R1_specific_volume(p, IF97_T_13);
    }
  }
  else
  {
    rho_min = 1.0 / R2_specific_volume(p23, T);
    rho_max = 1.0 / R1_specific_volume(p, IF97_T_13);
  }

  double tau = T_CRIT / T;
  double rho_find = 0.0;
  double rho_error = 1.0;
  while (rho_error > 1.0e-9)
  {
    rho_find = 0.5 * (rho_min + rho_max);
    double delta = rho_find / RHO_CRIT;
    double p_guess = R3_p(rho_find, T);

    if (p_guess > p) rho_max = rho_find;
    else             rho_min = rho_find;

    rho_error = fabs((rho_max - rho_min) / rho_find);
  }

  return rho_find;
}

double R3_specific_volume_from_pT(double p, double T) { return 1.0 / R3_rho_from_p_T_ITER(p, T); }

void R3_rho_T_x_from_p_h_ITER(double p, double h, double &rho, double &T, double &x)
{
  double T_min = IF97_T_13;
  double T_max = B23_T_from_p(p);

  if (p < P_CRIT)
  {
    double Ts = R4_T_sat_from_p(p);
    double rho_l_sat = R3_rho_l_sat_from_T_ITER(Ts);
    double h_l_sat = R3_specific_enthalpy(rho_l_sat, Ts);
    double rho_g_sat = R3_rho_g_sat_from_T_ITER(Ts);
    double h_g_sat = R3_specific_enthalpy(rho_g_sat, Ts);

    if ((h >= h_l_sat) && (h <= h_g_sat))
    {
      T = Ts;
      x = (h - h_l_sat) / (h_g_sat - h_l_sat);
      rho = rho_l_sat * (1.0 - x) + rho_g_sat * x;
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
        rho = R3_rho_from_p_T_ITER(p, T);
        double h_find = R3_specific_enthalpy(rho, T);

        if (h_find > h)   T_max = T;
        else              T_min = T;

        T_error = fabs((T_max - T_min) / T);
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
      rho = R3_rho_from_p_T_ITER(p, T);
      double h_find = R3_specific_enthalpy(rho, T);

      if (h_find > h)   T_max = T;
      else              T_min = T;

      T_error = fabs((T_max - T_min) / T);
    }

    x = -1.0;
  }
}

void R3_rho_T_x_from_p_s_ITER(double p, double s, double &rho, double &T, double &x)
{
  double T_min = IF97_T_13;
  double T_max = B23_T_from_p(p);

  if (p < P_CRIT)
  {
    double Ts = R4_T_sat_from_p(p);
    double rho_l_sat = R3_rho_l_sat_from_T_ITER(Ts);
    double s_l_sat = R3_specific_entropy(rho_l_sat, Ts);
    double rho_g_sat = R3_rho_g_sat_from_T_ITER(Ts);
    double s_g_sat = R3_specific_entropy(rho_g_sat, Ts);

    if ((s >= s_l_sat) && (s <= s_g_sat))
    {
      T = Ts;
      x = (s - s_l_sat) / (s_g_sat - s_l_sat);
      rho = rho_l_sat * (1.0 - x) + rho_g_sat * x;
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
        rho = R3_rho_from_p_T_ITER(p, T);
        double s_find = R3_specific_entropy(rho, T);

        if (s_find > s)   T_max = T;
        else              T_min = T;

        T_error = fabs((T_max - T_min) / T);
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
      rho = R3_rho_from_p_T_ITER(p, T);
      double s_find = R3_specific_entropy(rho, T);

      if (s_find > s)   T_max = T;
      else              T_min = T;

      T_error = fabs((T_max - T_min) / T);
    }

    x = -1.0;
  }
}

double R3_dp_ddelta(double delta, double tau)
{
  return RHO_CRIT * R_GAS * T_CRIT / tau * (2.0 * delta * R3_phi_delta(delta, tau) + delta * delta * R3_phi_delta_delta(delta, tau));
}

double R1_rho_from_pT(double p, double T) { return 1.0 / R1_specific_volume(p, T); }

double R1_drho_dp(double p, double T)
{
  double pi  = p / R1_pStar;
  double tau = R1_TStar / T;
  double gamma_pi = R1_gamma_pi(pi, tau);

  return -R1_gamma_pi_pi(pi, tau) / (R_GAS * T * gamma_pi * gamma_pi);
}

double R2_rho_from_pT(double p, double T) { return 1.0 / R2_specific_volume(p, T); }

double R2_drho_dp(double p, double T)
{
  double pi  = p / R2_pStar;
  double tau = R2_TStar / T;
  double denom = R2_gamma_0_pi(pi, tau) + R2_gamma_r_pi(pi, tau);

  return -(R2_gamma_0_pi_pi(pi, tau) + R2_gamma_r_pi_pi(pi, tau)) / (R_GAS * T * denom * denom);
}

double R5_rho_from_pT(double p, double T) { return 1.0 / R5_specific_volume(p, T); }

double R5_T_from_p_h_ITER(double p, double h)
{
  double T_min = IF97_T_25;
  double T_max = IF97_T_MAX;
  double T_find, T_error = 1.0;

  while (T_error > 1.0e-9)
  {
    T_find = 0.5 * (T_min + T_max);
    double h_find = R5_specific_enthalpy(p, T_find);

    if (h_find > h)   T_max = T_find;
    else              T_min = T_find;

    T_error = fabs((T_max - T_min) / T_find);
  }

  return T_find;
}

double R5_T_from_p_s_ITER(double p, double s)
{
  double T_min = IF97_T_25;
  double T_max = IF97_T_MAX;
  double T_find, T_error = 1.0;

  while (T_error > 1.0e-9)
  {
    T_find = 0.5 * (T_min + T_max);
    double s_find = R5_specific_entropy(p, T_find);

    if (s_find > s)   T_max = T_find;
    else              T_min = T_find;

    T_error = fabs((T_max - T_min) / T_find);
  }

  return T_find;
}
