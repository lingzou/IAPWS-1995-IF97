#include <cstdio>
#include <sstream>
#include <iomanip>
#include <iostream>

#include "unit_tests.h"
#include "IF97_helper.h"
#include "SurfaceTension.h"
#include "Viscosity.h"
#include "ThermalConductivity.h"

void Unit_Test_All()
{
  Unit_Test_B23();

  Unit_Test_R1();
  Unit_Test_R1_PH();
  Unit_Test_R1_PS();

  Unit_Test_R2();
  Unit_Test_R2Meta();
  Unit_Test_B2bc();
  Unit_Test_R2_PH();
  Unit_Test_R2_PS();

  Unit_Test_R3();

  Unit_Test_R4();

  Unit_Test_R5();

  genR3_sat_line();
  genR4_sat_line();

  Unit_Test_R3_rho_pT_ITER();
  Unit_Test_R3_Tx_ph_ITER();
  Unit_Test_R3_Tx_ps_ITER();

  Unit_Test_R5_T_ph_ITER();
  Unit_Test_R5_T_ps_ITER();

  Unit_Test_SurfTension();
  Unit_Test_Viscosity();

  Unit_Test_ThermCond_NoEnhancement();
  Unit_Test_ThermCond_R1();
  Unit_Test_ThermCond_R2();
  Unit_Test_ThermCond_R3();
}

void Unit_Test(std::string str)
{
  if(str.compare("B23") == 0)
    Unit_Test_B23();
  else if(str.compare("R1") == 0)
    Unit_Test_R1();
  else if(str.compare("R1_PH") == 0)
    Unit_Test_R1_PH();
  else if(str.compare("R1_PS") == 0)
    Unit_Test_R1_PS();
  else if(str.compare("R2") == 0)
    Unit_Test_R2();
  else if(str.compare("R2Meta") == 0)
    Unit_Test_R2Meta();
  else if(str.compare("B2bc") == 0)
    Unit_Test_B2bc();
  else if(str.compare("R2_PH") == 0)
    Unit_Test_R2_PH();
  else if(str.compare("R2_PS") == 0)
    Unit_Test_R2_PS();
  else if(str.compare("R3") == 0)
    Unit_Test_R3();
  else if(str.compare("R4") == 0)
    Unit_Test_R4();
  else if(str.compare("R5") == 0)
    Unit_Test_R5();
  else
  {
    std::cerr << "I do not understand this unit test: " << str << std::endl;
    exit(1);
  }
}

void Unit_Test_B23()
{
  FILE * ptr_File;
  ptr_File = fopen("UnitTest/B23.dat", "w");

  std::ostringstream out_ss;

  double T = 0.62315e3;
  double p = 0.165291643e8;

  out_ss << "#     p[Pa]               T[K]" << std::endl;
  out_ss << std::scientific << std::setprecision(8)
          << std::setw(20) << B23_p_from_T(T)
          << std::setw(20) << B23_T_from_p(p) << std::endl;

  fprintf(ptr_File, "%s", out_ss.str().c_str());
  fclose(ptr_File);
}

void Unit_Test_R1()
{
  FILE * ptr_File;
  ptr_File = fopen("UnitTest/R1.dat", "w");

  std::ostringstream out_ss;

  out_ss << "# p[Pa]     =                3.0e6              80.0e6               3.0e6" << std::endl;
  out_ss << "# T[K]      =                  300                 300                 500" << std::endl;
  out_ss << "v [m3/kg]   = " << std::scientific << std::setprecision(8)
      << std::setw(20) << R1_specific_volume(3.0e6, 300)
      << std::setw(20) << R1_specific_volume(80.0e6, 300)
      << std::setw(20) << R1_specific_volume(3.0e6, 500) << std::endl;

  out_ss << "h [J/kg]    = " << std::scientific << std::setprecision(8)
      << std::setw(20) << R1_specific_enthalpy(3.0e6, 300)
      << std::setw(20) << R1_specific_enthalpy(80.0e6, 300)
      << std::setw(20) << R1_specific_enthalpy(3.0e6, 500) << std::endl;

  out_ss << "u [J/kg]    = " << std::scientific << std::setprecision(8)
      << std::setw(20) << R1_specific_int_energy(3.0e6, 300)
      << std::setw(20) << R1_specific_int_energy(80.0e6, 300)
      << std::setw(20) << R1_specific_int_energy(3.0e6, 500) << std::endl;

  out_ss << "s [J/kg-K]  = " << std::scientific << std::setprecision(8)
      << std::setw(20) << R1_specific_entropy(3.0e6, 300)
      << std::setw(20) << R1_specific_entropy(80.0e6, 300)
      << std::setw(20) << R1_specific_entropy(3.0e6, 500) << std::endl;

  out_ss << "cp [J/kg-K] = " << std::scientific << std::setprecision(8)
      << std::setw(20) << R1_cp(3.0e6, 300)
      << std::setw(20) << R1_cp(80.0e6, 300)
      << std::setw(20) << R1_cp(3.0e6, 500) << std::endl;

  out_ss << "w [m/s]     = " << std::scientific << std::setprecision(8)
      << std::setw(20) << R1_sound_speed(3.0e6, 300)
      << std::setw(20) << R1_sound_speed(80.0e6, 300)
      << std::setw(20) << R1_sound_speed(3.0e6, 500) << std::endl;

  fprintf(ptr_File, "%s", out_ss.str().c_str());
  fclose(ptr_File);
}

void Unit_Test_R1_PH()
{
  // R1 backward T(p,h)
  FILE * ptr_File;
  ptr_File = fopen("UnitTest/R1_PH.dat", "w");

  std::ostringstream out_ss;

  out_ss << "T = " << std::scientific << std::setprecision(8) << std::setw(20)
      << R1_T_from_p_h(3.0e6, 500e3) << std::endl;
  out_ss << "T = " << std::scientific << std::setprecision(8) << std::setw(20)
      << R1_T_from_p_h(80.0e6, 500e3) << std::endl;
  out_ss << "T = " << std::scientific << std::setprecision(8) << std::setw(20)
      << R1_T_from_p_h(80.0e6, 1500e3) << std::endl;

  fprintf(ptr_File, "%s", out_ss.str().c_str());
  fclose(ptr_File);
}

void Unit_Test_R1_PS()
{
  // R1 backward T(p,s)
  FILE * ptr_File;
  ptr_File = fopen("UnitTest/R1_PS.dat", "w");

  std::ostringstream out_ss;

  out_ss << "T = " << std::scientific << std::setprecision(8) << std::setw(20)
      << R1_T_from_p_s(3.0e6, 0.5e3) << std::endl;
  out_ss << "T = " << std::scientific << std::setprecision(8) << std::setw(20)
      << R1_T_from_p_s(80.0e6, 0.5e3) << std::endl;
  out_ss << "T = " << std::scientific << std::setprecision(8) << std::setw(20)
      << R1_T_from_p_s(80.0e6, 3e3) << std::endl;

  fprintf(ptr_File, "%s", out_ss.str().c_str());
  fclose(ptr_File);
}

void Unit_Test_R2()
{
  // R2
  FILE * ptr_File;
  ptr_File = fopen("UnitTest/R2.dat", "w");

  std::ostringstream out_ss;

  out_ss << "# p[Pa]     =             0.0035e6            0.0035e6              30.0e6" << std::endl;
  out_ss << "# T[K]      =                  300                 700                 700" << std::endl;

  out_ss << "v [m3/kg]   = " << std::scientific << std::setprecision(8)
      << std::setw(20) << R2_specific_volume(0.0035e6, 300)
      << std::setw(20) << R2_specific_volume(0.0035e6, 700)
      << std::setw(20) << R2_specific_volume(30.0e6, 700) << std::endl;
  out_ss << "h [J/kg]    = " << std::scientific << std::setprecision(8)
      << std::setw(20) << R2_specific_enthalpy(0.0035e6, 300)
      << std::setw(20) << R2_specific_enthalpy(0.0035e6, 700)
      << std::setw(20) << R2_specific_enthalpy(30.0e6, 700) << std::endl;

  out_ss << "u [J/kg]    = " << std::scientific << std::setprecision(8)
      << std::setw(20) << R2_specific_int_energy(0.0035e6, 300)
      << std::setw(20) << R2_specific_int_energy(0.0035e6, 700)
      << std::setw(20) << R2_specific_int_energy(30.0e6, 700) << std::endl;

  out_ss << "s [J/kg-K]  = " << std::scientific << std::setprecision(8)
      << std::setw(20) << R2_specific_entropy(0.0035e6, 300)
      << std::setw(20) << R2_specific_entropy(0.0035e6, 700)
      << std::setw(20) << R2_specific_entropy(30.0e6, 700) << std::endl;

  out_ss << "cp [J/kg-K] = " << std::scientific << std::setprecision(8)
      << std::setw(20) << R2_cp(0.0035e6, 300)
      << std::setw(20) << R2_cp(0.0035e6, 700)
      << std::setw(20) << R2_cp(30.0e6, 700) << std::endl;

  out_ss << "w [m/s]     = " << std::scientific << std::setprecision(8)
      << std::setw(20) << R2_sound_speed(0.0035e6, 300)
      << std::setw(20) << R2_sound_speed(0.0035e6, 700)
      << std::setw(20) << R2_sound_speed(30.0e6, 700) << std::endl;

  fprintf(ptr_File, "%s", out_ss.str().c_str());
  fclose(ptr_File);
}

void Unit_Test_R2Meta()
{
  // R2Meta
  FILE * ptr_File;
  ptr_File = fopen("UnitTest/R2Meta.dat", "w");

  std::ostringstream out_ss;

  out_ss << "# p[Pa]     =                1.0e6               1.0e6               1.5e6" << std::endl;
  out_ss << "# T[K]      =                  450                 440                 450" << std::endl;

  out_ss << "v [m3/kg]   = " << std::scientific << std::setprecision(8)
      << std::setw(20) << R2Meta_specific_volume(1.0e6, 450)
      << std::setw(20) << R2Meta_specific_volume(1.0e6, 440)
      << std::setw(20) << R2Meta_specific_volume(1.5e6, 450) << std::endl;
  out_ss << "h [J/kg]    = " << std::scientific << std::setprecision(8)
      << std::setw(20) << R2Meta_specific_enthalpy(1.0e6, 450)
      << std::setw(20) << R2Meta_specific_enthalpy(1.0e6, 440)
      << std::setw(20) << R2Meta_specific_enthalpy(1.5e6, 450) << std::endl;

  out_ss << "u [J/kg]    = " << std::scientific << std::setprecision(8)
      << std::setw(20) << R2Meta_specific_int_energy(1.0e6, 450)
      << std::setw(20) << R2Meta_specific_int_energy(1.0e6, 440)
      << std::setw(20) << R2Meta_specific_int_energy(1.5e6, 450) << std::endl;

  out_ss << "s [J/kg-K]  = " << std::scientific << std::setprecision(8)
      << std::setw(20) << R2Meta_specific_entropy(1.0e6, 450)
      << std::setw(20) << R2Meta_specific_entropy(1.0e6, 440)
      << std::setw(20) << R2Meta_specific_entropy(1.5e6, 450) << std::endl;

  out_ss << "cp [J/kg-K] = " << std::scientific << std::setprecision(8)
      << std::setw(20) << R2Meta_cp(1.0e6, 450)
      << std::setw(20) << R2Meta_cp(1.0e6, 440)
      << std::setw(20) << R2Meta_cp(1.5e6, 450) << std::endl;

  out_ss << "w [m/s]     = " << std::scientific << std::setprecision(8)
      << std::setw(20) << R2Meta_sound_speed(1.0e6, 450)
      << std::setw(20) << R2Meta_sound_speed(1.0e6, 440)
      << std::setw(20) << R2Meta_sound_speed(1.5e6, 450) << std::endl;

  fprintf(ptr_File, "%s", out_ss.str().c_str());
  fclose(ptr_File);
}

void Unit_Test_B2bc()
{
  FILE * ptr_File;
  ptr_File = fopen("UnitTest/B2bc.dat", "w");

  std::ostringstream out_ss;

  out_ss << "#     p[Pa]               h[J/kg]" << std::endl;
  out_ss << std::scientific << std::setprecision(8)
          << std::setw(20) << B2bc_p_from_h(0.3516004323e7)
          << std::setw(20) << B2bc_h_from_p(100e6) << std::endl;

  fprintf(ptr_File, "%s", out_ss.str().c_str());
  fclose(ptr_File);
}

void Unit_Test_R2_PH()
{
  FILE * ptr_File;
  ptr_File = fopen("UnitTest/R2_PH.dat", "w");

  std::ostringstream out_ss;

  // R2 backward (p, h)
  out_ss << "T     = " << std::scientific << std::setprecision(8) << R2a_T_from_p_h(0.001e6, 3.0e6) << std::endl;
  out_ss << "T     = " << std::scientific << std::setprecision(8) << R2a_T_from_p_h(3.0e6, 3.0e6) << std::endl;
  out_ss << "T     = " << std::scientific << std::setprecision(8) << R2a_T_from_p_h(3.0e6, 4.0e6) << std::endl << std::endl;

  out_ss << "T     = " << std::scientific << std::setprecision(8) << R2b_T_from_p_h(5.0e6, 3.5e6) << std::endl;
  out_ss << "T     = " << std::scientific << std::setprecision(8) << R2b_T_from_p_h(5.0e6, 4.0e6) << std::endl;
  out_ss << "T     = " << std::scientific << std::setprecision(8) << R2b_T_from_p_h(25.0e6, 3.5e6) << std::endl << std::endl;

  out_ss << "T     = " << std::scientific << std::setprecision(8) << R2c_T_from_p_h(40.0e6, 2.7e6) << std::endl;
  out_ss << "T     = " << std::scientific << std::setprecision(8) << R2c_T_from_p_h(60.0e6, 2.7e6) << std::endl;
  out_ss << "T     = " << std::scientific << std::setprecision(8) << R2c_T_from_p_h(60.0e6, 3.2e6) << std::endl << std::endl;

  fprintf(ptr_File, "%s", out_ss.str().c_str());
  fclose(ptr_File);
}

void Unit_Test_R2_PS()
{
  FILE * ptr_File;
  ptr_File = fopen("UnitTest/R2_PS.dat", "w");

  std::ostringstream out_ss;

  // R2 backward (p, h)
  out_ss << "T     = " << std::scientific << std::setprecision(8) << R2a_T_from_p_s(0.1e6, 7.5e3) << std::endl;
  out_ss << "T     = " << std::scientific << std::setprecision(8) << R2a_T_from_p_s(0.1e6, 8.0e3) << std::endl;
  out_ss << "T     = " << std::scientific << std::setprecision(8) << R2a_T_from_p_s(2.5e6, 8.0e3) << std::endl << std::endl;

  out_ss << "T     = " << std::scientific << std::setprecision(8) << R2b_T_from_p_s(8.0e6, 6.0e3) << std::endl;
  out_ss << "T     = " << std::scientific << std::setprecision(8) << R2b_T_from_p_s(8.0e6, 7.5e3) << std::endl;
  out_ss << "T     = " << std::scientific << std::setprecision(8) << R2b_T_from_p_s(90.0e6, 6.0e3) << std::endl << std::endl;

  out_ss << "T     = " << std::scientific << std::setprecision(8) << R2c_T_from_p_s(20.0e6, 5.75e3) << std::endl;
  out_ss << "T     = " << std::scientific << std::setprecision(8) << R2c_T_from_p_s(80.0e6, 5.25e3) << std::endl;
  out_ss << "T     = " << std::scientific << std::setprecision(8) << R2c_T_from_p_s(80.0e6, 5.75e3) << std::endl << std::endl;

  fprintf(ptr_File, "%s", out_ss.str().c_str());
  fclose(ptr_File);
}

void Unit_Test_R3()
{
  FILE * ptr_File;
  ptr_File = fopen("UnitTest/R3.dat", "w");

  std::ostringstream out_ss;

  out_ss << "# rho[kg/m3]=                  500                 200                 500" << std::endl;
  out_ss << "# T[K]      =                  650                 650                 750" << std::endl;
  out_ss << "p [Pa]      = " << std::scientific << std::setprecision(8)
      << std::setw(20) << R3_p(500, 650)
      << std::setw(20) << R3_p(200, 650)
      << std::setw(20) << R3_p(500, 750) << std::endl;

  out_ss << "h [J/kg]    = " << std::scientific << std::setprecision(8)
      << std::setw(20) << R3_specific_enthalpy(500, 650)
      << std::setw(20) << R3_specific_enthalpy(200, 650)
      << std::setw(20) << R3_specific_enthalpy(500, 750) << std::endl;

  out_ss << "u [J/kg]    = " << std::scientific << std::setprecision(8)
      << std::setw(20) << R3_specific_int_energy(500, 650)
      << std::setw(20) << R3_specific_int_energy(200, 650)
      << std::setw(20) << R3_specific_int_energy(500, 750) << std::endl;

  out_ss << "s [J/kg-K]  = " << std::scientific << std::setprecision(8)
      << std::setw(20) << R3_specific_entropy(500, 650)
      << std::setw(20) << R3_specific_entropy(200, 650)
      << std::setw(20) << R3_specific_entropy(500, 750) << std::endl;

  out_ss << "cp [J/kg-K] = " << std::scientific << std::setprecision(8)
      << std::setw(20) << R3_cp(500, 650)
      << std::setw(20) << R3_cp(200, 650)
      << std::setw(20) << R3_cp(500, 750) << std::endl;

  out_ss << "w [m/s]     = " << std::scientific << std::setprecision(8)
      << std::setw(20) << R3_sound_speed(500, 650)
      << std::setw(20) << R3_sound_speed(200, 650)
      << std::setw(20) << R3_sound_speed(500, 750) << std::endl;

  fprintf(ptr_File, "%s", out_ss.str().c_str());
  fclose(ptr_File);
}

void Unit_Test_R4()
{
  FILE * ptr_File;
  ptr_File = fopen("UnitTest/R4.dat", "w");

  std::ostringstream out_ss;

  out_ss << "p_sat [Pa] = " << std::scientific << std::setprecision(8) << R4_p_sat_from_T(300.0) << std::endl;
  out_ss << "p_sat [Pa] = " << std::scientific << std::setprecision(8) << R4_p_sat_from_T(500.0) << std::endl;
  out_ss << "p_sat [Pa] = " << std::scientific << std::setprecision(8) << R4_p_sat_from_T(600.0) << std::endl << std::endl;

  out_ss << "T_sat [K] = " << std::scientific << std::setprecision(8) << R4_T_sat_from_p(0.1e6) << std::endl;
  out_ss << "T_sat [K] = " << std::scientific << std::setprecision(8) << R4_T_sat_from_p(1.0e6) << std::endl;
  out_ss << "T_sat [K] = " << std::scientific << std::setprecision(8) << R4_T_sat_from_p(10.0e6) << std::endl << std::endl;

  fprintf(ptr_File, "%s", out_ss.str().c_str());
  fclose(ptr_File);
}

void Unit_Test_R5()
{
  FILE * ptr_File;
  ptr_File = fopen("UnitTest/R5.dat", "w");

  std::ostringstream out_ss;

  out_ss << "# p[Pa]     =                0.5e6              30.0e6              30.0e6" << std::endl;
  out_ss << "# T[K]      =                 1500                1500                2000" << std::endl;
  out_ss << "v [m3/kg]   = " << std::scientific << std::setprecision(8)
      << std::setw(20) << R5_specific_volume(0.5e6, 1500)
      << std::setw(20) << R5_specific_volume(30.0e6, 1500)
      << std::setw(20) << R5_specific_volume(30.0e6, 2000) << std::endl;

  out_ss << "h [J/kg]    = " << std::scientific << std::setprecision(8)
      << std::setw(20) << R5_specific_enthalpy(0.5e6, 1500)
      << std::setw(20) << R5_specific_enthalpy(30.0e6, 1500)
      << std::setw(20) << R5_specific_enthalpy(30.0e6, 2000) << std::endl;

  out_ss << "u [J/kg]    = " << std::scientific << std::setprecision(8)
      << std::setw(20) << R5_specific_int_energy(0.5e6, 1500)
      << std::setw(20) << R5_specific_int_energy(30.0e6, 1500)
      << std::setw(20) << R5_specific_int_energy(30.0e6, 2000) << std::endl;

  out_ss << "s [J/kg-K]  = " << std::scientific << std::setprecision(8)
      << std::setw(20) << R5_specific_entropy(0.5e6, 1500)
      << std::setw(20) << R5_specific_entropy(30.0e6, 1500)
      << std::setw(20) << R5_specific_entropy(30.0e6, 2000) << std::endl;

  out_ss << "cp [J/kg-K] = " << std::scientific << std::setprecision(8)
      << std::setw(20) << R5_cp(0.5e6, 1500)
      << std::setw(20) << R5_cp(30.0e6, 1500)
      << std::setw(20) << R5_cp(30.0e6, 2000) << std::endl;

  out_ss << "w [m/s]     = " << std::scientific << std::setprecision(8)
      << std::setw(20) << R5_sound_speed(0.5e6, 1500)
      << std::setw(20) << R5_sound_speed(30.0e6, 1500)
      << std::setw(20) << R5_sound_speed(30.0e6, 2000) << std::endl;

  fprintf(ptr_File, "%s", out_ss.str().c_str());
  fclose(ptr_File);
}

void Unit_Test_R3_rho_pT_ITER()
{
  FILE * ptr_File;
  ptr_File = fopen("UnitTest/R3_rho_pT.dat", "w");

  double T_array[8] = {630.0, 637.0, 647.0, 626.16, 629.0, 697.0, 636.0, 854.0};
  double p_array[8] = {33.5e6, 20.e6, 21.0e6, 17.9e6, 17.37e6, 31.0e6, 97.0e6, 97.7e6};

  fprintf (ptr_File, "%20s%20s%20s\n", "T [K]", "p [Pa]", "rho [kg/m^3]");
  for (int i = 0; i < 8; i++)
    fprintf (ptr_File, "%20.8e%20.8e%20.8e\n", T_array[i], p_array[i], R3_rho_from_p_T_ITER(p_array[i], T_array[i]));

  fclose(ptr_File);
}

void Unit_Test_R3_Tx_ph_ITER()
{
  FILE * ptr_File;
  ptr_File = fopen("UnitTest/R3_Tx_ph.dat", "w");

  double h_array[8] = {1.64344943e6, 1.79130276e6, 2.48929822e6, 1.68821864e6, 2.56582491e6, 2.55197817e6, 1.61414193e6, 2.77581149e6};
  double p_array[8] = {33.5e6, 20.e6, 21.0e6, 17.9e6, 17.37e6, 31.0e6, 97.0e6, 97.7e6};

  fprintf (ptr_File, "%20s%20s%20s\n", "p [Pa]", "h [J/kg]", "T [K]");
  for (int i = 0; i < 8; i++)
  {
    double rho = 0.0, T = 0.0, x = 0.0;
    R3_rho_T_x_from_p_h_ITER(p_array[i], h_array[i], rho, T, x);
    fprintf (ptr_File, "%20.8e%20.8e%20.8e\n", p_array[i], h_array[i], T);
  }

  fclose(ptr_File);
}

void Unit_Test_R3_Tx_ps_ITER()
{
  FILE * ptr_File;
  ptr_File = fopen("UnitTest/R3_Tx_ps.dat", "w");

  double s_array[8] = {3.69009147e3, 3.95927258e3, 5.04168152e3, 3.80223220e3, 5.20290950e3, 5.05425538e3, 3.49976082e3, 5.06161454e3};
  double p_array[8] = {33.5e6, 20.e6, 21.0e6, 17.9e6, 17.37e6, 31.0e6, 97.0e6, 97.7e6};

  fprintf (ptr_File, "%20s%20s%20s\n", "p [Pa]", "s [J/kg-K]", "T [K]");
  for (int i = 0; i < 8; i++)
  {
    double rho = 0.0, T = 0.0, x = 0.0;
    R3_rho_T_x_from_p_s_ITER(p_array[i], s_array[i], rho, T, x);
    fprintf (ptr_File, "%20.8e%20.8e%20.8e\n", p_array[i], s_array[i], T);
  }

  fclose(ptr_File);
}

void Unit_Test_R5_T_ph_ITER()
{
  FILE * ptr_File;
  ptr_File = fopen("UnitTest/R5_T_ph.dat", "w");

  double h_array[8] = {4.22385130e6, 4.45997671e6, 4.95097566e6, 5.47017439e6, 6.00919444e6, 6.57122604e6, 7.15113591e6, 7.29686729e6};
  double p_array[8] = {1.0e3, 1.0e6, 5.0e6, 10.0e6, 20.0e6, 30.0e6, 40.0e6, 50.0e6};

  fprintf (ptr_File, "%20s%20s%20s\n", "p [Pa]", "h [J/kg]", "T [K]");
  for (int i = 0; i < 8; i++)
    fprintf (ptr_File, "%20.8e%20.8e%20.8e\n", p_array[i], h_array[i], R5_T_from_p_h_ITER(p_array[i], h_array[i]));

  fclose(ptr_File);
}

void Unit_Test_R5_T_ps_ITER()
{
  FILE * ptr_File;
  ptr_File = fopen("UnitTest/R5_T_ps.dat", "w");

  double s_array[8] = {1.17519976e4, 8.76987588e3, 8.40622438e3, 8.43285518e3, 8.42920714e3, 8.53640523e3, 8.67795842e3, 8.63838510e3};
  double p_array[8] = {1.0e3, 1.0e6, 5.0e6, 10.0e6, 20.0e6, 30.0e6, 40.0e6, 50.0e6};

  fprintf (ptr_File, "%20s%20s%20s\n", "p [Pa]", "s [J/kg-K]", "T [K]");
  for (int i = 0; i < 8; i++)
    fprintf (ptr_File, "%20.8e%20.8e%20.8e\n", p_array[i], s_array[i], R5_T_from_p_s_ITER(p_array[i], s_array[i]));

  fclose(ptr_File);
}

void Unit_Test_SurfTension()
{
  FILE * ptr_File;
  ptr_File = fopen("UnitTest/SurfTension.dat", "w");

  double TC_array[9] = {0.01, 75.0, 100.0, 145.0, 195.0, 225.0, 300.0, 350.0, 370.0}; // in [C]

  fprintf (ptr_File, "%20s%20s\n", "T [K]", "SurfTension [N/m]");
  for (int i = 0; i < 9; i++)
    fprintf (ptr_File, "%20.8e%20.8e\n", TC_array[i] + 273.15, surf_tension(TC_array[i] + 273.15));

  fclose(ptr_File);
}

void Unit_Test_Viscosity()
{
  FILE * ptr_File;
  ptr_File = fopen("UnitTest/Viscosity.dat", "w");

  double T_array[11] = {298.15, 298.15, 373.15, 433.15, 433.15, 873.15, 873.15, 873.15, 1173.15, 1173.15, 1173.15};
  double rho_array[11] = {998, 1200, 1000, 1, 1000, 1, 100, 600, 1, 100, 400};

  fprintf (ptr_File, "%20s%20s%20s\n", "rho [kg/m^3]", "T [K]", "Vics [Pa*s]");
  for (int i = 0; i < 11; i++)
    fprintf (ptr_File, "%20.8e%20.8e%20.9e\n", rho_array[i], T_array[i], viscosity(rho_array[i], T_array[i]));

  fclose(ptr_File);
}

void Unit_Test_ThermCond_NoEnhancement()
{
  FILE * ptr_File;
  ptr_File = fopen("UnitTest/ThermCondNoEnhancement.dat", "w");

  double T_array[4] = {298.15, 298.15, 298.15, 873.15};
  double rho_array[4] = {0.0, 998.0, 1200.0, 0.0};

  fprintf (ptr_File, "%20s%20s%20s\n", "rho [kg/m^3]", "T [K]", "lambda [W/(m-K)]");
  for (int i = 0; i < 4; i++)
    fprintf (ptr_File, "%20.8e%20.8e%20.8e\n", rho_array[i], T_array[i], thermal_conductivity_no_enhancement(rho_array[i], T_array[i]));

  fclose(ptr_File);
}

void Unit_Test_ThermCond_R1()
{
  FILE * ptr_File;
  ptr_File = fopen("UnitTest/ThermCondR1.dat", "w");

  double p_array[2] = {20.0e6, 50.0e6};
  double T_array[2] = {620.0, 620.0};

  for (int i = 0; i < 2; i++)
  {
    double rho = 1.0 / R1_specific_volume(p_array[i], T_array[i]);
    double rho_bar = rho / RHO_CRIT;
    double T_bar = T_array[i] / T_CRIT;
    double xi = correlation_length_TC(rho_bar, T_bar, zeta_R1(p_array[i], T_array[i]));
    double cp = R1_cp(p_array[i], T_array[i]);
    double cv = R1_cv(p_array[i], T_array[i]);

    fprintf (ptr_File, "%24s%20.8e\n", "p [Pa]", p_array[i]);
    fprintf (ptr_File, "%24s%20.8e\n", "T [K]",  T_array[i]);
    fprintf (ptr_File, "%24s%20.8e\n", "lambda", thermal_conductivity_R1(p_array[i], T_array[i]));
    fprintf (ptr_File, "%24s%20.8e\n", "lambda_0_bar", labmda0_bar(T_bar));
    fprintf (ptr_File, "%24s%20.8e\n", "lambda_1_bar", labmda1_bar(rho_bar, T_bar));
    fprintf (ptr_File, "%24s%20.8e\n", "rho [kg/m^3]", rho);
    fprintf (ptr_File, "%24s%20.8e\n", "drho_dp [kg/m^3/Pa]", R1_drho_dp(p_array[i], T_array[i]));
    fprintf (ptr_File, "%24s%20.8e\n", "drho_dp(R) [kg/m^3/Pa]", RHO_CRIT / P_CRIT * zeta_REF(rho_bar));
    fprintf (ptr_File, "%24s%20.8e\n", "xi [nm]", xi);
    fprintf (ptr_File, "%24s%20.8e\n", "cp [J/(kg-K)]", cp);
    fprintf (ptr_File, "%24s%20.8e\n", "cv [J/(kg-K)]", cv);
    fprintf (ptr_File, "%24s%20.8e\n", "Z(y)", Zy(xi / 0.4, rho_bar, cp / cv));
    fprintf (ptr_File, "%24s%20.8e\n\n", "mu", mu0_bar(T_bar) * mu1_bar(rho_bar, T_bar) * 1.0e-3);
  }

  fclose(ptr_File);
}

void Unit_Test_ThermCond_R2()
{
  FILE * ptr_File;
  ptr_File = fopen("UnitTest/ThermCondR2.dat", "w");

  double p_array[2] = {0.3e6, 50.0e6};
  double T_array[2] = {650.0, 800.0};

  for (int i = 0; i < 2; i++)
  {
    double rho = 1.0 / R2_specific_volume(p_array[i], T_array[i]);
    double rho_bar = rho / RHO_CRIT;
    double T_bar = T_array[i] / T_CRIT;
    double xi = correlation_length_TC(rho_bar, T_bar, zeta_R2(p_array[i], T_array[i]));
    double cp = R2_cp(p_array[i], T_array[i]);
    double cv = R2_cv(p_array[i], T_array[i]);

    fprintf (ptr_File, "%24s%20.8e\n", "p [Pa]", p_array[i]);
    fprintf (ptr_File, "%24s%20.8e\n", "T [K]",  T_array[i]);
    fprintf (ptr_File, "%24s%20.8e\n", "lambda", thermal_conductivity_R2(p_array[i], T_array[i]));
    fprintf (ptr_File, "%24s%20.8e\n", "lambda_0_bar", labmda0_bar(T_bar));
    fprintf (ptr_File, "%24s%20.8e\n", "lambda_1_bar", labmda1_bar(rho_bar, T_bar));
    fprintf (ptr_File, "%24s%20.8e\n", "rho [kg/m^3]", rho);
    fprintf (ptr_File, "%24s%20.8e\n", "drho_dp [kg/m^3/Pa]", R2_drho_dp(p_array[i], T_array[i]));
    fprintf (ptr_File, "%24s%20.8e\n", "drho_dp(R) [kg/m^3/Pa]", RHO_CRIT / P_CRIT * zeta_REF(rho_bar));
    fprintf (ptr_File, "%24s%20.8e\n", "xi [nm]", xi);
    fprintf (ptr_File, "%24s%20.8e\n", "cp [J/(kg-K)]", cp);
    fprintf (ptr_File, "%24s%20.8e\n", "cv [J/(kg-K)]", cv);
    fprintf (ptr_File, "%24s%20.8e\n", "Z(y)", Zy(xi / 0.4, rho_bar, cp / cv));
    fprintf (ptr_File, "%24s%20.8e\n\n", "mu", mu0_bar(T_bar) * mu1_bar(rho_bar, T_bar) * 1.0e-3);
  }

  fclose(ptr_File);
}

void Unit_Test_ThermCond_R3()
{
  FILE * ptr_File;
  ptr_File = fopen("UnitTest/ThermCondR3.dat", "w");

  double rho_array[2] = {222.0, 322.0};
  double T_array[2] = {647.35, 647.35};

  for (int i = 0; i < 2; i++)
  {
    double rho_bar = rho_array[i] / RHO_CRIT;
    double T_bar = T_array[i] / T_CRIT;
    double xi = correlation_length_TC(rho_bar, T_bar, zeta_R3(rho_bar, T_bar));
    double cp = R3_cp(rho_array[i], T_array[i]);
    double cv = R3_cv(rho_array[i], T_array[i]);

    fprintf (ptr_File, "%24s%20.8e\n", "rho [kg/m^3]", rho_array[i]);
    fprintf (ptr_File, "%24s%20.8e\n", "T [K]",  T_array[i]);
    fprintf (ptr_File, "%24s%20.8e\n", "lambda", thermal_conductivity_R3(rho_array[i], T_array[i]));
    fprintf (ptr_File, "%24s%20.8e\n", "lambda_0_bar", labmda0_bar(T_bar));
    fprintf (ptr_File, "%24s%20.8e\n", "lambda_1_bar", labmda1_bar(rho_bar, T_bar));
    fprintf (ptr_File, "%24s%20.8e\n", "drho_dp [kg/m^3/Pa]", RHO_CRIT / R3_dp_ddelta(rho_bar, 1.0 / T_bar));
    fprintf (ptr_File, "%24s%20.8e\n", "drho_dp(R) [kg/m^3/Pa]", RHO_CRIT / P_CRIT * zeta_REF(rho_bar));
    fprintf (ptr_File, "%24s%20.8e\n", "xi [nm]", xi);
    fprintf (ptr_File, "%24s%20.8e\n", "cp [J/(kg-K)]", cp);
    fprintf (ptr_File, "%24s%20.8e\n", "cv [J/(kg-K)]", cv);
    fprintf (ptr_File, "%24s%20.8e\n", "Z(y)", Zy(xi / 0.4, rho_bar, cp / cv));
    fprintf (ptr_File, "%24s%20.8e\n\n", "mu", mu0_bar(T_bar) * mu1_bar(rho_bar, T_bar) * 1.0e-3);
  }

  fclose(ptr_File);
}
