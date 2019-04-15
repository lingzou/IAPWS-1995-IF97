#include <cstdio>
#include <sstream>
#include <iomanip>
#include <iostream>

#include "unit_tests.h"
#include "IF97_helper.h"

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

  out_ss << "p_sat [Pa] = " << std::scientific << std::setprecision(8) << p_sat_from_T(300.0) << std::endl;
  out_ss << "p_sat [Pa] = " << std::scientific << std::setprecision(8) << p_sat_from_T(500.0) << std::endl;
  out_ss << "p_sat [Pa] = " << std::scientific << std::setprecision(8) << p_sat_from_T(600.0) << std::endl << std::endl;

  out_ss << "T_sat [K] = " << std::scientific << std::setprecision(8) << T_sat_from_p(0.1e6) << std::endl;
  out_ss << "T_sat [K] = " << std::scientific << std::setprecision(8) << T_sat_from_p(1.0e6) << std::endl;
  out_ss << "T_sat [K] = " << std::scientific << std::setprecision(8) << T_sat_from_p(10.0e6) << std::endl << std::endl;

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
