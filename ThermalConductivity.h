#ifndef THERMAL_CONDUCTIVITY_H
#define THERMAL_CONDUCTIVITY_H
#include <math.h>

/*
 * Reference [4]
 *
 * "Release on the IAPWS Formulation 2011 for the Thermal Conductivity of Ordinary Water Substance",
 *   IAPWS R15-11, The International Association for the Properties of Water and Steam (IAPWS), September, 2011.
 */

extern "C"
{
static const double PI = acos(-1);

// Table 1, Ref. [4], page 6
static const double L0[5] = {
   2.443221e-3,
   1.323095e-2,
   6.770357e-3,
  -3.454586e-3,
   4.096266e-4
};

// Table 2, Ref. [4], page 6
static const double L1[5][6] = {
  {  1.60397357, -0.646013523,  0.111443906,  0.102997357 , -0.0504123634,  0.00609859258},
  {  2.33771842, -2.78843778 ,  1.53616167 , -0.463045512 ,  0.0832827019, -0.00719201245},
  {  2.19650529, -4.54580785 ,  3.55777244 , -1.40944978  ,  0.275418278 , -0.0205938816 },
  { -1.21051378,  1.60812989 , -0.621178141,  0.0716373224,  0.0         ,  0.0          },
  { -2.7203370 ,  4.57586331 , -3.18369245 ,  1.1168348   , -0.19268305  ,  0.012913842  }
};

// Table 6, Ref. [4], page 11
static const double A[6][5] = {
  {  6.53786807199516,  6.52717759281799,   5.35500529896124,   1.55225959906681 ,   1.11999926419994  },
  { -5.61149954923348, -6.30816983387575,  -3.96415689925446,   0.464621290821181,   0.595748562571649 },
  {  3.39624167361325,  8.08379285492595,   8.91990208918795,   8.93237374861479 ,   9.88952565078920  },
  { -2.27492629730878, -9.82240510197603, -12.0338729505790 , -11.0321960061126  , -10.3255051147040   },
  { 10.2631854662709 , 12.1358413791395 ,   9.19494865194302,   6.16780999933360 ,   4.66861294457414  },
  {  1.97815050331519, -5.54349664571295,  -2.16866274479712,  -0.965458722086812,  -0.503243546373828 }
};

// Eqn. (16), page 5; Eqn. (17), page 6; Eqn. (18), page 7; Ref. [4]
double labmda0_bar(double T_bar);
double labmda1_bar(double rho_bar, double T_bar);
double labmda2_bar(double rho_bar, double T_bar, double cp, double cv, double mu_bar);

// Eqn. (24), page 7, Ref. [4]
double zeta_R1(double p, double T);
double zeta_R2(double p, double T);
double zeta_R3(double rho_bar, double T_bar);

// Eqn. (25), page 11, Ref. [4]
double zeta_REF(double rho_bar);

// Eqn. (22)-(21), page 7, Ref. [4]
double correlation_length_TC(double rho_bar, double T_bar, double zeta);
double Zy(double y, double rho_bar, double kappa);

// Eqn. (16) with lambda_2=1, page 5, Ref. [4]
double thermal_conductivity_no_enhancement(double rho, double T);

// Eqn. (16) with lambda_2, page 5, Ref. [4]
double thermal_conductivity_R1(double p, double T);
double thermal_conductivity_R2(double p, double T);
double thermal_conductivity_R3(double rho, double T);
double thermal_conductivity_R5(double p, double T);
}
#endif /*THERMAL_CONDUCTIVITY_H*/
