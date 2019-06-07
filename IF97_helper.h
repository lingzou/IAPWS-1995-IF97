#ifndef IF97_HELPER_H
#define IF97_HELPER_H

#include "IF97.h"
#include "IF97_interpolation.h"

struct IF97_Interpolation;

const static double R3_T_list[51] =
{
  624.15, 625.15, 626.15, 627.15, 628.15,
  629.15, 630.15, 631.15, 632.15, 633.15,
  634.15, 635.15, 636.15, 637.15, 638.15,
  639.15, 640.15, 641.15, 642.15, 643.15,
  644.15, 645.15, 646.15,
  646.2, 646.25, 646.3, 646.35, 646.4,
  646.45, 646.5, 646.55, 646.6, 646.65,
  646.7, 646.75, 646.8, 646.85, 646.9,
  646.95, 647, 647.05,
  647.06, 647.07, 647.08, 647.09,
  647.091, 647.092, 647.093, 647.094, 647.095, 647.0955
};

// R3_rho_l_sat_guess and R3_rho_g_sat_guess are really not guessed, they are from
// NIST webbook: https://webbook.nist.gov/chemistry/fluid/
const static double R3_rho_l_sat_guess[51] =
{
  570.64, 566.46, 562.15, 557.72, 553.14, 548.41, 543.5, 538.41, 533.11, 527.59, 521.82, 515.79, 509.45, 502.78, 495.74, 488.27, 480.29,
  471.67, 462.18, 451.43, 438.64, 422.26, 398.68, 397.17, 395.62, 394.01, 392.34, 390.62, 388.82, 386.95, 384.99, 382.92, 380.74, 378.41,
  375.9, 373.16, 370.12, 366.66, 362.56, 357.34, 349.5, 347.18, 344.3, 340.39, 333.9, 332.95, 331.88, 330.65, 329.16, 327.16, 325.7
};

const static double R3_rho_g_sat_guess[51] =
{
  116.1, 118.68, 121.37, 124.17, 127.09, 130.14, 133.33, 136.67, 140.19, 143.9, 147.82, 151.99, 156.43, 161.2, 166.35, 171.95, 178.11,
  184.98, 192.77, 201.84, 212.79, 226.84, 247.22, 248.56, 249.94, 251.38, 252.88, 254.44, 256.08, 257.8, 259.62, 261.54, 263.6, 265.81,
  268.21, 270.86, 273.82, 277.22, 281.29, 286.51, 294.37, 296.68, 299.56, 303.46, 309.84, 310.83, 311.93, 313.2, 314.73, 316.79, 318.27
};


int findRegion(double p, double T);
double rho_l_sat_from_T(double T);
double rho_g_sat_from_T(double T);

void genR3_sat_line();
void genR4_sat_line();
double R3_rho_from_p_T_ITER(double p, double T);
void R3_T_x_from_p_h_ITER(double p, double h, double &T, double &x);
void R3_T_x_from_p_s_ITER(double p, double s, double &T, double &x);
double R3_dp_ddelta(double delta, double tau);

double R1_drho_dp(double p, double T);
double R2_drho_dp(double p, double T);

double R5_T_from_p_h_ITER(double p, double h);
double R5_T_from_p_s_ITER(double p, double s);

#endif /*IF97_HELPER_H*/
