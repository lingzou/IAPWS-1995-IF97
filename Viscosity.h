#ifndef VISCOSITY_H
#define VISCOSITY_H

#include "IF97_helper.h"

static const double MU0_H[4] = {1.67752, 2.20462, 0.6366564, -0.241605};
double mu0_bar(double T_bar);

static const int MU1_POS[21][2] = {
    {0, 0},
    {1, 0},
    {2, 0},
    {3, 0},
    {0, 1},
    {1, 1},
    {2, 1},
    {3, 1},
    {5, 1},
    {0, 2},
    {1, 2},
    {2, 2},
    {3, 2},
    {4, 2},
    {0, 3},
    {1, 3},
    {0, 4},
    {3, 4},
    {4, 5},
    {3, 6},
    {5, 6}
};

static const double MU1_H[21] = {
   5.20094e-1,
   8.50895e-2,
  -1.08374   ,
  -2.89555e-1,
   2.22531e-1,
   9.99115e-1,
   1.88797   ,
   1.26613   ,
   1.20573e-1,
  -2.81378e-1,
  -9.06851e-1,
  -7.72479e-1,
  -4.89837e-1,
  -2.57040e-1,
   1.61913e-1,
   2.57399e-1,
  -3.25372e-2,
   6.98452e-2,
   8.72102e-3,
  -4.35673e-3,
  -5.93264e-4
};
double mu1_bar(double rho_bar, double T_bar); // Follow IF97 function arguments sequence

double viscosity(double rho, double T, bool critical_enhancement = false); // Follow IF97 function arguments sequence
double zeta(double rho_bar, double T_bar);
double correlation_length(double rho_bar, double T_bar); /*in nm */
double mu2_bar(double rho_bar, double T_bar);

double mu2_bar(double rho_bar, double T_bar, double xi);

#endif /*VISCOSITY_H*/