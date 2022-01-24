#ifndef VISCOSITY_H
#define VISCOSITY_H

#include "IF97_helper.h"

/*
 * Reference [3]
 *
 * "Release on the IAPWS Formulation 2008 for the Viscosity of Ordinary Water Substance",
 *   IAPWS R12-08, The International Association for the Properties of Water and Steam (IAPWS), September, 2008.
 */

extern "C"
{
// Table 1, page 5, Ref. [3]
static const double MU0_H[4] = {1.67752, 2.20462, 0.6366564, -0.241605};

// Eqn. (11), page 5, Ref. [3]
double mu0_bar(double T_bar);

// Table 2 (part of it), page 5, Ref. [3]
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

// Table 2 (value part), page 5, Ref. [3]
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

// Eqn. (12), page 5, Ref. [3]
double mu1_bar(double rho_bar, double T_bar); // Follow IF97 function arguments sequence

// Eqn. (10) without mu_2 critical enhancement part, which is not required by IF standard, page 5, Ref. [3]
double viscosity(double rho, double T); // Follow IF97 function arguments sequence (rho, T), instead of (T, rho)
}
#endif /*VISCOSITY_H*/
