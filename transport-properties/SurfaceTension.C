#include <math.h>
#include "SurfaceTension.h"

double surf_tension(double T)
{
  double tau = 1.0 - T / 647.096;

  return 235.8e-3 * pow(tau, 1.256) * (1.0 - 0.625 * tau);
}
