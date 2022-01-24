#ifndef SURFACE_TENSION_H
#define SURFACE_TENSION_H

/*
 * Reference [2]
 *
 * "Revised Release on Surface Tension of Ordinary Water Substance",
 *   IAPWS R1-76, The International Association for the Properties of Water and Steam (IAPWS), June, 2014.
 */

extern "C"
{
double surf_tension(double T);
}
#endif /*SURFACE_TENSION_H*/
