#ifndef IF97_H
#define IF97_H

const static double Rgas    = 0.461526;            // J/kg-K
const static double Tcrit   = 647.096;             // K
const static double Pcrit   = 22.064e6;            // Pa
const static double Rhocrit = 322.0;               // kg/m³
const static double Tmin    = 273.15;              // K
const static double Tmax    = 1073.15;             // K
const static double Pmin    = 0.000611213e6;       // Pa
const static double Pmax    = 100.0e6;             // Pa
const static double MW      = 0.018015268;         // kg/mol

static double B23_n[] = {
    0.34805185628969e3,
   -0.11671859879975e1,
    0.10192970039326e-2,
    0.57254459862746e3,
    0.13918839778870e2
};

void a_func();
double B23_p_from_T(double T);
double B23_T_from_p(double p);

#endif /*IF97_H*/
