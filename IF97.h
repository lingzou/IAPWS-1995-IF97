#ifndef IF97_H
#define IF97_H

const static double Rgas    = 0.461526e3;          // J/kg-K
const static double Tcrit   = 647.096;             // K
const static double Pcrit   = 22.064e6;            // Pa
const static double Rhocrit = 322.0;               // kg/m³
const static double Tmin    = 273.15;              // K
const static double Tmax    = 1073.15;             // K
const static double Pmin    = 0.000611213e6;       // Pa
const static double Pmax    = 100.0e6;             // Pa
const static double MW      = 0.018015268;         // kg/mol

static const double B23_n[] = {
    0.34805185628969e3,
   -0.11671859879975e1,
    0.10192970039326e-2,
    0.57254459862746e3,
    0.13918839778870e2
};

void a_func();
double B23_p_from_T(double T);
double B23_T_from_p(double p);

static const double R1Coef[34][3] = {
  {0,   -2,  0.14632971213167     },
  {0,   -1, -0.84548187169114     },
  {0,    0, -0.37563603672040e1   },
  {0,    1,  0.33855169168385e1   },
  {0,    2, -0.95791963387872     },
  {0,    3,  0.15772038513228     },
  {0,    4, -0.16616417199501e-1  },
  {0,    5,  0.81214629983568e-3  },
  {1,   -9,  0.28319080123804e-3  },
  {1,   -7, -0.60706301565874e-3  },
  {1,   -1, -0.18990068218419e-1  },
  {1,    0, -0.32529748770505e-1  },
  {1,    1, -0.21841717175414e-1  },
  {1,    3, -0.52838357969930e-4  },
  {2,   -3, -0.47184321073267e-3  },
  {2,    0, -0.30001780793026e-3  },
  {2,    1,  0.47661393906987e-4  },
  {2,    3, -0.44141845330846e-5  },
  {2,   17, -0.72694996297594e-15 },
  {3,   -4, -0.31679644845054e-4  },
  {3,    0, -0.28270797985312e-5  },
  {3,    6, -0.85205128120103e-9  },
  {4,   -5, -0.22425281908000e-5  },
  {4,   -2, -0.65171222895601e-6  },
  {4,   10, -0.14341729937924e-12 },
  {5,   -8, -0.40516996860117e-6  },
  {8,  -11, -0.12734301741641e-8  },
  {8,   -6, -0.17424871230634e-9  },
  {21, -29, -0.68762131295531e-18 },
  {23, -31,  0.14478307828521e-19 },
  {29, -38,  0.26335781662795e-22 },
  {30, -39, -0.11947622640071e-22 },
  {31, -40,  0.18228094581404e-23 },
  {32, -41, -0.93537087292458e-25 }
};

const static double R1_TStar = 1386.0;  // K
const static double R1_pStar = 16.53e6; // Pa

double R1_gamma(double p, double T);
double R1_gamma_pi(double p, double T);
double R1_gamma_tau(double p, double T);
double R1_gamma_pi_pi(double p, double T);
double R1_gamma_tau_tau(double p, double T);
double R1_gamma_pi_tau(double p, double T);

double R1_specific_volume(double p, double T);
double R1_specific_int_energy(double p, double T);
double R1_specific_entropy(double p, double T);
double R1_specific_enthalpy(double p, double T);
double R1_cp(double p, double T);
double R1_cv(double p, double T);
double R1_sound_speed(double p, double T);

#endif /*IF97_H*/
