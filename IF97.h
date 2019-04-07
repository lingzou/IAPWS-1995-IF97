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

static const double R1_ph_Coef[20][3] = {
  {0,  0, -0.23872489924521e3   },
  {0,  1,  0.40421188637945e3   },
  {0,  2,  0.11349746881718e3   },
  {0,  6, -0.58457616048039e1   },
  {0, 22, -0.15285482413140e-3  },
  {0, 32, -0.10866707695377e-5  },
  {1,  0, -0.13391744872602e2   },
  {1,  1,  0.43211039183559e2   },
  {1,  2, -0.54010067170506e2   },
  {1,  3,  0.30535892203916e2   },
  {1,  4, -0.65964749423638e1   },
  {1, 10,  0.93965400878363e-2  },
  {1, 32,  0.11573647505340e-6  },
  {2, 10, -0.25858641282073e-4  },
  {2, 32, -0.40644363084799e-8  },
  {3, 10,  0.66456186191635e-7  },
  {3, 32,  0.80670734103027e-10 },
  {4, 32, -0.93477771213947e-12 },
  {5, 32,  0.58265442020601e-14 },
  {6, 32, -0.15020185953503e-16 }
};

double R1_T_from_p_h(double p, double h);

static const double R1_ps_Coef[20][3] = {
  {0,  0,  0.17478268058307e3   },
  {0,  1,  0.34806930892873e2   },
  {0,  2,  0.65292584978455e1   },
  {0,  3,  0.33039981775489     },
  {0, 11, -0.19281382923196e-6  },
  {0, 31, -0.24909197244573e-22 },
  {1,  0, -0.26107636489332     },
  {1,  1,  0.22592965981586     },
  {1,  2, -0.64256463395226e-1  },
  {1,  3,  0.78876289270526e-2  },
  {1, 12,  0.35672110607366e-9  },
  {1, 31,  0.17332496994895e-23 },
  {2,  0,  0.56608900654837e-3  },
  {2,  1, -0.32635483139717e-3  },
  {2,  2,  0.44778286690632e-4  },
  {2,  9, -0.51322156908507e-9  },
  {2, 31, -0.42522657042207e-25 },
  {3, 10,  0.26400441360689e-12 },
  {3, 32,  0.78124600459723e-28 },
  {4, 32, -0.30732199903668e-30 }
};

double R1_T_from_p_s(double p, double s);


const static double R2_TStar = 540.0;  // K
const static double R2_pStar = 1.0e6; // Pa

static const double R2Coef0[9][2] = {
  { 0, -0.96927686500217e1  },
  { 1,  0.10086655968018e2  },
  {-5, -0.56087911283020e-2 },
  {-4,  0.71452738081455e-1 },
  {-3, -0.40710498223928    },
  {-2,  0.14240819171444e1  },
  {-1, -0.43839511319450e1  },
  { 2, -0.28408632460772    },
  { 3,  0.21268463753307e-1 },
};

static const double R2Coefr[43][3] = {
  {1,   0, -0.17731742473213e-2   },
  {1,   1, -0.17834862292358e-1   },
  {1,   2, -0.45996013696365e-1   },
  {1,   3, -0.57581259083432e-1   },
  {1,   6, -0.50325278727930e-1   },
  {2,   1, -0.33032641670203e-4   },
  {2,   2, -0.18948987516315e-3   },
  {2,   4, -0.39392777243355e-2   },
  {2,   7, -0.43797295650573e-1   },
  {2,  36, -0.26674547914087e-4   },
  {3,   0,  0.20481737692309e-7   },
  {3,   1,  0.43870667284435e-6   },
  {3,   3, -0.32277677238570e-4   },
  {3,   6, -0.15033924542148e-2   },
  {3,  35, -0.40668253562649e-1   },
  {4,   1, -0.78847309559367e-9   },
  {4,   2,  0.12790717852285e-7   },
  {4,   3,  0.48225372718507e-6   },
  {5,   7,  0.22922076337661e-5   },
  {6,   3, -0.16714766451061e-10  },
  {6,  16, -0.21171472321355e-2   },
  {6,  35, -0.23895741934104e2    },
  {7,   0, -0.59059564324270e-17  },
  {7,  11, -0.12621808899101e-5   },
  {7,  25, -0.38946842435739e-1   },
  {8,   8,  0.11256211360459e-10  },
  {8,  36, -0.82311340897998e1    },
  {9,  13,  0.19809712802088e-7   },
  {10,  4,  0.10406965210174e-18  },
  {10, 10, -0.10234747095929e-12  },
  {10, 14, -0.10018179379511e-8   },
  {16, 29, -0.80882908646985e-10  },
  {16, 50,  0.10693031879409      },
  {18, 57, -0.33662250574171      },
  {20, 20,  0.89185845355421e-24  },
  {20, 35,  0.30629316876232e-12  },
  {20, 48, -0.42002467698208e-5   },
  {21, 21, -0.59056029685639e-25  },
  {22, 53,  0.37826947613457e-5   },
  {23, 39, -0.12768608934681e-14  },
  {24, 26,  0.73087610595061e-28  },
  {24, 40,  0.55414715350778e-16  },
  {24, 58, -0.94369707241210e-6   }
};

double R2_gamma_0(double pi, double tau);
double R2_gamma_r(double pi, double tau);

double R2_gamma_0_pi(double pi, double tau);
double R2_gamma_0_pi_pi(double pi, double tau);
double R2_gamma_0_tau(double pi, double tau);
double R2_gamma_0_tau_tau(double pi, double tau);
double R2_gamma_0_pi_tau(double pi, double tau);

double R2_gamma_r_pi(double pi, double tau);
double R2_gamma_r_pi_pi(double pi, double tau);
double R2_gamma_r_tau(double pi, double tau);
double R2_gamma_r_tau_tau(double pi, double tau);
double R2_gamma_r_pi_tau(double pi, double tau);

double R2_specific_volume(double p, double T);
double R2_specific_int_energy(double p, double T);
double R2_specific_entropy(double p, double T);
double R2_specific_enthalpy(double p, double T);
double R2_cp(double p, double T);
double R2_cv(double p, double T);
double R2_sound_speed(double p, double T);

#endif /*IF97_H*/
