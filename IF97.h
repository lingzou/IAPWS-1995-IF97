#ifndef IF97_H
#define IF97_H
#include <stdlib.h>
#include <stdio.h>
/*
 *  Reference [1]
 *
 * "Revised Release on the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam",
 *    IAPWS R7-97, The International Association for the Properties of Water and Steam (IAPWS), August, 2007.
 */

extern "C"
{
// Eqn. (1)-(4), constants on Fig. 1 of Ref. [1], page 4-5
const static double R_GAS     = 0.461526e3;          // J/kg-K
const static double T_CRIT    = 647.096;             // K
const static double P_CRIT    = 22.064e6;            // Pa
const static double RHO_CRIT  = 322.0;               // kg/m³
const static double IF97_T_MIN     = 273.15;              // K
const static double IF97_T_13      = 623.15;              // K
const static double IF97_T_25      = 1073.15;             // K
const static double IF97_T_MAX     = 2273.15;             // K
const static double IF97_SAT_P_MIN = 611.213;             // Pa; See page 35, Ref. [1]
const static double IF97_P_MAX     = 100.0e6;             // Pa

/***************************************************************
 * Region 1
 ***************************************************************/
// Table 1, Ref. [1], page 6
static const double B23_n[5] = {
    0.34805185628969e3,
   -0.11671859879975e1,
    0.10192970039326e-2,
    0.57254459862746e3,
    0.13918839778870e2
};

// Eqn. (5) and (6), Refs. [1], page 5-6
double B23_p_from_T(double T);
double B23_T_from_p(double p);

// Table 2, Ref. [1], page 7
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

// Region 1 reference values, page 6, Ref. [1]
const static double R1_TStar = 1386.0;  // K
const static double R1_pStar = 16.53e6; // Pa

// Table 4, Ref. [1], page 8
double R1_gamma(double pi, double tau);
double R1_gamma_pi(double pi, double tau);
double R1_gamma_tau(double pi, double tau);
double R1_gamma_pi_pi(double pi, double tau);
double R1_gamma_tau_tau(double pi, double tau);
double R1_gamma_pi_tau(double pi, double tau);

// Table 3, Ref. [1], page 8
double R1_specific_volume(double p, double T);
double R1_specific_int_energy(double p, double T);
double R1_specific_entropy(double p, double T);
double R1_specific_enthalpy(double p, double T);
double R1_cp(double p, double T);
double R1_cv(double p, double T);
double R1_sound_speed(double p, double T);

// Table 6, Ref. [1], page 10
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

// Eqn. (11), Ref. [1], page 10
double R1_T_from_p_h(double p, double h);

// Table 8, Ref. [1], page 11
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

// Eqn. (13), Ref. [1], page 11
double R1_T_from_p_s(double p, double s);

/***************************************************************
 * Region 2
 ***************************************************************/
// Region 2 reference values, page 13
const static double R2_TStar = 540.0;  // K
const static double R2_pStar = 1.0e6; // Pa

// Table 10, Ref. [1], page 13
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

// Table 11, Ref. [1], page 14
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

// Eqn. (16)-(17), Ref. [1], page 13
double R2_gamma_0(double pi, double tau);
double R2_gamma_r(double pi, double tau);

// Table 13, Ref. [1], page 16
double R2_gamma_0_pi(double pi, double tau);
double R2_gamma_0_pi_pi(double pi, double tau);
double R2_gamma_0_tau(double pi, double tau);
double R2_gamma_0_tau_tau(double pi, double tau);
double R2_gamma_0_pi_tau(double pi, double tau);

// Table 14, Ref. [1], page 16
double R2_gamma_r_pi(double pi, double tau);
double R2_gamma_r_pi_pi(double pi, double tau);
double R2_gamma_r_tau(double pi, double tau);
double R2_gamma_r_tau_tau(double pi, double tau);
double R2_gamma_r_pi_tau(double pi, double tau);

// Table 12, Ref. [1], page 15
double R2_specific_volume(double p, double T);
double R2_specific_int_energy(double p, double T);
double R2_specific_entropy(double p, double T);
double R2_specific_enthalpy(double p, double T);
double R2_cp(double p, double T);
double R2_cv(double p, double T);
double R2_sound_speed(double p, double T);

/***************************************************************
 * Region 2 (Metastable)
 ***************************************************************/
// Table 10, Ref. [1], page 13 (with addjusted n^o_1 and n^o_2)
static const double R2MetaCoef0[9][2] = {
  { 0, -0.96937268393049e1  },
  { 1,  0.10087275970006e2  },
  {-5, -0.56087911283020e-2 },
  {-4,  0.71452738081455e-1 },
  {-3, -0.40710498223928    },
  {-2,  0.14240819171444e1  },
  {-1, -0.43839511319450e1  },
  { 2, -0.28408632460772    },
  { 3,  0.21268463753307e-1 },
};

// Table 16, Ref. [1], page 18
static const double R2MetaCoefr[13][3] = {
  {1,   0, -0.73362260186506e-2   },
  {1,   2, -0.88223831943146e-1   },
  {1,   5, -0.72334555213245e-1   },
  {1,  11, -0.40813178534455e-2   },
  {2,   1,  0.20097803380207e-2   },
  {2,   7, -0.53045921898642e-1   },
  {2,  16, -0.76190409086970e-2   },
  {3,   4, -0.63498037657313e-2   },
  {3,  16, -0.86043093028588e-1   },
  {4,   7,  0.75321581522770e-2   },
  {4,  10, -0.79238375446139e-2   },
  {5,   9, -0.22888160778447e-3   },
  {5,  10, -0.26456501482810e-2   }
};

// Eqn. (16) with adjusted n^o_1 and n^o_2, Ref. [1], page 13
double R2Meta_gamma_0(double pi, double tau);
// Eqn. (19), Ref. [1], page 18
double R2Meta_gamma_r(double pi, double tau);

// Table 13, Ref. [1], page 16
double R2Meta_gamma_0_pi(double pi, double tau);
double R2Meta_gamma_0_pi_pi(double pi, double tau);
double R2Meta_gamma_0_tau(double pi, double tau);
double R2Meta_gamma_0_tau_tau(double pi, double tau);
double R2Meta_gamma_0_pi_tau(double pi, double tau);

// Table 17, Ref. [1], page 19
double R2Meta_gamma_r_pi(double pi, double tau);
double R2Meta_gamma_r_pi_pi(double pi, double tau);
double R2Meta_gamma_r_tau(double pi, double tau);
double R2Meta_gamma_r_tau_tau(double pi, double tau);
double R2Meta_gamma_r_pi_tau(double pi, double tau);

// Table 12, Ref. [1], page 15
double R2Meta_specific_volume(double p, double T);
double R2Meta_specific_int_energy(double p, double T);
double R2Meta_specific_entropy(double p, double T);
double R2Meta_specific_enthalpy(double p, double T);
double R2Meta_cp(double p, double T);
double R2Meta_cv(double p, double T);
double R2Meta_sound_speed(double p, double T);

// Table 19, Ref. [1], page 22
static const double B2bc_n[5] = {
  0.90584278514723e3,
 -0.67955786399241e0,
  0.12809002730136e-3,
  0.26526571908428e4,
  0.45257578905948e1
};

// Eqn. (20)-(21), Ref. [1], page 21
double B2bc_p_from_h(double h);
double B2bc_h_from_p(double p);

// Table 20, Ref. [1], page 22
static const double R2a_ph_Coef[34][3] = {
  {0,  0,  0.10898952318288e4   },
  {0,  1,  0.84951654495535e3   },
  {0,  2, -0.10781748091826e3   },
  {0,  3,  0.33153654801263e2   },
  {0,  7, -0.74232016790248e1   },
  {0, 20,  0.11765048724356e2   },
  {1,  0,  0.18445749355790e1   },
  {1,  1, -0.41792700549624e1   },
  {1,  2,  0.62478196935812e1   },
  {1,  3, -0.17344563108114e2   },
  {1,  7, -0.20058176862096e3   },
  {1,  9,  0.27196065473796e3   },
  {1, 11, -0.45511318285818e3   },
  {1, 18,  0.30919688604755e4   },
  {1, 44,  0.25226640357872e6   },
  {2,  0, -0.61707422868339e-2  },
  {2,  2, -0.31078046629583     },
  {2,  7,  0.11670873077107e2   },
  {2, 36,  0.12812798404046e9   },
  {2, 38, -0.98554909623276e9   },
  {2, 40,  0.28224546973002e10  },
  {2, 42, -0.35948971410703e10  },
  {2, 44,  0.17227349913197e10  },
  {3, 24, -0.13551334240775e5   },
  {3, 44,  0.12848734664650e8   },
  {4, 12,  0.13865724283226e1   },
  {4, 32,  0.23598832556514e6   },
  {4, 44, -0.13105236545054e8   },
  {5, 32,  0.73999835474766e4   },
  {5, 36, -0.55196697030060e6   },
  {5, 42,  0.37154085996233e7   },
  {6, 34,  0.19127729239660e5   },
  {6, 44, -0.41535164835634e6   },
  {7, 28, -0.62459855192507e2   }
};

// Eqn. (22), Ref. [1], page 22
double R2a_T_from_p_h(double p, double h);

// Table 21, Ref. [1], page 23
static const double R2b_ph_Coef[38][3] = {
  {0,  0,  0.14895041079516e4   },
  {0,  1,  0.74307798314034e3   },
  {0,  2, -0.97708318797837e2   },
  {0, 12,  0.24742464705674e1   },
  {0, 18, -0.63281320016026     },
  {0, 24,  0.11385952129658e1   },
  {0, 28, -0.47811863648625     },
  {0, 40,  0.85208123431544e-2  },
  {1,  0,  0.93747147377932     },
  {1,  2,  0.33593118604916e1   },
  {1,  6,  0.33809355601454e1   },
  {1, 12,  0.16844539671904     },
  {1, 18,  0.73875745236695     },
  {1, 24, -0.47128737436186     },
  {1, 28,  0.15020273139707     },
  {1, 40, -0.21764114219750e-2  },
  {2,  2, -0.21810755324761e-1  },
  {2,  8, -0.10829784403677     },
  {2, 18, -0.46333324635812e-1  },
  {2, 40,  0.71280351959551e-4  },
  {3,  1,  0.11032831789999e-3  },
  {3,  2,  0.18955248387902e-3  },
  {3, 12,  0.30891541160537e-2  },
  {3, 24,  0.13555504554949e-2  },
  {4,  2,  0.28640237477456e-6  },
  {4, 12, -0.10779857357512e-4  },
  {4, 18, -0.76462712454814e-4  },
  {4, 24,  0.14052392818316e-4  },
  {4, 28, -0.31083814331434e-4  },
  {4, 40, -0.10302738212103e-5  },
  {5, 18,  0.28217281635040e-6  },
  {5, 24,  0.12704902271945e-5  },
  {5, 40,  0.73803353468292e-7  },
  {6, 28, -0.11030139238909e-7  },
  {7,  2, -0.81456365207833e-13 },
  {7, 28, -0.25180545682962e-10 },
  {9,  1, -0.17565233969407e-17 },
  {9, 40,  0.86934156344163e-14 }
};

// Eqn. (23), Ref. [1], page 23
double R2b_T_from_p_h(double p, double h);

// Table 22, Ref. [1], page 24
static const double R2c_ph_Coef[23][3] = {
  {-7,  0, -0.32368398555242e13 },
  {-7,  4,  0.73263350902181e13 },
  {-6,  0,  0.35825089945447e12 },
  {-6,  2, -0.58340131851590e12 },
  {-5,  0, -0.10783068217470e11 },
  {-5,  2,  0.20825544563171e11 },
  {-2,  0,  0.61074783564516e6  },
  {-2,  1,  0.85977722535580e6  },
  {-1,  0, -0.25745723604170e5  },
  {-1,  2,  0.31081088422714e5  },
  { 0,  0,  0.12082315865936e4  },
  { 0,  1,  0.48219755109255e3  },
  { 1,  4,  0.37966001272486e1  },
  { 1,  8, -0.10842984880077e2  },
  { 2,  4, -0.45364172676660e-1 },
  { 6,  0,  0.14559115658698e-12},
  { 6,  1,  0.11261597407230e-11},
  { 6,  4, -0.17804982240686e-10},
  { 6, 10,  0.12324579690832e-6 },
  { 6, 12, -0.11606921130984e-5 },
  { 6, 16,  0.27846367088554e-4 },
  { 6, 20, -0.59270038474176e-3 },
  { 6, 22,  0.12918582991878e-2 }
};

// Eqn. (24), Ref. [1], page 23
double R2c_T_from_p_h(double p, double h);

// Table 25, Ref. [1], page 26
static const double R2a_ps_Coef[46][3] = {
  {-1.50, -24, -0.39235983861984e6  },
  {-1.50, -23,  0.51526573827270e6  },
  {-1.50, -19,  0.40482443161048e5  },
  {-1.50, -13, -0.32193790923902e3  },
  {-1.50, -11,  0.96961424218694e2  },
  {-1.50, -10, -0.22867846371773e2  },
  {-1.25, -19, -0.44942914124357e6  },
  {-1.25, -15, -0.50118336020166e4  },
  {-1.25,  -6,  0.35684463560015    },
  {-1.00, -26,  0.44235335848190e5  },
  {-1.00, -21, -0.13673388811708e5  },
  {-1.00, -17,  0.42163260207864e6  },
  {-1.00, -16,  0.22516925837475e5  },
  {-1.00,  -9,  0.47442144865646e3  },
  {-1.00,  -8, -0.14931130797647e3  },
  {-0.75, -15, -0.19781126320452e6  },
  {-0.75, -14, -0.23554399470760e5  },
  {-0.50, -26, -0.19070616302076e5  },
  {-0.50, -13,  0.55375669883164e5  },
  {-0.50,  -9,  0.38293691437363e4  },
  {-0.50,  -7, -0.60391860580567e3  },
  {-0.25, -27,  0.19363102620331e4  },
  {-0.25, -25,  0.42660643698610e4  },
  {-0.25, -11, -0.59780638872718e4  },
  {-0.25,  -6, -0.70401463926862e3  },
  { 0.25,   1,  0.33836784107553e3  },
  { 0.25,   4,  0.20862786635187e2  },
  { 0.25,   8,  0.33834172656196e-1 },
  { 0.25,  11, -0.43124428414893e-4 },
  { 0.50,   0,  0.16653791356412e3  },
  { 0.50,   1, -0.13986292055898e3  },
  { 0.50,   5, -0.78849547999872    },
  { 0.50,   6,  0.72132411753872e-1 },
  { 0.50,  10, -0.59754839398283e-2 },
  { 0.50,  14, -0.12141358953904e-4 },
  { 0.50,  16,  0.23227096733871e-6 },
  { 0.75,   0, -0.10538463566194e2  },
  { 0.75,   4,  0.20718925496502e1  },
  { 0.75,   9, -0.72193155260427e-1 },
  { 0.75,  17,  0.20749887081120e-6 },
  { 1.00,   7, -0.18340657911379e-1 },
  { 1.00,  18,  0.29036272348696e-6 },
  { 1.25,   3,  0.21037527893619    },
  { 1.25,  15,  0.25681239729999e-3 },
  { 1.50,   5, -0.12799002933781e-1 },
  { 1.50,  18, -0.82198102652018e-5 }
};

// Eqn. (25), Ref. [1], page 25
double R2a_T_from_p_s(double p, double s);

// Table 26, Ref. [1], page 27
static const double R2b_ps_Coef[44][3] = {
  {-6,  0,  0.31687665083497e6  },
  {-6, 11,  0.20864175881858e2  },
  {-5,  0, -0.39859399803599e6  },
  {-5, 11, -0.21816058518877e2  },
  {-4,  0,  0.22369785194242e6  },
  {-4,  1, -0.27841703445817e4  },
  {-4, 11,  0.99207436071480e1  },
  {-3,  0, -0.75197512299157e5  },
  {-3,  1,  0.29708605951158e4  },
  {-3, 11, -0.34406878548526e1  },
  {-3, 12,  0.38815564249115    },
  {-2,  0,  0.17511295085750e5  },
  {-2,  1, -0.14237112854449e4  },
  {-2,  6,  0.10943803364167e1  },
  {-2, 10,  0.89971619308495    },
  {-1,  0, -0.33759740098958e4  },
  {-1,  1,  0.47162885818355e3  },
  {-1,  5, -0.19188241993679e1  },
  {-1,  8,  0.41078580492196    },
  {-1,  9, -0.33465378172097    },
  { 0,  0,  0.13870034777505e4  },
  { 0,  1, -0.40663326195838e3  },
  { 0,  2,  0.41727347159610e2  },
  { 0,  4,  0.21932549434532e1  },
  { 0,  5, -0.10320050009077e1  },
  { 0,  6,  0.35882943516703    },
  { 0,  9,  0.52511453726066e-2 },
  { 1,  0,  0.12838916450705e2  },
  { 1,  1, -0.28642437219381e1  },
  { 1,  2,  0.56912683664855    },
  { 1,  3, -0.99962954584931e-1 },
  { 1,  7, -0.32632037778459e-2 },
  { 1,  8,  0.23320922576723e-3 },
  { 2,  0, -0.15334809857450    },
  { 2,  1,  0.29072288239902e-1 },
  { 2,  5,  0.37534702741167e-3 },
  { 3,  0,  0.17296691702411e-2 },
  { 3,  1, -0.38556050844504e-3 },
  { 3,  3, -0.35017712292608e-4 },
  { 4,  0, -0.14566393631492e-4 },
  { 4,  1,  0.56420857267269e-5 },
  { 5,  0,  0.41286150074605e-7 },
  { 5,  1, -0.20684671118824e-7 },
  { 5,  2,  0.16409393674725e-8 }
};

// Eqn. (26), Ref. [1], page 26
double R2b_T_from_p_s(double p, double s);

// Table 27, Ref. [1], page 28
static const double R2c_ps_Coef[30][3] = {
  {-2, 0,  0.90968501005365e3   },
  {-2, 1,  0.24045667088420e4   },
  {-1, 0, -0.59162326387130e3   },
  { 0, 0,  0.54145404128074e3   },
  { 0, 1, -0.27098308411192e3   },
  { 0, 2,  0.97976525097926e3   },
  { 0, 3, -0.46966772959435e3   },
  { 1, 0,  0.14399274604723e2   },
  { 1, 1, -0.19104204230429e2   },
  { 1, 3,  0.53299167111971e1   },
  { 1, 4, -0.21252975375934e2   },
  { 2, 0, -0.31147334413760     },
  { 2, 1,  0.60334840894623     },
  { 2, 2, -0.42764839702509e-1  },
  { 3, 0,  0.58185597255259e-2  },
  { 3, 1, -0.14597008284753e-1  },
  { 3, 5,  0.56631175631027e-2  },
  { 4, 0, -0.76155864584577e-4  },
  { 4, 1,  0.22440342919332e-3  },
  { 4, 4, -0.12561095013413e-4  },
  { 5, 0,  0.63323132660934e-6  },
  { 5, 1, -0.20541989675375e-5  },
  { 5, 2,  0.36405370390082e-7  },
  { 6, 0, -0.29759897789215e-8  },
  { 6, 1,  0.10136618529763e-7  },
  { 7, 0,  0.59925719692351e-11 },
  { 7, 1, -0.20677870105164e-10 },
  { 7, 3, -0.20874278181886e-10 },
  { 7, 4,  0.10162166825089e-9  },
  { 7, 5, -0.16429828281347e-9  }
};

// Eqn. (27), Ref. [1], page 27
double R2c_T_from_p_s(double p, double s);

/***************************************************************
 * Region 3
 ***************************************************************/
// Table 30, Ref. [1], page 30
static const double R3Coef[40][3] = {
  {0,   0,    0.10658070028513e1  },
  {0,   0,   -0.15732845290239e2  },
  {0,   1,    0.20944396974307e2  },
  {0,   2,   -0.76867707878716e1  },
  {0,   7,    0.26185947787954e1  },
  {0,  10,   -0.28080781148620e1  },
  {0,  12,    0.12053369696517e1  },
  {0,  23,   -0.84566812812502e-2 },
  {1,   2,   -0.12654315477714e1  },
  {1,   6,   -0.11524407806681e1  },
  {1,  15,    0.88521043984318    },
  {1,  17,   -0.64207765181607    },
  {2,   0,    0.38493460186671    },
  {2,   2,   -0.85214708824206    },
  {2,   6,    0.48972281541877e1  },
  {2,   7,   -0.30502617256965e1  },
  {2,  22,    0.39420536879154e-1 },
  {2,  26,    0.12558408424308    },
  {3,   0,   -0.27999329698710    },
  {3,   2,    0.13899799569460e1  },
  {3,   4,   -0.20189915023570e1  },
  {3,  16,   -0.82147637173963e-2 },
  {3,  26,   -0.47596035734923    },
  {4,   0,    0.43984074473500e-1 },
  {4,   2,   -0.44476435428739    },
  {4,   4,    0.90572070719733    },
  {4,  26,    0.70522450087967    },
  {5,   1,    0.10770512626332    },
  {5,   3,   -0.32913623258954    },
  {5,  26,   -0.50871062041158    },
  {6,   0,   -0.22175400873096e-1 },
  {6,   2,    0.94260751665092e-1 },
  {6,  26,    0.16436278447961    },
  {7,   2,   -0.13503372241348e-1 },
  {8,  26,   -0.14834345352472e-1 },
  {9,   2,    0.57922953628084e-3 },
  {9,  26,    0.32308904703711e-2 },
  {10,  0,    0.80964802996215e-4 },
  {10,  1,   -0.16557679795037e-3 },
  {11, 26,   -0.44923899061815e-4 },
};

// Table 32, Ref. [1], page 32
double R3_phi(double delta, double tau);
double R3_phi_delta(double delta, double tau);
double R3_phi_delta_delta(double delta, double tau);
double R3_phi_tau(double delta, double tau);
double R3_phi_tau_tau(double delta, double tau);
double R3_phi_delta_tau(double delta, double tau);

// Table 31, Ref. [1], page 31
double R3_p(double rho, double T);
double R3_specific_int_energy(double rho, double T);
double R3_specific_entropy(double rho, double T);
double R3_specific_enthalpy(double rho, double T);
double R3_cv(double rho, double T);
double R3_cp(double rho, double T);
double R3_sound_speed(double rho, double T);

/***************************************************************
 * Region 4
 ***************************************************************/
// Table 34, Ref. [1], page 34
static const double R4Coef[10] = {
  0.11670521452767e4,
 -0.72421316703206e6,
 -0.17073846940092e2,
  0.12020824702470e5,
 -0.32325550322333e7,
  0.14915108613530e2,
 -0.48232657361591e4,
  0.40511340542057e6,
 -0.23855557567849,
  0.65017534844798e3
};

// Eqn. (30), page 33; Eqn. (31), page 35, Ref. [1]
double R4_p_sat_from_T(double T);
double R4_T_sat_from_p(double p);
void checkTSatValid(double T);
void checkPSatValid(double p);

/***************************************************************
 * Region 5
 ***************************************************************/
// Table 37, Ref. [1], page 37
static const double R5Coef0[6][2] = {
  { 0, -0.13179983674201e2  },
  { 1,  0.68540841634434e1  },
  {-3, -0.24805148933466e-1 },
  {-2,  0.36901534980333    },
  {-1, -0.31161318213925e1  },
  { 2, -0.32961626538917    }
};

// Table 38, Ref. [1], page 37
static const double R5Coefr[6][3] = {
  {1, 1,  0.15736404855259e-2},
  {1, 2,  0.90153761673944e-3},
  {1, 3, -0.50270077677648e-2},
  {2, 3,  0.22440037409485e-5},
  {2, 9, -0.41163275453471e-5},
  {3, 7,  0.37919454822955e-7}
};

// Eqn. (33), page 36; Eqn. (34), page 37, Ref. [1]
double R5_gamma_0(double pi, double tau);
double R5_gamma_r(double pi, double tau);

// Table 40, Ref. [1], page 39
double R5_gamma_0_pi(double pi, double tau);
double R5_gamma_0_pi_pi(double pi, double tau);
double R5_gamma_0_tau(double pi, double tau);
double R5_gamma_0_tau_tau(double pi, double tau);
double R5_gamma_0_pi_tau(double pi, double tau);

// Table 41, Ref. [1], page 39
double R5_gamma_r_pi(double pi, double tau);
double R5_gamma_r_pi_pi(double pi, double tau);
double R5_gamma_r_tau(double pi, double tau);
double R5_gamma_r_tau_tau(double pi, double tau);
double R5_gamma_r_pi_tau(double pi, double tau);

// Table 39, Ref. [1], page 38
double R5_specific_volume(double p, double T);
double R5_specific_int_energy(double p, double T);
double R5_specific_entropy(double p, double T);
double R5_specific_enthalpy(double p, double T);
double R5_cp(double p, double T);
double R5_cv(double p, double T);
double R5_sound_speed(double p, double T);
}
#endif /*IF97_H*/
