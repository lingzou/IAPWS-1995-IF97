/* References:
 * [1]. Revised Release on the IAPWS Formulation 1995 for the Thermodynamic Properties of Ordinary Water Substance for General and Scientific Use.
 *        IAPWS R6-95(2018)
 * [2]. The IAPWS Formulation 1995 for the Thermodynamic Properties of Ordinary Water Substance for General and Scientific Use,
 *        W. Wagner and A. Prub, J. Phys. Chem. Ref. Data, Vol. 31, No. 2, 2002
 */
#ifndef IAPWS1995REV_H
#define IAPWS1995REV_H

namespace IAPWS1995Rev
{
  //Reference constants, Ref. [1]
  static const double T_critical  = 647.096;    //K         //Eqn (1) page 3, Ref. [1]
  static const double p_critical  = 22.064E6;   //Pa        //Eqn (2.2b) page 398, Eqn (6.2) page 428, Ref. [2]
  static const double Rho_critical= 322.0;      //kg/m^3    //Eqn (2) page 3, Ref. [1]
  static const double R    = 0.46151805E3;      //J/kg-K    //Eqn (3) page 3, Ref. [1]
  //Note the different unit used in the code, [J/kg-K], and in the document, [kJ/kg-K].

  // Constants
  // Table 1, page 8, Ref. [1] (and Table 6.1, page 429, Ref. [2])
  static const double n0[8]    = {
              //-8.32044648201E0,   //IAPWS95, Ref. [2]
               //6.6832105268E0,    //IAPWS95, Ref. [2]
              -8.3204464837497,     //IAPWS95-Rev, Ref. [1] // The only two changes between Ref. [1] and [2].
               6.6832105275932,     //IAPWS95-Rev, Ref. [1]
               3.00632,
               0.012436,
               0.97315,
               1.27950,
               0.96956,
               0.24873
            };

  // Table 1, page 8, Ref. [1]
  static const double gamma0[5]    = {
               1.28728967,
               3.53734222,
               7.74073708,
               9.24437796,
              27.5075105
            };

  // Table 2, page 9, Ref. [1]
  static const double n1[7]    = {
               0.12533547935523E-1,   // 1
               0.78957634722828E1,    // 2
              -0.87803203303561E1,    // 3
               0.31802509345418E0,    // 4
              -0.26145533859358E0,    // 5
              -0.78199751687981E-2,   // 6
               0.88089493102134E-2    // 7
            };
  static const double n2[44]    = {
              -0.66856572307965E0,    // 8
               0.20433810950965E0,    // 9
              -0.66212605039687E-4,   // 10
              -0.19232721156002E0,    // 11
              -0.25709043003438E0,    // 12
               0.16074868486251E0,    // 13
              -0.40092828925807E-1,   // 14
                0.39343422603254E-6,  // 15
              -0.75941377088144E-5,   // 16
               0.56250979351888E-3,   // 17
              -0.15608652257135E-4,   // 18
               0.11537996422951E-8,   // 19
               0.36582165144204E-6,   // 20
              -0.13251180074668E-11,  // 21
              -0.62639586912454E-9,   // 22
              -0.10793600908932E0,    // 23
               0.17611491008752E-1,   // 24
               0.22132295167546E0,    // 25
              -0.40247669763528E0,    // 26
               0.58083399985759E0,    // 27
               0.49969146990806E-2,   // 28
              -0.31358700712549E-1,   // 29
              -0.74315929710341E0,    // 30
               0.47807329915480E0,    // 31
               0.20527940895948E-1,   // 32
              -0.13636435110343E0,    // 33
               0.14180634400617E-1,   // 34
               0.83326504880713E-2,   // 35
              -0.29052336009585E-1,   // 36
               0.38615085574206E-1,   // 37
              -0.20393486513704E-1,   // 38
              -0.16554050063734E-2,   // 39
               0.19955571979541E-2,   // 40
               0.15870308324157E-3,   // 41
              -0.16388568342530E-4,   // 42
               0.43613615723811E-1,   // 43
               0.34994005463765E-1,   // 44
              -0.76788197844621E-1,   // 45
               0.22446277332006E-1,   // 46
              -0.62689710414685E-4,   // 47
              -0.55711118565645E-9,   // 48
              -0.19905718354408E0,    // 49
               0.31777497330738E0,    // 50
              -0.11841182425981E0,    // 51
            };
  static const double n3[3]    = {
              -0.31306260323435E2,    // 52
               0.31546140237781E2,    // 53
              -0.25213154341695E4,    // 54
            };
  static const double n4[2]    = {
              -0.14874640856724E0,    // 55
               0.31806110878444E0     // 56
            };


  static const double d1[7]   = { 1,  1,  1,  2,  2,  3,  4  };  // 1 to 7
  static const double d2[44]  = {
             1,  1,  1,  2,  2,  3,  4,     // 8 to 14
             4,  5,  7,  9, 10, 11, 13,     //15 to 21
            15,  1,  2,  2,  2,  3,  4,     //22 to 28
             4,  4,  5,  6,  6,  7,  9,     //29 to 35
             9,  9,  9,  9, 10, 10, 12,     //36 to 42
             3,  4,  4,  5, 14,  3,  6,     //43 to 49
             6,  6        };                //50 to 51
  static const double d3[3]  = {  3,  3,  3 };  //52 to 54


  static const double t1[7]  = {  -0.5E0, 0.875E0, 1.E0, 0.5E0, 0.75E0, 0.375E0, 1.E0  };  //1 to 7
  static const double t2[44]  = {   4,  6, 12,  1,  5,  4,  2,    //8 to 14
            13,  9,  3,  4, 11,  4, 13,    //15 to 21
             1,  7,  1,  9, 10, 10,  3,    //22 to 28
             7, 10, 10,  6, 10, 10,  1,    //29 to 35
             2,  3,  4,  8,  6,  9,  8,    //36 to 42
            16, 22, 23, 23, 10, 50, 44,    //43 to 49
            46, 50        };  //50 to 51
  static const double t3[3]  = {    0,  1,  4      };  //52 to 54

  static const double c2[44]  = {   1,  1,  1,  1,  1,  1,  1,    //8 to 14
             1,  1,  1,  1,  1,  1,  1,    //15 to 21
             1,  2,  2,  2,  2,  2,  2,    //22 to 28
             2,  2,  2,  2,  2,  2,  2,    //29 to 35
             2,  2,  2,  2,  2,  2,  2,    //36 to 42
             3,  3,  3,  3,  4,  6,  6,    //43 to 49
             6,  6        };  //50 to 51

  static const double alpha3[3] = {  20,    20,    20     };  //52 to 54
  static const double beta3[3]  = {  150,   150,   250    };  //52 to 54
  static const double gamma3[3] = {  1.21,  1.21,  1.25   };  //52 to 54
  static const double epsil3[3] = {  1.0,   1.0,   1.0    };  //52 to 54

  static const double a4[2]     = {  3.5,  3.5      };  //55 to 56
  static const double b4[2]     = {  0.85, 0.95     };  //55 to 56
  static const double BB4[2]    = {  0.2,  0.2      };  //55 to 56
  static const double CC4[2]    = {  28.0, 32.0     };  //55 to 56
  static const double DD4[2]    = {  700,  800      };  //55 to 56
  static const double AA4[2]    = {  0.32, 0.32     };  //55 to 56
  static const double beta4[2]  = {  0.3,  0.3      };  //55 to 56
  //End of Table 2 constants

  //Equantion (5), page 3, and Table 4, page 12, Ref. [1]
  //phi_0 and its derivatives
  double func_phi_0(const double Temperature, const double rho);
  double func_d_phi0_d_delta(const double Temperature, const double rho);
  double func_d_phi0_d_delta2(const double Temperature, const double rho);
  double func_d_phi0_d_tau(const double Temperature, const double rho);
  double func_d_phi0_d_tau2(const double Temperature, const double rho);
  double func_d_phi0_d_tau_d_delta(const double Temperature, const double rho);

  void phi_0_and_its_derivatives(const double Temperature, const double rho,
            double & phi0,
            double & d_phi0_d_delta,
            double & d_phi0_d_delta2,
            double & d_phi0_d_tau,
            double & d_phi0_d_tau2,
            double & d_phi0_d_tau_d_delta);

  //Equantion (6), page 4, and Table 5, page 13, Ref. [1]
  //phi_r and its derivatives
  double func_phi_r(const double Temperature, const double rho);
  double func_phi_r(const double Temperature, const double rho);
  double func_d_phir_d_delta(const double Temperature, const double rho);
  double func_d_phir_d_delta2(const double Temperature, const double rho);
  double func_d_phir_d_tau(const double Temperature, const double rho);
  double func_d_phir_d_tau2(const double Temperature, const double rho);
  double func_d_phir_d_tau_d_delta(const double Temperature, const double rho);

  void phi_r_and_its_derivatives(const double Temperature, const double rho,
            double & phir,
            double & d_phir_d_delta,
            double & d_phir_d_delta2,
            double & d_phir_d_tau,
            double & d_phir_d_tau2,
            double & d_phir_d_tau_d_delta);


  //Thermodynamics properties
  double Pressure(const double Temperature, const double rho);
  double Entropy(const double Temperature, const double rho);
  double IntEnergy(const double Temperature, const double rho);
  double Enthalpy(const double Temperature, const double rho);
  double GibbsFreeEnergy(const double Temperature, const double rho);
  double Cv(const double Temperature, const double rho);
  double Cp(const double Temperature, const double rho);
  double SoundSpeed(const double Temperature, const double rho);
  double JouleThomson_Coeff(const double Temperature, const double rho);
  double IsoThrottling_Coeff(const double Temperature, const double rho);
  double Isentropic_TemperaturePressure_Coeff(const double Temperature, const double rho);
  double Second_Virial_Coeff(const double Temperature, const double rho);
  double Third_Virial_Coeff(const double Temperature, const double rho);

  // Aux equations for saturation line and properties
  // See discussions at the end of section 2.3.1 of Ref. [2],
  //   the difference between the results of these aux equations and the above equations are extremely small,
  //   they are not thermodynamically consistent.

  // constants, page 399 and 400 of Ref. [2]
  static const double a[6]  = { -7.85951783,  1.84408259,  -11.7866497,  22.6807411,  -15.9618719,  1.80122502  };
  static const double coeff_P[6]  = { 1.0,    1.5,    3.0,    3.5,    4.0,    7.5    };
  static const double b[6]  = { 1.99274064,   1.09965342,  -0.510839303,  -1.75493479,  -45.5170352,  -6.74694450E5  };
  static const double coeff_L[6]  = { 1.0/3.0,    2.0/3.0,  5.0/3.0,  16.0/3.0,  43.0/3.0,  110.0/3.0  };
  static const double c[6]  = { -2.0315024,    -2.6830294,  -5.38626492,  -17.2991605,  -44.7586581,  -63.9201063  };
  static const double coeff_V[6]  = { 2.0/6.0,    4.0/6.0,  8.0/6.0,  18.0/6.0,  37.0/6.0,  71.0/6.0  };

  static const double alpha_0   = 1.0e3;  // [J/kg]
  static const double psi_0     = alpha_0 / T_critical; // [J/kg-K]
  static const double d_alpha   = -1135.905627715;
  static const double d_psi     = 2319.5246;
  static const double d[5]      = { -5.65134998E-8,  2690.66631,  127.287297,  -135.003439,  0.981825814  };
  static const double coeff_alpha[5]  = { -19.0,    1.0,    4.5,    5.0,    54.5    };

  void temperatureInRange(const double T);  // helper function
  void pressureInRange(const double p);     // helper function
  double psat_from_T(const double T);       // eqn (2.5), page 398, Ref. [2]
  double Tsat_from_p(const double p);       // helper function, iterative inverse from psat_from_T
  double dpsat_dT(const double T);          // eqn (2.5a), page 399, Ref. [2]
  double rho_l_sat_from_T(const double T);  // eqn (2.6), page 399, Ref. [2]
  double rho_g_sat_from_T(const double T);  // eqn (2.7), page 399, Ref. [2]
  double alpha_ratio(const double T);       // eqn (2.9), page 400, Ref. [2]
  double psi_ratio(const double T);         // eqn (2.9a), page 400, Ref. [2]
  double e_l_sat_from_T(const double T);    // eqn (2.12), page 400, Ref. [2]
  double e_g_sat_from_T(const double T);    // eqn (2.13), page 400, Ref. [2]
  double h_l_sat_from_T(const double T);    // eqn (2.10), page 400, Ref. [2]
  double h_g_sat_from_T(const double T);    // eqn (2.11), page 400, Ref. [2]
  double s_l_sat_from_T(const double T);    // eqn (2.14), page 401, Ref. [2]
  double s_g_sat_from_T(const double T);    // eqn (2.15), page 401, Ref. [2]
}

#endif //IAPWS1995REV_H
