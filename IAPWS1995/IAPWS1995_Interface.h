#ifndef IAPWS1995REV_INTERFACE_H
#define IAPWS1995REV_INTERFACE_H

namespace IAPWS1995Rev
{
  /***************************************************************
   * (p, T)-based properties
   ***************************************************************/
  double rho_l_from_pT(double p, double T);
  double e_l_from_pT(double p, double T);
  double h_l_from_pT(double p, double T);
  double s_l_from_pT(double p, double T);
  double cv_l_from_pT(double p, double T);
  double cp_l_from_pT(double p, double T);
  double c_l_from_pT(double p, double T);
  double k_l_from_pT(double p, double T);
  double mu_l_from_pT(double p, double T);

  double rho_g_from_pT(double p, double T);
  double e_g_from_pT(double p, double T);
  double h_g_from_pT(double p, double T);
  double s_g_from_pT(double p, double T);
  double cv_g_from_pT(double p, double T);
  double cp_g_from_pT(double p, double T);
  double c_g_from_pT(double p, double T);
  double k_g_from_pT(double p, double T);
  double mu_g_from_pT(double p, double T);
}
#endif //IAPWS1995REV_INTERFACE_H
