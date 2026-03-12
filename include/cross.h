#include "w_kernels.hpp"

#ifndef __CROSS_SECTION__
#define __CROSS_SECTION__



struct densityMatrix {
  Wkernels::Mat4 u, l, s;
  double dsigmaT_dt, dsigmaL_dt;   // [nb / GeV²]
  double Pl = 0.0, SL = 0.0, ST = 0.0;   // polarisations
};

class cross {
 protected:

  static constexpr double Mp = 0.93827;
  static constexpr double alem  = 1.0 / 137.035999084;
  double beamEnergy;

  TRandom3 rnd;
  
 public:

  cross(){beamEnergy = 10.6;}
  cross(double be){beamEnergy = be;}
  
  ~cross(){}
  
  double dsigma_3fold(double xB, double Q2, double y,
		      double eps, double dsigmaT_dt, double dsigmaL_dt)
  {
    double pref = alem/(2.0*M_PI);
    double kin  = (y*y)/(1.0-eps)*(1.0-xB)/xB/Q2;
    return pref * kin * (dsigmaT_dt + eps*dsigmaL_dt);
  }

  double dsigma_7fold(double WUU,double WLU,double WUL,
		      double WLL,double WUT,double WLT,
		      double Pl,double SL,double ST,
		      double sigma3)
  {
    double S = WUU + Pl*WLU + SL*WUL + Pl*SL*WLL
      + ST*WUT + Pl*ST*WLT;
    //printf(" S = %f\n",S);
    return sigma3 * S / (4.0*M_PI*M_PI);
  }


  double weight2(double q2, double xb, double mt, double pol){
    std::vector<double> v;
    for(int i = 0; i < 2000; i++){
      double fi_d = rnd.Uniform(-M_PI, M_PI);
      double fi_p = rnd.Uniform(-M_PI, M_PI);
      double t_p = acos(rnd.Uniform(-1,1));
      double t_d = acos(rnd.Uniform(-1,1));
      double w = weight(q2,xb,mt,t_p,fi_p,t_d,fi_d,pol);
      v.push_back(w);
    }

    double summ = 0.0;
    for(int i = 0; i < (int) v.size(); i++) summ += v[i];
    return summ/((int) v.size());
  }
  
  double weight(double q2, double xb, double mt,
		double theta_p, double phi_p, double theta_d, double phi_d, double polar){

    densityMatrix dm = loadMatrix(xb,q2,mt);
    double y = Y(q2,xb);
    double eps = Eps(q2,xb);
    double sigma3 = dsigma_3fold(xb,q2,y,eps,dm.dsigmaT_dt,dm.dsigmaL_dt);
    auto WUU = Wkernels::UU(dm.u,eps, phi_p, phi_d);
    auto WLU = Wkernels::LU(dm.u,eps, phi_p, phi_d);

    double cosTheta = std::cos(theta_d);
    double sinTheta = std::sin(theta_d);

    double W_LU = cosTheta*cosTheta*WLU.LL+std::sqrt(2)*cosTheta*sinTheta*WLU.LT+sinTheta*sinTheta*WLU.TT;
    double W_UU = cosTheta*cosTheta*WUU.LL+std::sqrt(2)*cosTheta*sinTheta*WUU.LT+sinTheta*sinTheta*WUU.TT;

    double pol = polar<0?-1.0:1.0;

    double w  = 1000*dsigma_7fold(W_UU,W_LU,0,0,0,0,
			     pol,dm.SL,dm.ST, sigma3);
    //printf("ST = %f\n",dm.ST);
    //printf("t = %f %f %f %f %f --- %f %f (%f %f) sigma 3 = %f , w = %f \n",mt,q2,xb, WUU.LL, WLU.LL, W_UU,W_LU, cosTheta, sinTheta,sigma3,w);
    return w;
    
  }

  double getW_LU(double q2, double xb, double mt,
		double theta_p, double phi_p, double theta_d, double phi_d, double polar){
    densityMatrix dm = loadMatrix(xb,q2,mt);
    double eps = Eps(q2,xb);
    auto WLU = Wkernels::LU(dm.u,eps, phi_p, phi_d);
    double cosTheta = std::cos(theta_d);
    double sinTheta = std::sin(theta_d);
    //printf("%f %f\n",WLU.LT, WLU.LL);
    double W_LU = cosTheta*cosTheta*WLU.LL+std::sqrt(2)*cosTheta*sinTheta*WLU.LT+sinTheta*sinTheta*WLU.TT;
    return W_LU;
  }
  
  double Y(double q2,double xb){
    double nu = q2/(2.0*Mp*xb);
    return nu/beamEnergy;
  }
  double gamma_bj(double xB, double Q2)
  { return 2.0 * xB * Mp / std::sqrt(Q2); }
  //-------------------------------------------------------------------------
  double epsilon(double y, double γ)
  {
    double y2γ2 = y*y*γ*γ; return (1.0-y-0.25*y2γ2)/(1.0-y+0.5*y*y+0.25*y2γ2);
  }
  
  double Eps(double q2, double xb){
    double gg = gamma_bj(xb,q2);
    double y = Y(q2,xb);
    return epsilon(y,gg);
  }
  
  densityMatrix loadMatrix(double xB, double Q2,
  				 double t /*GeV²*/)
  {
    densityMatrix ph{};
    
    /* ---------------------------------------------------------------
     *  TODO:
     *  1.  Set the u,l,s helicity matrices for the chosen (xB,Q²,t)
     *      or read them from a file.
     *  2.  Provide dσT/dt and dσL/dt (nb/GeV²).
     * --------------------------------------------------------------*/
    ph.dsigmaT_dt = 1.0;   // placeholder
    ph.dsigmaL_dt = 1.0;   // placeholder
    
    ph.u[Wkernels::h('0')][Wkernels::h('0')][Wkernels::h('+')][Wkernels::h('+')] = {1.00*exp(0.4*t), 0.0000};   // real example
    ph.u[Wkernels::h('0')][Wkernels::h('0')][Wkernels::h('0')][Wkernels::h('0')] = {0.00, 0.0000};   // real example
    ph.u[Wkernels::h('0')][Wkernels::h('0')][Wkernels::h('0')][Wkernels::h('+')] = {0.00, 0.2000};   // real example
    ph.u[Wkernels::h('0')][Wkernels::h('0')][Wkernels::h('-')][Wkernels::h('+')] = {0.00, 0.0000};
    
    return ph;
  }
  
};
#endif
