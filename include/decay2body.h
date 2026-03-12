#include "fizika.h"
//#include "coordinate.h"

#ifndef __DECAY_2_BODY__
#define __DECAY_2_BODY__

class decay2body{

  static constexpr double Mp = 0.93827;
  fizika::vector3 ex,ey,ez;
 public:
  decay2body(){
    ex.setXYZ(1.0,0.0,0.0);
    ey.setXYZ(0.0,1.0,0.0);
    ez.setXYZ(0.0,0.0,1.0);
  }
  virtual ~decay2body(){};

  /*******************************************************************************
   * FNCTION:
  *******************************************************************************/
  std::vector<fizika::lorentz4>  decay( fizika::lorentz4 &ref, fizika::lorentz4 &parent, double m1, double m2, double costheta, double phi){
    
    double M = parent.m(); double M2 = M*M;
    double term = (M2 - (m1 + m2)*(m1 + m2)) * (M2 - (m1 - m2)*(m1 - m2));
    double p = (term > 0 ? std::sqrt(term)/(2.0*M) : 0.0);

    fizika::vector3 uy = ref.vect().cross(parent.vect()).unit();//parent.vect().cross(ref.vect()).unit();
    fizika::vector3 uz = parent.vect().unit();
    fizika::vector3 ux = uy.cross(uz);

    fizika::transformer::basis E{{1,0,0}, {0,1,0}, {0,0,1}};
    fizika::transformer::basis U{
      {ux.x(),ux.y(),ux.z()}, {uy.x(),uy.y(),uy.z()}, {uz.x(),uz.y(),uz.z()}
    };
    fizika::transformer trans(E,U);
    
    double sintheta = sqrt(1-costheta*costheta);
    double px = p*sintheta*cos(phi);
    double py = p*sintheta*sin(phi);
    double pz = p*costheta;

    fizika::vector3 vU(px,py,pz); fizika::vector3 vE = trans.toE(vU);
    fizika::vector3 vp1(vE.x(),vE.y(),vE.z());

    std::vector<fizika::lorentz4> list;
    
    fizika::lorentz4 vl1(  vE.x(), vE.y(), vE.z(), sqrt(vE.mag2()+m1*m1));
    fizika::lorentz4 vl2( -vE.x(),-vE.y(),-vE.z(), sqrt(vE.mag2()+m2*m2));
    
    fizika::vector3 boost = parent.boostVector();
    vl1.boost( boost );  vl2.boost( boost );
    list.push_back(vl1); list.push_back(vl2);
    return list;
  }
  /*******************************************************************************
   * FNCTION:
  *******************************************************************************/
  fizika::lorentz4  eprime(double beamE, double Q2, double xB){
    double nu = Q2 / (2.0 * Mp * xB);

      double enu = nu;
      // 2) scattered‐electron energy E' = E_beam − ν
      double Eprime = beamE - nu;
      //printf("nu = %f, eprime = %f\n",nu,Eprime);
      double y  = nu/beamE;
      //double gam     = gamma_bj(xB,Q2);
      //double eps   = epsilon(y,gam);

      double sin2t2 = Q2/(4.0*beamE*Eprime);
      double sint2  = sqrt(sin2t2);
      double t2     = asin(sint2);
      //printf("theta = %f\n",(t2*2.0)*57.29);
      double cosTh  = cos(t2*2.0);
      // 3) scattering angle θ_e from Q2 = 2 E E' (1 - cos θ)
      //double cosTh = 1.0 - Q2 / (2.0 * beamE * Eprime);
      //printf("cosTh = %f, eprime = %f\n",cosTh, Eprime);
      // guard against small numerical slip:
      if (cosTh > +1) cosTh = +1;
      if (cosTh < -1) cosTh = -1;
      double sinTh = std::sqrt(1.0 - cosTh*cosTh);
      //printf("cosTh = %f\n",cosTh);
      // 4) build four‐vector (px, py, pz, E)
      //    assume scattering in the x–z plane: py = 0
      double px = Eprime * sinTh;
      double pz = Eprime * cosTh;    
      return fizika::lorentz4(px, 0.0, pz, sqrt(Eprime*Eprime+0.0005*0.0005));
  }

  double getXb(fizika::lorentz4 &e, fizika::lorentz4 &ep){
    double nu = e.e()-ep.e();
    fizika::lorentz4 q2 = e-ep;
    return -q2.m2()/(2.0*Mp*nu);
  }
  /*******************************************************************************
   * FNCTION:
  *******************************************************************************/
  fizika::lorentz4 toFrame(fizika::lorentz4 &ref, fizika::lorentz4 &parent, fizika::lorentz4 &p1){

    fizika::vector3 ux = parent.vect().cross(ref.vect()).unit();
    fizika::vector3 uz = parent.vect().unit();
    fizika::vector3 uy = uz.cross(ux);
        
    fizika::lorentz4     d1 = p1;
    d1.boost(-parent.boostVector());

    //printf("before = %9.5f, after %9.5f\n",p1.p(),d1.p());

    fizika::transformer::basis E{{1,0,0}, {0,1,0}, {0,0,1}};
    fizika::transformer::basis U{
      {ux.x(),ux.y(),ux.z()}, {uy.x(),uy.y(),uy.z()}, {uz.x(),uz.y(),uz.z()}
    };
    //CoordinateTransformer trans(E, U);
    fizika::transformer trans(E,U);
    
    fizika::vector3 vE{d1.px(),d1.py(),d1.pz()};
    fizika::vector3 vU = trans.toU(vE);

    return fizika::lorentz4(vU.x(),vU.y(),vU.z(),d1.e());
  }
  
};

#endif
