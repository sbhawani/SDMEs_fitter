
#include "fizika.h"
#include "decay2body.h"

#ifndef __REACTION__
#define __REACTION__

namespace sim {
  
  enum class kinematics {
    Q2,XB,NU,EPS,Y
  };
  
  struct particle {
    int pid; int status;
    fizika::lorentz4 v;fizika::vector3  vrtx;
    particle(int p, int s,fizika::lorentz4 __v) : pid(p), status(s), v(__v), vrtx(0.0,0.0,0.0){} 
  };
    
  class event {
    
  public:
    std::vector<particle> pts;
    std::vector<double>   params;
    double beamPol = 1;
    double targetPol = 0;
    int beamPID = 11;
    int targetPID = 2212;
    double beamEnergy = 10.6;
    
  public:
    event(){}
    virtual ~event(){}
    void reset(){params.clear(); pts.clear();}
    void add(particle p){pts.push_back(p);};
    void add(double p){params.push_back(p);}
    fizika::lorentz4 vector(int row){ return pts[row].v;};
    
    void show(int row){
      printf("%3d  0 %3d %6d %3d  0  %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f\n",
	     row+1,
	     pts[row].status,pts[row].pid, pts[row].status,
	     pts[row].v.px(),pts[row].v.py(),pts[row].v.pz(),
	     pts[row].v.e(),pts[row].v.m(),
	     pts[row].vrtx.x(),pts[row].vrtx.y(),pts[row].vrtx.z()
	     );
    }

    bool hasNaN(){
      for(int i = 0; i < (int) pts.size(); i++){
	if(std::isnan(pts[i].v.px())) return true;
	if(std::isnan(pts[i].v.py())) return true;
	if(std::isnan(pts[i].v.pz())) return true;
      }
      return false;
    }
    
    void show(){
      printf("%5d %4d %4d %5.2f %5.2f %5d %8.5f %5d %3d %e ",
	     (int) pts.size(), 1, 0, beamPol, targetPol,
	     beamPID, beamEnergy, targetPID, 1, 0.6
	     );
      for(int i = 0; i < (int) params.size(); i++) printf("%9.5f ",params[i]);
      printf("\n");
      for(int i = 0; i < (int) pts.size(); i++) show(i);
    }

    void rotateZ(double phi){
      for(int i = 0; i < (int) pts.size(); i++){
	pts[i].v.rotateZ(phi);
      }
    }
  };
  
  class reaction {
  public:
    static constexpr double Mp = 0.93827;

    double Mv = 1.02;
  
    double rq2, rxb, rbeam, rpol;
  
    fizika::lorentz4 e_in,e_out;
    fizika::lorentz4 p_in,p_out;
    fizika::lorentz4 q,vec;
    fizika::lorentz4 decayOne, decayTwo;

    double eps,  y,  enu;
    double decayTheta, decayPhi;
    double prodTheta, prodPhi;
    double decayDaughterMassOne,decayDaughterMassTwo;
    int    particleIDproduce = 113;
    int    particleIDdecayOne = 211;
    int    particleIDdecayTwo = -211;
    
    TRandom3 rng;

    decay2body decay;
  
  public:

    reaction(){}
    reaction(double beame_, double massv_){
      rbeam = beame_; Mv = massv_; 
    }
    //-------------------------------------------------------------------------
    void set(double beame_, double massv_){ rbeam = beame_; Mv = massv_;}
    //-------------------------------------------------------------------------
    void setDecay(double dm1, double dm2){
      decayDaughterMassOne=dm1; decayDaughterMassTwo = dm2;
    }

    void setDecayIds(int parent, int d1, int d2){
      particleIDproduce = parent; particleIDdecayOne = d1; particleIDdecayTwo = d2;
    }
    //-------------------------------------------------------------------------
    void generate(double q2_, double xb_){
      rq2 = q2_; rxb = xb_; produce__(rbeam,rq2,rxb);
      produce__2();produce__3();
      rpol = rng.Uniform(-1, 1);
      //printf("POLARIZATION = %f\n",rpol);
    }
    //-------------------------------------------------------------------------
    void produce__(double beame_, double q2_, double xb_){
      rng.SetSeed(0);
      p_in.setXYZE(0,0,0, Mp);
      e_in.setXYZE(0,0, beame_, sqrt(beame_*beame_+0.0005*0.0005));
      e_out = decay.eprime(beame_,q2_,xb_);//scatteredElectron(beame_,q2_,xb_);
      q = e_in-e_out;
    }
    //-------------------------------------------------------------------------
    void produce__2(){
      double cos_t = rng.Uniform(-1,1);double phi = rng.Uniform(-M_PI, M_PI);
      //phi = 0.0; //--- debug
      prodTheta   = acos(cos_t); prodPhi = phi; fizika::lorentz4 cm = q + p_in;
      std::vector<fizika::lorentz4> products = decay.decay(e_out,cm,Mp,Mv,cos_t,phi);
      p_out = products[0]; vec = products[1];
    }
    //-------------------------------------------------------------------------
    void produce__3(){
      double cos_t = rng.Uniform(-1,1); double phi = rng.Uniform(-M_PI, M_PI);
      //cos_t = 0.5; phi = 0.0; // -- debug
      decayTheta   = acos(cos_t); decayPhi     = phi;
      std::vector<fizika::lorentz4> products = decay.decay(q,vec,decayDaughterMassOne,
							   decayDaughterMassTwo,cos_t,phi);
      decayOne = products[0]; decayTwo = products[1];
    }
    //-------------------------------------------------------------------------
    double get(kinematics type){
      switch (type) {
      case kinematics::Q2:
	return rq2;
      case kinematics::XB:
	return rxb;
      default:
	return -1;
      }
    }
    double Q2(){return rq2;}
    double E(){return rbeam;}
    double xB(){return rxb;}
    double pol(){return rpol;}
    
    double Y(){
      double nu = rq2/(2.0*Mp*rxb);
      return nu/rbeam;
    }

    double Eps(){
      double gg = gamma_bj(rxb,rq2);
      double y = Y();
      return epsilon(y,gg);
    }
    //-------------------------------------------------------------------------
    void show(){
      printf("header: production [%9.5f, %9.5f], decay [%9.5f %9.5f]\n",prodTheta,
	     prodPhi,decayTheta,decayPhi);
      e_in.print("e  in"); e_out.print("e out"); p_in.print("p  in");
      p_out.print("p out"); vec.print("v out");
    }
    //-------------------------------------------------------------------------
    void getEvent(event &ev){
      ev.beamEnergy = rbeam;
      ev.reset(); ev.add(rq2); ev.add(rxb); ev.add(prodTheta); ev.add(prodPhi);
      ev.add(decayTheta); ev.add(decayPhi); ev.add(particle(  11,0,e_in));
      ev.add(particle(2212,0,p_in));        ev.add(particle(  11,1,e_out));
      ev.add(particle(2212,1,p_out));       ev.add(particle( particleIDproduce,2,vec));
      
      ev.add(particle(particleIDdecayTwo,1,decayOne));    ev.add(particle( particleIDdecayOne,1,decayTwo));
      ev.beamPol = rpol<0?-1.0:1.0;
      double phi = rng.Uniform(-M_PI, M_PI);
      ev.rotateZ(phi);
    }
    //-------------------------------------------------------------------------
    double gamma_bj(double xB, double Q2)
    { return 2.0 * xB * Mp / std::sqrt(Q2); }
    //-------------------------------------------------------------------------
    double epsilon(double y, double γ)
    {
      double y2γ2 = y*y*γ*γ; return (1.0-y-0.25*y2γ2)/(1.0-y+0.5*y*y+0.25*y2γ2);
    }
    //-------------------------------------------------------------------------
  };

  class candidate {
  public:
    reaction react;
    double  weight;
    candidate(double beam, double massv, double m1, double m2){
      react.set(beam,massv);weight=0.0; react.setDecay(m1,m2);
    }
  
  };
}
#endif
