#ifndef INTEGRATE_PARTICLE_H
#define INTEGRATE_PARTICLE_H
#include <Tools/cptr.h>

namespace lagrangianP {
  class momentumIntegrator: public Loci::CPTR_type {
  public:
    virtual 
      void integrateMomentum(NewParticle &p,
			     float dt,
			     const fluid_info &fluidInfo,
			     const vector3d<float> &accel, // external accelerations
			     double &Re_out) const = 0 ;
  } ;

  class CliftGauvinIntegrator : public momentumIntegrator {
    densityFunction rhop ;
  public:
  CliftGauvinIntegrator(const densityFunction &rpin) : rhop(rpin) {} 
    virtual 
      void integrateMomentum(NewParticle &p,
			     float dt,
			     const fluid_info &fluidInfo,
			     const vector3d<float> &accel, // external accelerations
			     double &Re_out) const ;
  } ;

  class CliftGauvinDropletIntegrator : public momentumIntegrator {
    densityFunction rhop ;
    tensionFunction tensionp ;
  public:
  CliftGauvinDropletIntegrator(const densityFunction &rpin,
			       const tensionFunction &tpin) : 
    rhop(rpin), tensionp(tpin) {} 
    virtual 
      void integrateMomentum(NewParticle &p,
			     float dt,
			     const fluid_info &fluidInfo,
			     const vector3d<float> &accel, // external accelerations
			     double &Re_out) const ;
  } ;


  class dropletBreakupMethod: public Loci::CPTR_type {
  public:
    virtual 
      void breakup(NewParticle &p,  float dt,
		   const fluid_info &fluidInfo) const = 0 ;
  } ;

  class noBreakupModel : public dropletBreakupMethod {
  public:
    virtual 
      void breakup(NewParticle &p,  float dt,
		   const fluid_info &fluidInfo) const ;
  } ;

  class simpleWeberBreakup : public dropletBreakupMethod {
    densityFunction rhop ;
    tensionFunction tensionp ;
    double Weber_crit ; // critical weber number
    double viscp ; // droplet viscosity
  public:
    simpleWeberBreakup(const particleBinEoS &particleBinInfo,
		       const Loci::options_list &olarg,
		       const Loci::options_list &ol,
		       const std::map<std::string,Loci::options_list> &speciesDB) ;
    virtual 
      void breakup(NewParticle &p,  float dt,
		   const fluid_info &fluidInfo) const ;
  } ;
}

#endif
