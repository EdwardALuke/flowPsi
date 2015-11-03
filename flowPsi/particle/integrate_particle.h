//#############################################################################
//#
//# Copyright 2015, Mississippi State University
//#
//# This file is part of the flowPsi computational fluid dynamics solver.
//#
//# The flowPsi solver is free software: you can redistribute it and/or modify
//# it under the terms of the GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The flowPsi solver is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# GNU General Public License for more details.
//#
//# You should have received a copy of the GNU General Public License
//# along with the flowPsi solver.  If not, see <http://www.gnu.org/licenses>
//#
//#############################################################################

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
