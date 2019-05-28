//#############################################################################
//#
//# Copyright 2015-2019, Mississippi State University
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

// This header is now the main place to define (and also add new)
// particle types. Different particle types may be defined
// to represent different kind of particles. See the comments in
// "particle_base.h" for the details of the current design.
#ifndef PARTICLE_H
#define PARTICLE_H

#include <Loci.h>
#include <iostream>
#include <list>
#include <string>
#include "particle_base.h"
#include "particle_model.h"
// here we define a concrete particle type
namespace lagrangianP {
  int check_list_size(const options_list &o, std::string option) ;
  bool check_list_units(const options_list &o, std::string option, std::string unit) ;  
  enum particle_info {NEW_PARTICLE=1,POOL_PARTICLE=2} ;
  
  struct auxiliary_info {
    int turbophoresisIndex ;
    int pfluidIndex ;
    int RlIndex ;     // temporal damping for particle timestep
    int TsigmaIndex ; // sigma of turbulent fluctuations
  } ;
    
  struct fluid_info {
    vector3d<float> cellcenter ;
    float rfluid ;
    vector3d<float> ufluid,umax,umin ;
    tensor3d<float> ugrad ;
    float viscfluid ;
    float fvolpp ; // fluid volume per particle
  } ;


  struct NewParticle {
    int id ;
    Entity cell ;
    short bin ;  // particle property bin
    short info ; // basic info
    float number_in_parcel ; // Number of physical particles this particle
    // represents
    float mass ;
    Loci::Array<vector3d<float>,2> pos ; // position history
    Loci::Array<vector3d<float>,3> vel ; // velocity history
    vector3d<float> ufluidp ; // Turbulent fluctuation model 

    double get_density(const densityFunction &rhop) const {
      return rhop.get_density() ;
    }

    double get_diameter(const densityFunction &rhop) const {
      double density = rhop.get_density() ;
      return max(pow(6.0*mass/(pi*density),1./3.),1e-20) ;
    }
    void update_mass(double diameter, const particleEoSList &pm) {
      double rho = pm[bin].rhop.get_density() ;
      double m = pi*pow(diameter,3)*rho/6.0 ;
      
      mass = m ;
    }
    NewParticle(double diameter, vect3d p, vect3d pv,
                double parcel, int b, const particleEoSList &pm) {
      id = 0 ;
      cell = -1 ;
      bin = b ;
      info = 1 ; // first bit means newly created
      mass = 1.0 ;
      pos[0] = vector3d<float>(p.x,p.y,p.z) ;
      pos[1] = pos[0] ;
      vel[0] = vector3d<float>(pv.x,pv.y,pv.z) ;
      vel[1] = vel[0] ;
      vel[2] = vel[0] ;
      ufluidp = vector3d<float>(0,0,0) ;
      number_in_parcel = parcel ;
      update_mass(diameter,pm) ;
    }
    
    NewParticle(Entity e=-1, double m=0.0,
                vec3d p=vec3d(0,0,0), vec3d pv=vec3d(0,0,0),
                double parcel=1.0,int b=0) :
      id(0),cell(e),bin(b),
      number_in_parcel(parcel) {
      info = 1 ; // first bit means newly created
      mass = m ;
      pos[0] = vector3d<float>(p.x,p.y,p.z) ;
      pos[1] = pos[0] ;
      vel[0] = vector3d<float>(pv.x,pv.y,pv.z) ;
      vel[1] = vel[0] ;
      vel[2] = vel[0] ;
      ufluidp = vector3d<float>(0,0,0) ;
    }
    
    int sizeState() { return 18;}
    // Now we need to provide the interface required in
    // the base particle definition
    int
    get_id() const {return id ;}
    void
    set_id(int i) { id = i ;}

    Entity
    get_cell() const {return cell ;}
    void
    set_cell(Entity c) { cell = c ;}

    vec3d get_position() const {
      return vec3d(pos[0].x,pos[0].y,pos[0].z) ;
    }
    void
    set_position(const vec3d& p)
    { pos[0] = vector3d<float>(p.x,p.y,p.z) ;
      pos[1] = pos[0]; }

    vec3d get_velocity() const {
      return vec3d(vel[0].x,vel[0].y,vel[0].z) ;
    }
    void set_velocity(const vec3d& v)
    { vel[0] = vector3d<float> (v.x,v.y,v.z) ;
      vel[1] = vel[0] ;}

    void pack(unsigned char* buffer, size_t& position) const {
      NewParticle* b =
        reinterpret_cast<NewParticle*>(buffer+position) ;
      *b = *this ;
      position += pack_size() ;
    }

    void unpack(unsigned char* buffer, size_t& position) {
      NewParticle* b =
        reinterpret_cast<NewParticle*>(buffer+position) ;
      *this = *b ;
      position += pack_size() ;
    }

    size_t
    pack_size() const {
      return sizeof(NewParticle) ;
    }
    
    double get_mass() const {
      return mass ;
    }

    
  } ; // end of NewParticle
  
  // I/O for BasicParticle
  inline std::ostream&
  operator<<(std::ostream& s, const NewParticle& p) {
    s << "particle info    : " << std::endl ;
    s << "         id      : " << p.id << std::endl ;
    s << "         cell    : " << p.cell << std::endl ;
    s << "         mass    : " << p.mass << std::endl ;
    s << "         position: " << p.pos[0] << std::endl ;
    s << "         velocity: " << p.vel[0] << std::endl ;
    s << "      parcel size: " << p.number_in_parcel << std::endl ;
    s << "              bin: " << p.bin << std::endl ;
    return s ;
  }
  // not intended for actual use, only defined to pass the compilation
  inline std::istream&
  operator>>(std::istream& s, NewParticle& p) {
    s >> p.id >> p.cell
      >> p.pos[0] >> p.pos[1] 
      >> p.vel[0] >> p.vel[1] >> p.vel[2]
      >> p.mass
      >> p.number_in_parcel >> p.bin ;
    return s ;
  }
  
  // alias the NewParticle
#define NEW_PARTICLE_FORMULATION
#ifdef NEW_PARTICLE_FORMULATION
  typedef NewParticle Particle ;
#endif
  
  template <class T> struct summation_list {
    void operator()(T &l1, const T &l2) {
      std::copy(l2.begin(),l2.end(),std::back_inserter(l1)) ;
    }
  } ;


} // end of namespace lagrangianP

namespace Loci {
  // Loci traits stuff for NewParticle
  template<>
  struct data_schema_traits<lagrangianP::NewParticle> {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      CompoundDatatypeP cmpd =
        CompoundFactory(lagrangianP::NewParticle()) ;
      LOCI_INSERT_TYPE(cmpd, lagrangianP::NewParticle, id) ;
      LOCI_INSERT_TYPE(cmpd, lagrangianP::NewParticle, cell) ;
      LOCI_INSERT_TYPE(cmpd, lagrangianP::NewParticle, bin) ;
      LOCI_INSERT_TYPE(cmpd, lagrangianP::NewParticle, info) ;
      LOCI_INSERT_TYPE(cmpd, lagrangianP::NewParticle, number_in_parcel) ;
      LOCI_INSERT_TYPE(cmpd, lagrangianP::NewParticle, pos) ;
      LOCI_INSERT_TYPE(cmpd, lagrangianP::NewParticle, vel) ;
      LOCI_INSERT_TYPE(cmpd, lagrangianP::NewParticle, mass) ;
      LOCI_INSERT_TYPE(cmpd, lagrangianP::NewParticle, ufluidp) ;
      return DatatypeP(cmpd) ;
    }
  } ;
  
  // here goes the stuff for list<NewParticle>
  class listNewParticleSchemaConverter {
    std::list<lagrangianP::NewParticle>& lp ;
  public:
    explicit listNewParticleSchemaConverter
    (std::list<lagrangianP::NewParticle>& lpr):lp(lpr) {}

    int
    getSize() const {
      return lp.size() ;
    }

    void
    getState(lagrangianP::NewParticle* buf, int& size) {
      size = getSize() ;
      int ii = 0 ;
      for(std::list<lagrangianP::NewParticle>::const_iterator
            ci=lp.begin();ci!=lp.end();++ci)
        buf[ii++] = *ci ;
    }
    
    void
    setState(lagrangianP::NewParticle* buf, int size) {
      lp.clear() ;
      for(int i=0;i<size;++i)
        lp.push_back(buf[i]) ;
    }
  } ;

  template<>
  struct data_schema_traits<std::list<lagrangianP::NewParticle> > {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
    typedef lagrangianP::NewParticle Converter_Base_Type ;
    typedef listNewParticleSchemaConverter Converter_Type ;
  } ;

  // Loci traits stuff for fluid_info
  template<>
  struct data_schema_traits<lagrangianP::fluid_info> {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      CompoundDatatypeP cmpd =
        CompoundFactory(lagrangianP::fluid_info()) ;
      LOCI_INSERT_TYPE(cmpd, lagrangianP::fluid_info, cellcenter) ;
      LOCI_INSERT_TYPE(cmpd, lagrangianP::fluid_info, rfluid) ;
      LOCI_INSERT_TYPE(cmpd, lagrangianP::fluid_info, ufluid) ;
      LOCI_INSERT_TYPE(cmpd, lagrangianP::fluid_info, umax) ;
      LOCI_INSERT_TYPE(cmpd, lagrangianP::fluid_info, umin) ;
      LOCI_INSERT_TYPE(cmpd, lagrangianP::fluid_info, ugrad) ;
      LOCI_INSERT_TYPE(cmpd, lagrangianP::fluid_info, viscfluid) ;
      LOCI_INSERT_TYPE(cmpd, lagrangianP::fluid_info, fvolpp) ;
      return DatatypeP(cmpd) ;
    }
  } ;
  
  template<>
  struct data_schema_traits<lagrangianP::auxiliary_info> {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      CompoundDatatypeP cmpd =
        CompoundFactory(lagrangianP::auxiliary_info()) ;
      LOCI_INSERT_TYPE(cmpd, lagrangianP::auxiliary_info, turbophoresisIndex) ;
      LOCI_INSERT_TYPE(cmpd, lagrangianP::auxiliary_info, pfluidIndex) ;
      return DatatypeP(cmpd) ;
    }
  } ;
}

namespace lagrangianP {
  // here we define a cell space class used to retrieve
  // fluid-cell related values for the particle code use
  struct CellSpace: public Singleton<CellSpace> {
    friend class Singleton<CellSpace> ;
  private:
    // we will disable the ability to copy the space
    CellSpace(const CellSpace&) ;
    CellSpace& operator=(const CellSpace&) ;
  protected:
    // flag indicates whether this is in a parallel environment or not
    bool parallel_run ;
    // constructor needs to be protected to use the Singleton stuff
    CellSpace():fluidInfo(0),fluidAuxInfo(0) {
      parallel_run = Loci::exec_current_fact_db->is_distributed_start() ;
    }
  public:
    // some cell space stores cached in
    Loci::storeRepP fluidInfo ;
    Loci::storeRepP fluidAuxInfo ;

    // this method freshes the cell space, which means to retrieve
    // the data associated with stores from the cell partition
    // through the particle space mapping information.
    void
    refresh_cell_space(Loci::storeRepP fluidInfo_rep,
                       Loci::storeRepP fluidAuxInfo_rep,
                       const ParticleSpace& pspace) ;
  } ;
    
} // end of namespace lagrangianP

#include "integrate_particle.h"
#endif
