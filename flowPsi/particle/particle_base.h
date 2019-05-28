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

/******************************************************************
NOTE: This file is developed for the purpose of separating the
original particle code into two parts: one that handles the basic
particle operations (e.g., particle walks, most of the parallel
communications, etc.), i.e., the particle geometry and computer
science stuffs; and the other one that handles the numerical
operations on particles (e.g., position integration, interaction
with fluid, etc.), i.e., the particle physics part.

Our intention is that the basic particle part does not care the
particle physics module and the basic part also cannot visit
the physics part (of course, the physics module will need to
visit the basic part since it depends on it).

For this purpose, we rely on the users to supply a particle
type of "ParticleType", the base module only deals with this
templated type and the type std::list<ParticleType> and other
containers parameterized by the particle type. However, the
base module do require that the particle type provides the
following methods defined:

 int get_id() const ; // return the particle id

 void set_id(int id) ; // set the particle id

 Entity get_cell() const ; // return the reference cell

 void set_cell(Entity c) ; // set the reference cell

 const vec3d& get_position() const ; // return particle position

 void set_position(const vec3d& p) ; // set particle position

 const vec3d& get_velocity() const ; // return particle velocity

 void set_velocity(const vec3d& v) ; // set particle velocity

method to pack the particle into a contiguous buffer.
"buffer" is the buffer to pack in, "position" is the
starting position in the buffer to begin packing,
after this method, "position" will be set to the
following position in the buffer after the packed stuff

 void pack(unsigned char* buffer, size_t& position) const ;

method to unpack a contiguous buffer into a particle record
"buffer" is the memory from which to unpack, "position"
provides the initial starting address in the buffer
and will be increased to the position immediately after
the retrieved buffer contents to facilitate next unpack

 void unpack(unsigned char* buffer, size_t& position) ;

 size_t pack_size() const ; // return a pack size

also operators "<<" and ">>" are expected to be defined
******************************************************************/

#ifndef PARTICLE_BASE_H
#define PARTICLE_BASE_H

#include "particle_config.h"
#include "util.h"
#include "par_util.h"
#include "store_traverser.h"
#include <Loci.h>
#include <Tools/except.h>
#include <iostream>
#include <set>
#include <map>
#include <list>
#include <vector>
#include <string>
#include <utility>
#include <iterator>
#include <functional>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

// redef some basic types for short names
typedef Loci::vector3d<double> vec3d ;
typedef Loci::vector2d<double> vec2d ;
#ifndef M_PI
#define M_PI            3.14159265358979323846
#endif

namespace lagrangianP {

  using Loci::Area ;
  
  typedef Loci::real_t real ;
  typedef unsigned char byte_t ;
  
  const real pi = M_PI ;
  const real EPSILON = 1e-30 ;
  
  typedef Loci::vector3d<real> vect3d ;
  typedef Loci::tensor3d<real> tens3d ;

  using Loci::vector3d ;
  using Loci::tensor3d ;
  using Loci::norm ;
  using Loci::dot ;
  using Loci::cross ;

  inline void approx_NN(const std::vector<Loci::kdTree::coord3d> &target_pnts,
                        const std::vector<int> &target_ids,
                        const std::vector<Loci::kdTree::coord3d> &search_pnts,
                        std::vector<int> &close_pt,
                        MPI_Comm comm) {
    using namespace Loci ;
    using namespace Loci::kdTree ;
    using std::vector ;
    int p = 0 ;
    int myid = 0 ;
    /* Get the number of processors and my processor identification */
    MPI_Comm_size(comm,&p) ;
    MPI_Comm_rank(comm,&myid) ;

    int lM = search_pnts.size() ;

    if(p==1) {
      kd_tree kd(target_pnts,target_ids) ;
      for(int i=0;i<lM;++i) 
        close_pt[i] = kd.find_closest(search_pnts[i]) ;
      return ;
    }

    // Sample 100000 points total max, but at least 5 per processor.
    // If target_pnts is smaller than this sample all the points
    const int sample_size = max(min(max(100000/p,5),(int)target_pnts.size()),1);
    int samp_freq = max(int(target_pnts.size()/sample_size),1) ;
    int nsamples = target_pnts.size()/samp_freq ;

    vector<double> rmin(lM,1e65) ;

    {// First sample the targets
      vector<coord3d> sample_pts(nsamples) ;
      vector<int> sample_ids(nsamples) ;
      for(int i=0;i<nsamples;++i) {
        sample_pts[i] = target_pnts[i*samp_freq] ;
        sample_ids[i] = target_ids[i*samp_freq] ;
      }

      int tsz = nsamples ;
      vector<int> rcounts(p) ;
      MPI_Allgather(&tsz,1,MPI_INT,&rcounts[0],1,MPI_INT,comm) ;
      vector<int> rdispls(p) ;
      rdispls[0] = 0 ;
      for(int i=1;i<p;++i) {
        rdispls[i] = rdispls[i-1]+rcounts[i-1] ;
      }
      int rsize = rdispls[p-1]+rcounts[p-1] ;
      
      vector<coord3d> tpnts(rsize);
      vector<int> tids(rsize) ;
      MPI_Allgatherv((void *)&sample_ids[0],tsz,MPI_INT,
                     &tids[0],&rcounts[0],&rdispls[0],
                     MPI_INT,comm) ;
      for(int i=0;i<p;++i)
        rcounts[i] *= 3 ;
      rdispls[0] = 0 ;
      for(int i=1;i<p;++i) {
        rdispls[i] = rdispls[i-1]+rcounts[i-1] ;
      }
      rsize = rdispls[p-1]+rcounts[p-1] ;
      
      
      MPI_Allgatherv((void *)&sample_pts[0],tsz*3,MPI_DOUBLE,
                     &tpnts[0],&rcounts[0],&rdispls[0],
                     MPI_DOUBLE,comm) ;
      
      kd_tree kd(tpnts,tids) ;

      for(int i=0;i<lM;++i) 
        close_pt[i] = kd.find_closest(search_pnts[i],rmin[i]) ;
    }
    
  }
  
   template <class T> class tmp_array {
    int sz ;
    T data[25] ;
    T * p ;
    void alloc(int size) {
      sz = size ;
      p = data ;
      if(sz > 25)
        p = new T[sz] ;
    }
    void free() {
      if(sz > 25)
        delete[] p ;
    }
    tmp_array() { alloc(0) ; }
  public:
    tmp_array(int size) {
      alloc(size) ;
    }
    tmp_array(const tmp_array &ta) {
      alloc(ta.sz) ;
      for(int i=0;i<sz;++i)
        p[i] = ta.p[i] ;
    }
    tmp_array &operator=(const tmp_array &ta) {
      free() ;
      alloc(ta.sz) ;
      for(int i=0;i<sz;++i)
        p[i] = ta.p[i] ;
      return *this ;
    }
    ~tmp_array() { free(); }
    T & operator[](int i) { return p[i] ; }
    T & operator[](int i) const { return p[i] ; }
    operator T *() { return p ; }
    operator const T *() const { return p ; }
  } ;

  inline double get_val_in_units(std::istream &s, const char *units) {
    double val ;
    s >> val ;
    while(s.peek() == ' ' || s.peek() == '\t')
      s.get() ;
    if(isalpha(s.peek())) {
      std::string units_in ;
      while(s.peek() != EOF && s.peek() != '\n' && s.peek() !='\r')
        units_in += s.get() ;
      Loci::UNIT_type units_value(Loci::UNIT_type::MKS,"general",val,units_in) ;
      if(!units_value.is_compatible(units)) {
        std::ostringstream oss ;
        oss << "Expecting units compatible with " << units
            << "!" << std::endl
            << "incompatible with input of '"<<units<<"'" << std::endl ;
        throw Loci::StringError(oss.str()) ;
      }
      val = units_value.get_value_in(units) ;
    }
    return val ;
  }
  
  struct TimeValue {
    double val ;
    operator double() const { return val ; }
    TimeValue &operator=(const TimeValue &v) { val = v.val ; return *this ;}
    TimeValue &operator=(double v) { val = v ; return *this; }
    TimeValue &operator=(float v) { val = v ; return *this; }
    TimeValue &operator=(int v) { val = v ; return *this; }
  } ;

  inline std::ostream & operator<<(std::ostream &s, const TimeValue &ss)
  {
    s << ss.val << " second" << std::endl ;
    return s ;
  }
    
  inline std::istream &operator>>(std::istream &s, TimeValue &ss) {
    ss.val = get_val_in_units(s,"second") ;
    return s ;
  }        
  
  // here we will define the < operator for the ParticleType.
  // this default one is based on the reference cell comparison
  template<class ParticleType> bool
  operator<(const ParticleType& pl, const ParticleType& pr) {
    return pl.get_cell() < pr.get_cell() ;
  }

  // we will be using list<particle> with Loci::blackbox so
  // define I/Os to pass the Loci traits requirements
  // these are not intended for actual use
  template<typename ParticleType> std::ostream&
  operator<<(std::ostream& s, const std::list<ParticleType>& lp) {
    for(typename std::list<ParticleType>::const_iterator
          ci=lp.begin();ci!=lp.end();++ci)
      s << *ci << std::endl ;
    return s ;
  }

  template<typename ParticleType> std::istream&
  operator>>(std::istream& s, std::list<ParticleType>& lp) {
    size_t sz ; s >> sz ;
    lp.resize(sz) ;
    for(typename std::list<ParticleType>::iterator
          ci=lp.begin();ci!=lp.end();++ci)
      s >> *ci ;
    return s ;
  }

  // structures that hold mesh cell topology for particle
  // walking. these ones are defined mainly for the incremental
  // caching strategy used for the particle walking algorithm

  // topogical data structure for interior faces in the mesh
  // each records the face id in Loci::fact_db's global numbering,
  // the center of the face, the normal of the face, and the
  // neighboring cell id in Loci::fact_db's global numbering.
  // the neighboring cell is defined to be the one on the opposite
  // side of the recorded face normal here
  struct MeshInteriorFaceTopo {
    Entity id ;                 // in global numbering
    vec3d center ;
    vec3d normal ;
    Entity neighbor ;
    int iblank; // this is the iblank value of the neighbor cell

    MeshInteriorFaceTopo(Entity fid=-1,
                         const vec3d& fc=vec3d(0,0,0),
                         const vec3d& fn=vec3d(0,0,0),
                         Entity n=-1,
                         int ib=0):
      id(fid),center(fc),normal(fn),neighbor(n),iblank(ib) {}
  } ;
  // io operators only defined to pass the compilation
  inline std::istream&
  operator>>(std::istream& s, MeshInteriorFaceTopo& f) {
    return s ;
  }
  
  inline std::ostream&
  operator<<(std::ostream& s, const MeshInteriorFaceTopo& f) {
    return s ;
  }

  inline std::istream&
  operator>>(std::istream& s,
             std::vector<MeshInteriorFaceTopo>& f) {
    return s ;
  }
  
  inline std::ostream&
  operator<<(std::ostream& s,
             const std::vector<MeshInteriorFaceTopo>& f) {
    return s ;
  }

  // topogical data structure for boundary faces in the mesh
  // each records the face id in Loci::fact_db's global numbering,
  // the center of the face, the normal of the face. we don't
  // need the neighbor cell here since it will be the boundary.
  // instead we are interested to know what type of the boundary
  // it is
  struct MeshBoundaryFaceTopo {
    Entity id ;                 // in global numbering
    vec3d center ;
    vec3d normal ;
    int wall_type ;
    int iblank; // the iblank value for a boundary is NOT actually
                // needed.  however our particle walking code
                // is setup in a generic way (paramterized by the
                // mesh topology), thus setting a dummy iblank
                // value here for the boundary topology will simplify
                // the walking code 

    MeshBoundaryFaceTopo(Entity fid=-1,
                         const vec3d& fc=vec3d(0,0,0),
                         const vec3d& fn=vec3d(0,0,0),
                         int wt=0):
      id(fid),center(fc),normal(fn),wall_type(wt),iblank(0) {}
  } ;
  // io operators only defined to pass the compilation
  inline std::istream&
  operator>>(std::istream& s, MeshBoundaryFaceTopo& f) {
    return s ;
  }
  
  inline std::ostream&
  operator<<(std::ostream& s, const MeshBoundaryFaceTopo& f) {
    return s ;
  }
  
  inline std::istream&
  operator>>(std::istream& s,
             std::vector<MeshBoundaryFaceTopo>& f) {
    return s ;
  }
  
  inline std::ostream&
  operator<<(std::ostream& s,
             const std::vector<MeshBoundaryFaceTopo>& f) {
    return s ;
  }
  
  template<class T> class ValConverter {
    T & ref ;
  public:
    ValConverter(T &iref) : ref(iref) {}
    int getSize() const {
      return 1 ;
    }
    void getState(double *buf, int &size) {
      buf[0] = ref.val ;
      size = 1 ;
    }
    void setState(double *buf, int size) {
      ref.val = buf[0] ;
    }
  } ;
  
} // end of namespace lagrangianP

// Here comes the nasty definitions for Loci traits :(

namespace Loci {
  // creating custom compound datatype for mesh topology
  // data-structure used in Loci containers
  // first the interior face data structure
  template<>
  struct data_schema_traits<lagrangianP::MeshInteriorFaceTopo> {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      CompoundDatatypeP cmpd =
        CompoundFactory(lagrangianP::MeshInteriorFaceTopo()) ;
      LOCI_INSERT_TYPE(cmpd, lagrangianP::MeshInteriorFaceTopo, id) ;
      LOCI_INSERT_TYPE(cmpd, lagrangianP::MeshInteriorFaceTopo, center) ;
      LOCI_INSERT_TYPE(cmpd, lagrangianP::MeshInteriorFaceTopo, normal) ;
      LOCI_INSERT_TYPE(cmpd, lagrangianP::MeshInteriorFaceTopo, neighbor) ;
      LOCI_INSERT_TYPE(cmpd, lagrangianP::MeshInteriorFaceTopo, iblank) ;
      return DatatypeP(cmpd) ;
    }
  } ;
  // create type for std::vector<lagrangianP::MeshInteriorFaceTopo>
  typedef std::vector<lagrangianP::MeshInteriorFaceTopo> Vmift ;

  class VmiftSchemaConverter {
    Vmift& RefObj ;
  public:
    explicit VmiftSchemaConverter(Vmift& new_obj):RefObj(new_obj) {}
    int getSize() const {
      return RefObj.size() ;
    }
    void getState(lagrangianP::MeshInteriorFaceTopo* buf, int& size) {
      size = getSize() ;
      int ii = 0 ;
      for(Vmift::const_iterator ci=RefObj.begin();
          ci!=RefObj.end();++ci)
        buf[ii++] = *ci ;
    }
    void setState(lagrangianP::MeshInteriorFaceTopo* buf, int size) {
      RefObj.resize(size) ;
      for(int i=0;i<size;++i)
        RefObj[i] = buf[i] ;
    }
  } ;

  template<>
  struct data_schema_traits<Vmift> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
    typedef lagrangianP::MeshInteriorFaceTopo Converter_Base_Type ;
    typedef VmiftSchemaConverter Converter_Type ;
  } ;

  // here goes the stuffs for the boundary face data structure
  template<>
  struct data_schema_traits<lagrangianP::MeshBoundaryFaceTopo> {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      CompoundDatatypeP cmpd =
        CompoundFactory(lagrangianP::MeshBoundaryFaceTopo()) ;
      LOCI_INSERT_TYPE(cmpd, lagrangianP::MeshBoundaryFaceTopo, id) ;
      LOCI_INSERT_TYPE(cmpd, lagrangianP::MeshBoundaryFaceTopo, center) ;
      LOCI_INSERT_TYPE(cmpd, lagrangianP::MeshBoundaryFaceTopo, normal) ;
      LOCI_INSERT_TYPE(cmpd, lagrangianP::MeshBoundaryFaceTopo, wall_type) ;
      LOCI_INSERT_TYPE(cmpd, lagrangianP::MeshBoundaryFaceTopo, iblank) ;
      return DatatypeP(cmpd) ;
    }
  } ;
  // create type for std::vector<lagrangianP::MeshBoundaryFaceTopo>
  typedef std::vector<lagrangianP::MeshBoundaryFaceTopo> Vmbft ;

  class VmbftSchemaConverter {
    Vmbft& RefObj ;
  public:
    explicit VmbftSchemaConverter(Vmbft& new_obj):RefObj(new_obj) {}
    int getSize() const {
      return RefObj.size() ;
    }
    void getState(lagrangianP::MeshBoundaryFaceTopo* buf, int& size) {
      size = getSize() ;
      int ii = 0 ;
      for(Vmbft::const_iterator ci=RefObj.begin();
          ci!=RefObj.end();++ci)
        buf[ii++] = *ci ;
    }
    void setState(lagrangianP::MeshBoundaryFaceTopo* buf, int size) {
      RefObj.resize(size) ;
      for(int i=0;i<size;++i)
        RefObj[i] = buf[i] ;
    }
  } ;
  
  template<>
  struct data_schema_traits<Vmbft> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
    typedef lagrangianP::MeshBoundaryFaceTopo Converter_Base_Type ;
    typedef VmbftSchemaConverter Converter_Type ;
  } ;

  template<> struct data_schema_traits<lagrangianP::TimeValue> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
    typedef double Converter_Base_Type ;
    typedef lagrangianP::ValConverter<lagrangianP::TimeValue> Converter_Type ;
  } ;
} // end of namespace Loci


  // Here again we will need to define some data-structures and
  // methods that are used throughout the particle base part
  // these also provide prototypes for the physics module
namespace lagrangianP {
  // the first thing we need to develop is to create a
  // mini space (sort of) for managing the particles
  // this also includes storing and managing the mesh
  // geometry data and the communication structures that
  // are useful in the particle code.
  // We will make this a Singleton object because we want
  // to control the instantiation to be only one object
  // at this moment
  template<typename T>
  class Singleton {
  public:
    static T& instance() {
      static T the_singleton_object ; // T must have a private or
      // protected constructor
      return the_singleton_object ;
    }
    virtual ~Singleton() {}
  } ;
    
  struct ParticleSpace: public Singleton<ParticleSpace> {
    friend class Singleton<ParticleSpace> ;
  private:
    // we will disable the ability to copy the space
    ParticleSpace(const ParticleSpace&) ;
    ParticleSpace& operator=(const ParticleSpace&) ;
  protected:
    // constructor needs to be protected to use the Singleton stuff
    ParticleSpace() ;
  public:
    // first structure holds the Geometry data
    struct Geometry {
    private:
      // disable copy
      Geometry(const Geometry&) ;
      Geometry& operator=(const Geometry&) ;
      // this one generates/copies the local Geometry cache
      void
      init_local_geom
      (const const_store<std::vector<MeshInteriorFaceTopo> >& mift,
       const const_store<std::vector<MeshBoundaryFaceTopo> >& mbft
#ifdef CHECK_PARTICLE_WALK
       ,const const_store<double>& cs,
       const const_store<Loci::Array<double,6> >& cb
#endif
       ) ;
      // returns the number of unoccupied slots in the buffer
      size_t
      available_slots() const {
        return cache_size - cells_cached.size() ;
      }
      // this clears the cache, but doesn't release the memory
      void
      clear() { cells_cached = EMPTY ;}
      // reinitialize the cache based on the new passed in size
      void resize(size_t sz) ;
      // this method compacts the cache based on a passed-in
      // entitySet that describes the valid cells still in the
      // cache. this method then moves all the valid cell records
      // to the beginning of the buffer and regenerates the
      // numbering maps
      void compact(const entitySet& valid) ;
      // this method does a simple caching, i.e., directly
      // get the data requested and place it in the buffer
      // starting from the index given. NOTE: it is assumed that
      // before calling this method, there will be enough space
      // in the cache to hold the new data. therefore it is
      // caller's responsibility to ensure that the buffer will
      // not overflow in this method
      void simple_cache(const entitySet& request,int start) ;
      // pointer to the Loci::fact_db
      Loci::fact_db* factsP ;
    public:
      // local Geometry buffer (fixed and unchanged)
      Loci::store<vec3d> local_cellcenter ;
      Loci::store<std::vector<MeshInteriorFaceTopo> >
      local_cell_int_face_topo ;
      Loci::store<std::vector<MeshBoundaryFaceTopo> >
      local_cell_bnd_face_topo ;
#ifdef CHECK_PARTICLE_WALK
      // these are used when particle walk check is needed
      store<double> cell_sphere ;
      store<Loci::Array<double,6> > cell_box ;
#endif
      // remote Geometry buffer (the cache)
      // cellcenters are temporarily not used
      //Loci::store<vec3d> dynamic_cellcenter ;
      Loci::store<std::vector<MeshInteriorFaceTopo> >
      dynamic_cell_int_face_topo ;
      Loci::store<std::vector<MeshBoundaryFaceTopo> >
      dynamic_cell_bnd_face_topo ;
      // things used to manage the Geometry cache
      size_t cache_size ;
      entitySet cache_domain ;
      entitySet cells_cached ;
      Map l2g ;               // local -> global conversion in the cache
      dMap g2l ;              // global -> local in the cache
      bool parallel_run ;       // indicate whether this is a parallel code
      bool do_recache ;

    public:
      Geometry()
        :factsP(0),cache_size(0),cache_domain(EMPTY),
         cells_cached(EMPTY),parallel_run(false),do_recache(false) {}
      ~Geometry() {
        local_cell_int_face_topo.allocate(EMPTY) ;
        local_cell_bnd_face_topo.allocate(EMPTY) ;
        //dynamic_cellcenter.allocate(EMPTY) ;
        dynamic_cell_int_face_topo.allocate(EMPTY) ;
        dynamic_cell_bnd_face_topo.allocate(EMPTY) ;
        l2g.allocate(EMPTY) ; g2l.allocate(EMPTY) ;          
      }
      // this is a method that caches the local face topology
      void
      cache_local
      (const const_store<std::vector<MeshInteriorFaceTopo> >& mift,
       const const_store<std::vector<MeshBoundaryFaceTopo> >& mbft,
#ifdef CHECK_PARTICLE_WALK
       const const_store<double>& cs,
       const const_store<Loci::Array<double,6> >& cb,
#endif
       Loci::fact_db* fp) {
        factsP = fp ;
        parallel_run = factsP->is_distributed_start() ;
        
        
        init_local_geom(mift, mbft
#ifdef CHECK_PARTICLE_WALK
                        ,cs, cb
#endif
                        ) ;
        // then resize the cache to be equal to the
        // local Geometry cache in the parallel run
        if(parallel_run)
          resize(local_cellcenter.domain().size()) ;
      }
      void turn_on_recache() {
        do_recache = true ;
      }
      // this method will cause all the local geometry data
      // to be recached (this is basically to account for the
      // problems of moving mesh)
      void recache_local() ;
      // this method takes in a requested cell set and caches the
      // remote stores locally in the cache. this method performs
      // an incremental caching algorithm
      void cache(const entitySet& request) ;
    } ;
    // second structure holds all particles' mapping info.
    struct MappingInfo {
      // flag to indicate whether we are in a parallel
      // or sequential environment
      bool parallel_run ;
      Loci::entitySet cells ; // all reference cells in
      // fact_db's initial global numbering
      // it is equal to g2l.domain()
      // provided here for fast access

      Loci::entitySet cells_l ; // cells mapped in a compact
      // local numbering
      // because we allocate the l2g Map
      // larger than currently used records,
      // this variable is provided for fast
      // access of the real used local
      // numbering, it is equal to
      // l2g.domain() - local_unused

      dMap g2l ;          // original global -> local numbering
      Map l2g ;         // local -> original global numbering
      // this set is used to record the unused index set in the local
      // numbering scheme, this is primarily an optimization used to
      // minimize the updating steps in the mapping information. this
      // also records the over allocated local numbering sets reserved
      // for fast append.
      Loci::entitySet local_unused ;
      // these are caches of the Loci fact_db's map and their
      // repective domain, cached here for convenience. and in
      // case of sequential run, we just build two identity maps
      dMap loci_g2l ;
      Map loci_l2g ;
      entitySet loci_g2l_dom ;
      entitySet loci_l2g_dom ;
      // constructor
      MappingInfo(bool pv=false)
        :parallel_run(pv),cells(EMPTY),
         cells_l(EMPTY),local_unused(EMPTY) {
        g2l.allocate(EMPTY) ; l2g.allocate(EMPTY) ;
      }
      ~MappingInfo() {
        g2l.allocate(EMPTY) ;
        l2g.allocate(EMPTY) ;
      }
      // this method is used to cache (or build) loci maps
      // it must be called before using this object, usually
      // the ParticleSpace is responsible to call it during
      // its initialization.
      void
      cache_loci_maps(Loci::fact_db*) ;
      // methods to update the mapping info based on a
      // sequence of particles
      template<class ParticleSequenceIterator> void
      update_mapping_info(ParticleSequenceIterator begin,
                          ParticleSequenceIterator end) ;
      // methods to update the mapping info based on a new
      // reference cells set, this method uses a more advanced
      // incremental updating algorithm.
      void
      update_mapping_info(const entitySet& new_cells) ;
      // this one adds the mapping info for the cells
      // contained in the sequence, NOTE: unlike the other
      // version, it does not remove any current mapping info
      template<class ParticleSequenceIterator> void
      add_mapping_info(ParticleSequenceIterator begin,
                       ParticleSequenceIterator end) ;
      void
      add_mapping_info(const entitySet& new_cells) ;
    } ;
    // third structure holds the communication structure.
    struct Communicator {
      friend struct ParticleSpace ;
      // first we define a communication structure for point 2 point
      // communication used in the particle program. it records
      // the domain (in entitySet) that a process needs to send to
      // others and needs to receive from others.
      struct P2pCommInfo {
        int proc ;            // process id to talk/listen to
        Loci::entitySet entities_g ;      // in global numbering
        Loci::entitySet entities_l ;      // in local numbering
        // depending on context
        // it is either the original
        // Loci::fact_db' local numbering
        // or the new local numbering
        // used here for particles
        int size ;                        // size of the entitySet
        P2pCommInfo(int p=0, const Loci::entitySet& eg=EMPTY,
                    const Loci::entitySet& el=EMPTY):proc(p),
                                                     entities_g(eg),
                                                     entities_l(el) {
          size = eg.size() ;
        }
        // ordering of the P2pCommInfo is based on the proc id
        bool
        operator<(const P2pCommInfo& p) const {
          return proc < p.proc ;
        }
        // merge operator
        P2pCommInfo&
        operator+=(const P2pCommInfo& p) {
          entities_g += p.entities_g ;
          entities_l += p.entities_l ;
          size = entities_g.size() ;
          return *this ;
        }
      } ;
      // recording communication patterns for the cell domain
      std::vector<P2pCommInfo> send ;
      std::vector<P2pCommInfo> recv ;
      // mapping info used to do pack/unpack remap during the
      // communication using the send/recv list
      // "pack" is the map used when packing the request entity set
      // during sending operations, it is usually a
      // local -> global mapping and has a type of Loci::Map.
      // the "unpack" is the map used when unpacking the request
      // entity set during the receiving operations. it is usually
      // a global -> local mapping and has a type of Loci::dMap.
      // either of them can be empty, however, users have to
      // make sure that the use of the communicator object must be
      // compatible with the pack/unpack map it contains.
      Map pack ;
      dMap unpack ;
      // these are shortcuts to the entity distribution stored
      // in the send/recv list. the "send_alloc" is the union of
      // "entities_l" in the send vector, the "recv_alloc" is the
      // union of the "entities_l" in the recv vector. they
      // represent the local numbering of the respond entitities,
      // and the local numbering of the requested entities. they
      // can be used to obtain an allocation domain.
      Loci::entitySet send_alloc ;
      Loci::entitySet recv_alloc ;
      // constructor
      Communicator():send_alloc(EMPTY),recv_alloc(EMPTY) {}
      // methods that clears a communicator object
      void
      clear() {
        send.clear(), recv.clear() ;
        send_alloc = EMPTY, recv_alloc = EMPTY ;
      }
    private:
      // the "+=" operator is defined as "merge" operations
      // it means merge in the contents to the left communicator
      // from the right side communicator. for example,
      // A += B means merge the contents in communicator "B"
      // into the communicator "A". Note: merge is simple because
      // it is a monotonically increasing operation, we will
      // never need to clean any info from the communicators.
      // on the contrary, "subtraction" is different in that it
      // may cause removing of information in the communicators,
      // and is thus harder to implement. also this "+=" is defined
      // primarily for performance reasons (incremental updates),
      // the "-=" has a less certain performance benefits. also
      // the meaning of "-=" is somewhat harder to define, for
      // example, what happens if we are subtracting something
      // that is not originally existed?
      Communicator&
      operator+=(const Communicator& c) ;
    } ;

  private:
    // there are used internally

    // pointer to Loci fact_db
    Loci::fact_db* factsP ;
    // a storage used to hold the face history for particle tracking
    // this is used for robustness issues
    std::map<int, std::set<int> >* face_history ;
    // the threshold to start particle tracking with
    // face history tracing. it means that after this number of
    // steps in the walking algorithm, we start to record faces
    // for those particles that are still walking. the default is
    // set to be 9 (meaning to begin face tracing after 9 steps)
    int face_history_tracing_start ;
    // this is used to control what is considered a "long" walk
    // i.e., we give warning messages after particles have been
    // walking for this many steps and still haven't been located
    // default is set to 500
    int max_walking_step_warning ;
    // this is how the particles are distributed, current
    // supported methods are Hilbert curve ("hilbert") or
    // an ORB method ("orb")
    std::string particle_distribution_method;
    // this dictates the dot product threshold in the walk
    // default is 0, but it is generally set to be
    // some fraction of the minimum length in the mesh
    double dot_threshold ;
    // this records one of the particle redistribution criteria
    // default is 1.1, meaning if the max is 1.1 times more than
    // the mean, we will initiate a redistribution.
    double particle_redistribution_mean_max_threshold ;
    // this records at what rate to periodically redistribute the
    // particles regardless of the per particle number per process.
    // the default is set to be 250. meaning if particles are not
    // redistributed in the last 250 geometric location calls,
    // we then redistribute it
    int particle_redistribution_freq ;
    // major structures
    Geometry geom ;
    MappingInfo cell_mapping ;
    // the pull_comm is used to receive data from the
    // fluid partition. the push_comm is used to send (usually
    // through reduction) data to the fluid partition.
    Communicator pull_comm, push_comm ;
    // this records how many particles are located
    // in the last call of the "locate" method. this
    // essentially represents the length of the particle
    // list
    size_t particle_number ;
    // this is the number used internally for assigning particle
    // ids within the ParticleSpace, during the registration
    // of new particles, the new particles will be relabeled
    // according this internal id scheme.
    long particle_id_alloc ;
    // flag to indicate whether or not we are in a parallel environment
    bool parallel_run ;
    // this variable keeps track of the total number
    // of steps used in locating all the particles last time.
    int last_location_steps ;

    // this method computes a new pull and push communicator
    // based on the passed in "domain", which indicates the
    // entitySet to receive in the pull communicator and the
    // entitySet to send in the push communicator. the "domain"
    // is required to be in Loci::fact_db's initial global
    // numbering since this is the only scheme that everyone
    // understands.
    void
    update_comm(const entitySet& domain,
                Communicator& pull, Communicator& push) ;
    // this method does not destroy the existing comms, just
    // add the domain into their send/recv lists
    void
    add_comm(const entitySet& new_domain,
             Communicator& pull, Communicator& push) ;
    // these are interface based on particle sequence
    template<class ParticleSequenceIterator> void
    update_comm(ParticleSequenceIterator begin,
                ParticleSequenceIterator end,
                Communicator& pull, Communicator& push) ;
    template<class ParticleSequenceIterator> void
    add_comm(ParticleSequenceIterator begin,
             ParticleSequenceIterator end,
             Communicator& pull, Communicator& push) ;
    /////////////////////////////////////////////////////////
    // these are methods that locate particles

    // a small auxiliary struct used to help to select
    // which face to cross in the particle walking process
    struct ParticleWalkFaceTopoAux {
      size_t idx ;
      double dt ;
      int iblank; // iblank value of the neighboring cell
      ParticleWalkFaceTopoAux(size_t i=0,double d=0,int ib=0)
        :idx(i),dt(d),iblank(ib) {}
      bool
      operator<(const ParticleWalkFaceTopoAux& pwfta) const {
        return dt < pwfta.dt ;
      }
    } ;
    // methods to determine if a particle is in a
    // cell, in case not, then compute which cell to go
    template<class MeshFaceTopo, class ParticleType> bool
    walk_particle(const ParticleType& p,
                  const std::vector<MeshFaceTopo>& face_topo,
                  size_t& face_topo_idx) ;
    template<class MeshFaceTopo, class ParticleType> bool
    walk_particle_wht(const ParticleType& p,
                      const std::vector<MeshFaceTopo>& face_topo,
                      size_t& face_topo_idx) ;
    // small function to handle particles on the bouncing wall
    template<class ParticleType> void
    bounce_particle(ParticleType& p,
                    const vec3d& an, const vec3d& fc) ;
    template<class ParticleType> void
    bounce_particle_wht(ParticleType& p,
                        const vec3d& an, const vec3d& fc) ;
    // small struct used in particle location
    struct LocateListInfo {
      entitySet ref_cells ;
      size_t len ;
      LocateListInfo(const entitySet& es=EMPTY,
                     size_t l=0):ref_cells(es),len(l) {}
    } ;
    // move particles for one step
    template<class ParticleType,
             class WalkFunInt, class WalkFunBnd, class BounceFun>
    std::pair<LocateListInfo, LocateListInfo>
    move_particles(std::list<ParticleType>& to_locate,
                   std::list<ParticleType>& located,
                   WalkFunInt walk_fun_int,
                   WalkFunBnd walk_fun_bnd,
                   BounceFun bounce_fun) ;
    // this version also gathers the bounced particles
    // (but does not process them)
    // in addition, this version also gathers any particles
    // that enter an iblanked cell for later processing
    template<class ParticleType,
             class WalkFunInt, class WalkFunBnd>
    std::pair<LocateListInfo, LocateListInfo>
    move_particles(std::list<ParticleType>& to_locate,
                   std::list<ParticleType>& located,
                   std::vector<ParticleType>& bounced,
                   std::vector<int>& bounce_face,
                   std::vector<ParticleType>& iblanked,
                   WalkFunInt walk_fun_int,
                   WalkFunBnd walk_fun_bnd) ;
    // freeze and thaw a particle list, return the number of
    // particles in the list
    template<class ParticleType> size_t
    freeze_particles(std::list<ParticleType>& lp) ;
    template<class ParticleType> size_t
    thaw_particles(std::list<ParticleType>& lp) ;
    // different interface
    template<class ParticleSequenceIterator> size_t
    freeze_particles(ParticleSequenceIterator b,
                     ParticleSequenceIterator e) ;
    template<class ParticleSequenceIterator> size_t
    thaw_particles(ParticleSequenceIterator b,
                   ParticleSequenceIterator e) ;    
  public:
    // method to initialize the local geometry data
    // in the "geom" object. this is needed to be a
    // separate method and must be called once since
    // the "geom" object needs to know the
    // "cell_int_face_topo" and "cell_bnd_face_topo"
    // facts in the Loci fact database, which may be
    // only generated after the "geom" object has been
    // created
    void
    geom_cache_local
    (const const_store<std::vector<MeshInteriorFaceTopo> >& mift,
     const const_store<std::vector<MeshBoundaryFaceTopo> >& mbft
#ifdef CHECK_PARTICLE_WALK
     ,const const_store<double>& cs,
     const const_store<Loci::Array<double,6> >& cb
#endif
     ) {
      geom.cache_local(mift,mbft,
#ifdef CHECK_PARTICLE_WALK
                       cs,cb,
#endif
                       factsP) ;
    }
    const MappingInfo&
    get_cell_mapping() const { return cell_mapping ;}
    size_t
    get_particle_number() const { return particle_number ;}
    // get the push and pull Communicator
    const Communicator&
    get_pull_comm() const { return pull_comm ; }
    const Communicator&
    get_push_comm() const { return push_comm ; }
    // this method collects all the cell in a particle sequence
    template<class ParticleSequenceIterator> entitySet
    collect_cells(ParticleSequenceIterator begin,
                  ParticleSequenceIterator end) const ;
    // here is the most important method, to provide
    // the capability to locate a list of particles
    template<class ParticleType> void
    locate_particles(std::list<ParticleType>& lp) ;
    // this is a different locate particle method.
    // it locates a list of particles, but does not
    // perform any boundary functions, instead it
    // also returns a list of particles that "bounce"
    // onto a boundary face, a separate list also
    // contains the identity of the boundary faces that
    // the particles bounce onto. The order of the
    // bounce_particles and the bounce_faces are maintained
    // as the same, i.e., the same relative position in the 
    // two lists specifies the particle and its bounce face.
    // NOTE: the returned bounced particles are deregistered
    // from the ParticleSpace and are no longer managed by
    // the ParticleSpace. The particles are distributed
    // according to their bounce face partitions, i.e.,
    // bounced particles are distributed to the process
    // that own the bounce faces. The bounce faces are
    // therefore kept in the local numbering scheme of the
    // Loci's initial global partition, however the cells
    // that the particles are referring to inside the
    // particle data-structure are in the global numbering.
    // This is so since we do not know if the cells and
    // the bounce faces are in the same partition.
    //
    // currently, we also collect particles that enter iblanked
    // cells.  the same logic is used in the collected vector
    // of particles in iblaned cell: these particles are removed
    // from the management of the current ParticleSpace and then
    // redistributed to the cell partition (the iblanked cells they
    // are in), however the iblanked cell they are in are also
    // converted to Loci's local numbering, this is mainly to 
    // facilitate the later processing of these particles
    template<class ParticleType> void
    locate_particles(std::list<ParticleType>& lp,
                     std::vector<ParticleType>& bounce_particles,
                     std::vector<int>& bounce_faces,
                     std::vector<ParticleType>& iblanked_particles) ;
    // this is a internal routine used by the above function
    // and also by the "register_particles" function to check
    // if any newly registered particls are not properly located.
    // this is mainly just to separate the core part of the
    // location code so that it may be reused in different parts.
    // returns the total number of steps walked.
    template<class ParticleType> int
    locate_particles_core(std::list<ParticleType>& to_locate,
                          std::list<ParticleType>& located,
                          std::vector<ParticleType>& bounce_particles,
                          std::vector<int>& bounce_faces,
                          std::vector<ParticleType>& iblanked_particles,
                          entitySet& new_cells,
                          size_t& result_located_len);

    // newly added particles must be registered in
    // the ParticleSpace using this method.
    // the new particles are included in the sequence
    // defined within [first, last). and the newly
    // registered particles will be created in the
    // list "rp". the particles in the sequence are
    // assumed to have their reference cells in the
    // Loci fact_db's initial global numbering
    template<class ParticleType, class ForwardIterator> void
    register_particles(ForwardIterator first,
                       ForwardIterator last,
                       std::list<ParticleType>& rp) ;
    // this is another registering method used for new particles
    // stored in a Loci::store<std::list<ParticleType>
    // it constructs a new list of particles in the second parameter
    // the particles contained in the store are assumed to have
    // their cells in the Loci fact_db's local numbering
    template<class ParticleType> void
    register_particles(Loci::storeRepP sp,
                       std::list<ParticleType>& rp) ;
    // this is another registering method primary used to read
    // in a particle restart file, whose name is passed in as a string
    // after this method, particles in the restart file will be
    // reconstructed in the passed in list "rp"
    // 
    // THIS IS NOW DEPRECATED, use the restart_particles routine!
    template<class ParticleType> void
    register_particles(std::string restart_file,
                       std::list<ParticleType>& rp) ;
    // this is a slow restart procedure used only for comparison purpose
    template<class ParticleType> void
    restart_particles_slow(std::string restart_file,
                           std::list<ParticleType>& rp);
    // this is a new restart procedure that is based on file numbers.
    template<class ParticleType> void
    restart_particles(std::string restart_file,
                      std::list<ParticleType>& rp);
    
    // THIS METHOD IS BEING PHASING OUT...
    template<class ParticleType> void
    insert_new_particles(std::vector<ParticleType> &particles_to_insert,
                       std::list<ParticleType>& rp) ;
    // some methods to set the thresholds
    void
    set_face_history_tracing_start(int s) {
      face_history_tracing_start = s ;
    }
    void
    set_max_walking_step_warning(int s) {
      max_walking_step_warning = s ;
    }
    void
    set_particle_distribution_method(std::string s) {
      particle_distribution_method = s;
    }
    void
    set_dot_threshold(double s) {
      dot_threshold = s ;
    }
    void
    set_particle_dist_mm_threshold(double s) {
      particle_redistribution_mean_max_threshold = s ;
    }
    void
    set_particle_dist_freq(int s) {
      particle_redistribution_freq = s ;
    }
    void turn_on_recache() {
      geom.turn_on_recache() ;
    }
    int
    get_last_location_steps() const {
      return last_location_steps ;
    }
      
  } ; // end of ParticleSpace
  
    
  ////////////////////////////////////////////////////////////////
  // this section implements some ParticleSpace::MappingInfo
  // methods that have templates
  ////////////////////////////////////////////////////////////////
  template<class ParticleSequenceIterator> void
  ParticleSpace::MappingInfo::update_mapping_info
  (ParticleSequenceIterator begin, ParticleSequenceIterator end) {
    // first get all the reference cells
    Entity cur_ref_cell = -1 ;
    std::vector<Entity> tmp ;
    for(;begin!=end;++begin) {
      if(begin->get_cell() != cur_ref_cell) {
        cur_ref_cell = begin->get_cell() ;
        tmp.push_back(cur_ref_cell) ;
      }
    }
    entitySet new_cells =
      Loci::create_intervalSet(tmp.begin(),tmp.end()) ;
    update_mapping_info(new_cells) ;
  }

  template<class ParticleSequenceIterator> void
  ParticleSpace::MappingInfo::add_mapping_info
  (ParticleSequenceIterator begin, ParticleSequenceIterator end) {
    // first get all the reference cells
    Entity cur_ref_cell = -1 ;
    std::vector<Entity> tmp ;
    for(;begin!=end;++begin) {
      if(begin->get_cell() != cur_ref_cell) {
        cur_ref_cell = begin->get_cell() ;
        tmp.push_back(cur_ref_cell) ;
      }
    }
    entitySet new_cells =
      Loci::create_intervalSet(tmp.begin(),tmp.end()) ;
    // since we don't want to remove the current mapping info,
    // we just merge the new_cells into the current cells
    new_cells += cells ;
    update_mapping_info(new_cells) ;
  }

  inline void
  ParticleSpace::MappingInfo::add_mapping_info
  (const entitySet& new_cells) {
    entitySet cs = new_cells + cells ;
    update_mapping_info(cs) ;
  }

  ////////////////////////////////////////////////////////////////
  // this section implements ParticleSpace methods, they need
  // to stay in the header file because they have templates
  ////////////////////////////////////////////////////////////////
  template<class ParticleSequenceIterator> entitySet
  ParticleSpace::collect_cells(ParticleSequenceIterator begin,
                               ParticleSequenceIterator end) const {
    Entity cur_ref_cell = -1 ;
    std::vector<Entity> tmp ;
    for(;begin!=end;++begin) {
      if(begin->get_cell() != cur_ref_cell) {
        cur_ref_cell = begin->get_cell() ;
        tmp.push_back(cur_ref_cell) ;
      }
    }
    return Loci::create_intervalSet(tmp.begin(),tmp.end()) ;
  }
  
  template<class ParticleSequenceIterator> void
  ParticleSpace::update_comm(ParticleSequenceIterator begin,
                             ParticleSequenceIterator end,
                             Communicator& pull, Communicator& push) {
    // first get all the reference cells
    Entity cur_ref_cell = -1 ;
    std::vector<Entity> tmp ;
    for(;begin!=end;++begin) {
      if(begin->get_cell() != cur_ref_cell) {
        cur_ref_cell = begin->get_cell() ;
        tmp.push_back(cur_ref_cell) ;
      }
    }
    entitySet new_cells =
      Loci::create_intervalSet(tmp.begin(),tmp.end()) ;
    update_comm(new_cells, pull, push) ;
  }
  
  template<class ParticleSequenceIterator> void
  ParticleSpace::add_comm(ParticleSequenceIterator begin,
                          ParticleSequenceIterator end,
                          Communicator& pull, Communicator& push) {
    // first get all the reference cells
    Entity cur_ref_cell = -1 ;
    std::vector<Entity> tmp ;
    for(;begin!=end;++begin) {
      if(begin->get_cell() != cur_ref_cell) {
        cur_ref_cell = begin->get_cell() ;
        tmp.push_back(cur_ref_cell) ;
      }
    }
    entitySet new_cells =
      Loci::create_intervalSet(tmp.begin(),tmp.end()) ;
    add_comm(new_cells, pull, push) ;
  }

  // given a set of face topology info. and a particle,
  // this function decides whether the particle is "outside"
  // the face set, it returns "true" if the particle is
  // outside the supplied face set, "false" if inside.
  // in case of outside, it also sets the parameter
  // "face_topo_idx", which is an index number to the
  // supplied "face_topo" vector that indicates which
  // face the particle should cross in the current walk.
  // NOTE: we use a template<class MeshFaceTopo> to
  // indicate the type of the supplied face topology.
  // This is mainly to reuse the code since we have two
  // types of face topology, the interior and boundary types.
  template<class MeshFaceTopo, class ParticleType> bool
  ParticleSpace::walk_particle
  (const ParticleType& p,
   const std::vector<MeshFaceTopo>& face_topo, size_t& face_topo_idx) {
    bool walk = false ;
    std::vector<ParticleWalkFaceTopoAux> walk_aux ;
    for(size_t i=0;i<face_topo.size();++i) {
      int iblank = face_topo[i].iblank;
      vec3d v1 = p.get_position() - face_topo[i].center ;
      //      normalize(v1) ;
      const vec3d& v2 = face_topo[i].normal ;
      double dt = dot(v1,v2) ;
      if(dt < -dot_threshold) {
        walk_aux.push_back(ParticleWalkFaceTopoAux(i,dt,iblank)) ;
        walk = true ;
      }
    }
    if(walk) {
      // choose the face that generated the most negative dot product
      std::sort(walk_aux.begin(),walk_aux.end()) ;
      // after sorting, we will want to avoid any path that has
      // a iblank neighboring cell, unless we do not have any
      // other choices, so we will do a search through the sorted
      // list of potential faces to cross and pick the one that
      // does not have a iblanked neighbor
      size_t choice = 0;
      for(size_t i=0;i<walk_aux.size();++i)
        if(walk_aux[i].iblank < 2) {
          choice = i;
          break;
        }
      face_topo_idx = walk_aux[choice].idx ;
    }
    return walk;
  }

  // version with face history checking
  template<class MeshFaceTopo, class ParticleType> bool
  ParticleSpace::walk_particle_wht
  (const ParticleType& p,
   const std::vector<MeshFaceTopo>& face_topo, size_t& face_topo_idx) {
    bool walk = false ;
    std::vector<ParticleWalkFaceTopoAux> walk_aux(face_topo.size()) ;
    std::vector<bool> face_mark(face_topo.size(),false) ;
    // retrieve the face history if there is one there
    std::map<int,std::set<int> >::iterator
      mi = face_history->find(p.get_id()) ;
    std::set<int>* visited_faces = 0 ;
    if(mi != face_history->end()) {
      visited_faces = &(mi->second) ;
      for(size_t i=0;i<face_topo.size();++i) {
        int face = face_topo[i].id ;
        if(visited_faces->find(face) != visited_faces->end())
          face_mark[i] = true ;
      }
    }
    for(size_t i=0;i<face_topo.size();++i) {
      // if a face is marked before, we don't cross it
      if(face_mark[i])
        continue ;
      int iblank = face_topo[i].iblank;
      vec3d v1 = p.get_position() - face_topo[i].center ;
      //      normalize(v1) ;
      const vec3d& v2 = face_topo[i].normal ;
      double dt = dot(v1,v2) ;
      if(dt < -dot_threshold) {
        walk_aux.push_back(ParticleWalkFaceTopoAux(i,dt,iblank)) ;
        walk = true ;
      }
    }
    if(walk) {
      std::sort(walk_aux.begin(),walk_aux.end()) ;
      // see comments int the above function for the logic of this part
      size_t choice = 0;
      for(size_t i=0;i<walk_aux.size();++i)
        if(walk_aux[i].iblank < 2) {
          choice = i;
          break;
        }
      face_topo_idx = walk_aux[choice].idx ;
      // record the face history since we will cross
      // the face face_topo[face_topo_idx].id
      if(visited_faces)
        visited_faces->insert(face_topo[face_topo_idx].id) ;
      else
        (*face_history)[p.get_id()].
          insert(face_topo[face_topo_idx].id) ;
    }
    return walk ;
  }

  // functions that process the particles on the bouncing wall
  template<class ParticleType> inline void
  ParticleSpace::bounce_particle(ParticleType& p,
                                 const vec3d& an, const vec3d& fc) {
    vec3d pos = p.get_position() ;
    vec3d vel = p.get_velocity() ;
    
    p.set_position(pos - ( (2.0 * dot(pos-fc, an)) * an)) ;
    p.set_velocity(vel - ( (2.0 * dot(vel, an)) * an)) ;
  }

  template<class ParticleType> inline void
  ParticleSpace::bounce_particle_wht(ParticleType& p,
                                     const vec3d& an, const vec3d& fc) {
    vec3d pos = p.get_position() ;
    vec3d vel = p.get_velocity() ;
    
    p.set_position(pos - ( (2.0 * dot(pos-fc, an)) * an)) ;
    p.set_velocity(vel - ( (2.0 * dot(vel, an)) * an)) ;

    // we will need to erase history for bounced particles
    face_history->erase(p.get_id()) ;
  }
  
  // a function that takes a list<particle> and locates it.
  // this is the major function that does the particle location
  // along with all of the necessary parallel communications
  //
  // This function only does one step of locating, i.e., a
  // particle is either in the current cell, or it is moved
  // to a neighbor cell for further processing.
  // 
  // after the function, the passed in parameter "to_locate"
  // contains particles that are not yet located and are moved
  // to another cell. "located" contains all the particles that
  // are located.
  // It returns a pair of list info, the first is the list info
  // for the particles in the "to_locate" list upon
  // completing this function, the second info is for
  // the particles in the "located" list upon
  // completing this function. the return object is created
  // primarily to speed up the list and entitySet processing,
  // since list.size() may take linear time and collecting
  // a large entitySet may take a long time
  // 
  // the parameter "WalkFunInt" is the type for the
  // "walk_particle" function used inside when applied to
  // MeshInteriorFaceTopo. the parameter "WalkFunBnd" is
  // the type for the "walk_particle" function used inside
  // when applied to MeshBoundaryFaceTopo. it is parameterized
  // because we have several versions of it, one normal version,
  // the other one with face history checking. doing so can
  // save us some coding work (instead of having several versions
  // of "move_particles" function, we can just have one)
  //
  // the parameter "BounceFun" is the type for the function
  // that processes particles that bounce onto the wall boundary.
  // doing so has the same reason as using "WalkFun"
  template<class ParticleType, class WalkFunInt,
           class WalkFunBnd, class BounceFun>
  std::pair<ParticleSpace::LocateListInfo,
            ParticleSpace::LocateListInfo>
  ParticleSpace::move_particles
  (std::list<ParticleType>& to_locate, std::list<ParticleType>& located,
   WalkFunInt walk_fun_int, WalkFunBnd walk_fun_bnd, BounceFun bounce_fun) {

    entitySet to_locate_ref_cells = EMPTY ;
    size_t to_locate_len = 0 ;
    entitySet located_ref_cells = EMPTY ;
    size_t located_len = 0 ;

    store<std::vector<MeshInteriorFaceTopo> >
      cell_int_face_topo(geom.dynamic_cell_int_face_topo.Rep()) ;
    store<std::vector<MeshBoundaryFaceTopo> >
      cell_bnd_face_topo(geom.dynamic_cell_bnd_face_topo.Rep()) ;

    typename std::list<ParticleType>::iterator li = to_locate.begin() ;

    Entity cur_ref_cell = -1 ;
    Entity c = -1 ;
    
    while(li != to_locate.end()) {
      // get the cell number
      if(li->get_cell() != cur_ref_cell) {
        cur_ref_cell = li->get_cell() ;
        c = geom.g2l[cur_ref_cell] ;
      }
      size_t face_topo_idx ;
      // first walk the particle with the interior faces
      // we use the parameterized walking function
      if( (this->*walk_fun_int)(*li,
                                cell_int_face_topo[c], face_topo_idx)) {
        Entity neighbor =
          cell_int_face_topo[c][face_topo_idx].neighbor ;
        li->set_cell(neighbor) ;
        ++li ;
        // update return object
        to_locate_ref_cells += neighbor ;
        ++to_locate_len ;
        
        continue ;
      }
      // if the particle doesn't walk through interior
      // faces, we check to see if it will hit any boundary face
      // this walk function does not need to use the
      // parameterized one because we don't need to do
      // face history checking for any boundary faces
      // because we know particles must have not crossed
      // any boundary face before.
      //      if( dot(li->get_velocity(),cell_bnd_face_topo[c][face_topo_idx].normal) >
      if( (this->*walk_fun_bnd)(*li,
                                cell_bnd_face_topo[c], face_topo_idx)) {
        // if we reach here, that means the particle is only
        // outside of the boundary face (since we've run through
        // all the interior faces before). in this case, the
        // particle falls outside of the computational domain
        // we process it according to the boundary wall type.
        int wall_bc = cell_bnd_face_topo[c][face_topo_idx].wall_type ;
        if(wall_bc < 2) {
          // remove the particle
          typename std::list<ParticleType>
            ::iterator li2 = li ;
          ++li2 ;
          to_locate.erase(li) ;
          li = li2 ;
        } else {
          // bounce wall
          const vec3d& an =
            cell_bnd_face_topo[c][face_topo_idx].normal ;
          const vec3d& fc =
            cell_bnd_face_topo[c][face_topo_idx].center ;

          (this->*bounce_fun)(*li, an, fc) ;

          ++li ;

          // update return object
          to_locate_ref_cells += cur_ref_cell ;
          ++to_locate_len ;          
        }
        continue ;
      }
      // if we reach here, the particle is still within the current
      // reference cell, we need to move it to the located list
      typename std::list<ParticleType>
        ::iterator li_bak = li ;
      ++li_bak ;
      located.splice(located.end(), to_locate, li) ;
      li = li_bak ;
      // update the return object
      located_ref_cells += cur_ref_cell ;
      ++located_len ;
    }
    return std::make_pair
      (LocateListInfo(to_locate_ref_cells,to_locate_len),
       LocateListInfo(located_ref_cells,located_len)) ;
  }

  template<class ParticleType,
           class WalkFunInt, class WalkFunBnd>
  std::pair<ParticleSpace::LocateListInfo,
            ParticleSpace::LocateListInfo>
  ParticleSpace::move_particles(std::list<ParticleType>& to_locate,
                                std::list<ParticleType>& located,
                                std::vector<ParticleType>& bounced,
                                std::vector<int>& bounce_face,
                                std::vector<ParticleType>& iblanked,
                                WalkFunInt walk_fun_int,
                                WalkFunBnd walk_fun_bnd) {
    entitySet to_locate_ref_cells = EMPTY ;
    size_t to_locate_len = 0 ;
    entitySet located_ref_cells = EMPTY ;
    size_t located_len = 0 ;

    store<std::vector<MeshInteriorFaceTopo> >
      cell_int_face_topo(geom.dynamic_cell_int_face_topo.Rep()) ;
    store<std::vector<MeshBoundaryFaceTopo> >
      cell_bnd_face_topo(geom.dynamic_cell_bnd_face_topo.Rep()) ;

    typename std::list<ParticleType>::iterator li = to_locate.begin() ;

    Entity cur_ref_cell = -1 ;
    Entity c = -1 ;
    
    while(li != to_locate.end()) {
      // get the cell number
      if(li->get_cell() != cur_ref_cell) {
        cur_ref_cell = li->get_cell() ;
        c = geom.g2l[cur_ref_cell] ;
      }
      size_t face_topo_idx ;
      // first walk the particle with the interior faces
      // we use the parameterized walking function
      if( (this->*walk_fun_int)(*li,
                                cell_int_face_topo[c], face_topo_idx)) {
        Entity neighbor =
          cell_int_face_topo[c][face_topo_idx].neighbor ;
        li->set_cell(neighbor) ;
        // check to see if the particle is inside an iblanked cell:
        // but we'll give a warning if the iblank value of the 
        // containing cell is beyond 2:
        int iblank = cell_int_face_topo[c][face_topo_idx].iblank;
        if(iblank > 2) {
          if(Loci::MPI_rank == 0) {
            std::cout << "*****ParticleSpace Message*****" << std::endl ;
            std::cout << "A particle enters a cell with "  << std::endl;
            std::cout << "iblank value: " << iblank << std::endl ;
            std::cout << "--the particle is being removed!" << std::endl;
            std::cout << "*******************************" << std::endl ;
          }
          typename std::list<ParticleType>::iterator li_bak = li;
          ++li_bak;
          to_locate.erase(li);
          li = li_bak;
        } else if(iblank > 1) {
          typename std::list<ParticleType>::iterator li_bak = li;
          ++li_bak;
          iblanked.push_back(*li);
          to_locate.erase(li);
          li = li_bak;
        } else {
          ++li ;
          // update return object
          to_locate_ref_cells += neighbor ;
          ++to_locate_len ;
        }
        
        continue ;
      }
      // if the particle doesn't walk through interior
      // faces, we check to see if it will hit any boundary face
      // this walk function does not need to use the
      // parameterized one because we don't need to do
      // face history checking for any boundary faces
      // because we know particles must have not crossed
      // any boundary face before.
      if( (this->*walk_fun_bnd)(*li,
                                cell_bnd_face_topo[c], face_topo_idx)) {
        // if we reach here, that means the particle is only
        // outside of the boundary face (since we've run through
        // all the interior faces before). in this case, the
        // particle falls outside of the computational domain
        // we record it and remove it from the to_locate list
        int bf = cell_bnd_face_topo[c][face_topo_idx].id ;
        typename std::list<ParticleType>
          ::iterator li_bak = li ;
        ++li_bak ;
        bounced.push_back(*li) ;
        bounce_face.push_back(bf) ;
        to_locate.erase(li) ;
        li = li_bak ;

        continue ;
      }
      // if we reach here, the particle is still within the current
      // reference cell, we need to move it to the located list
      typename std::list<ParticleType>
        ::iterator li_bak = li ;
      ++li_bak ;
      located.splice(located.end(), to_locate, li) ;
      li = li_bak ;
      // update the return object
      located_ref_cells += cur_ref_cell ;
      ++located_len ;
    }
    return std::make_pair
      (LocateListInfo(to_locate_ref_cells,to_locate_len),
       LocateListInfo(located_ref_cells,located_len)) ;
  }
  
  template<class ParticleType> size_t
  ParticleSpace::freeze_particles(std::list<ParticleType>& lp) {
    return freeze_particles(lp.begin(), lp.end()) ;
  }
  
  template<class ParticleType> size_t
  ParticleSpace::thaw_particles(std::list<ParticleType>& lp) {
    return thaw_particles(lp.begin(), lp.end()) ;
  }
  
  template<class ParticleSequenceIterator> size_t
  ParticleSpace::freeze_particles(ParticleSequenceIterator b,
                                  ParticleSequenceIterator e) {
    // freeze is global -> local number mapping
    size_t sz = 0 ;
    if(parallel_run) {
      // freeze is only meaningful in a parallel environment
      for(;b!=e;++b,++sz) {
        Entity cell = b->get_cell() ;
        b->set_cell(cell_mapping.g2l[cell]) ;
      }
    } else
      sz = std::distance(b,e) ;

    return sz ;
  }

  template<class ParticleSequenceIterator> size_t
  ParticleSpace::thaw_particles(ParticleSequenceIterator b,
                                ParticleSequenceIterator e) {
    // thaw is local -> global number mapping
    size_t sz = 0 ;
    if(parallel_run) {
      // thaw is only meaningful in a parallel environment
      for(;b!=e;++b,++sz) {
        Entity cell = b->get_cell() ;
        b->set_cell(cell_mapping.l2g[cell]) ;
      }
    } else
      sz = std::distance(b,e) ;

    return sz ;
  }

  template<class ParticleType> void
  ParticleSpace::locate_particles(std::list<ParticleType>& lp) {
    ////////////////////////////////////////////////////
    // the local geometry cache is flushed and recached
    // if the mesh data is somehow changing between
    // particle location steps.
    ////////////////////////////////////////////////////
    if(geom.do_recache) {
      if(Loci::MPI_rank == 0) {
        std::cout << "*****ParticleSpace Message*****" << std::endl ;
        std::cout << " Re-caching all geometric data" << std::endl ;
        std::cout << "*******************************" << std::endl ;
      }
      geom.recache_local() ;
    }
    
    // this variable is used to track how long since the
    // last particle distribution happened.
    static int distribution_tick = 0 ;
    
    std::list<ParticleType>& to_locate = lp ;
    std::list<ParticleType> located ;

    // the first thing we need to do is to thaw the particle to cell map
    thaw_particles(lp) ;

    int current_walking_step = 0 ;

    // this records all the new reference cells
    // when the location is complete
    entitySet new_cells ;

    // initial cell requests are those in the current cell_mapping
    entitySet request_cells = cell_mapping.cells ;

    // this records the final located particle number
    size_t result_located_len = 0 ;

    // start the location process
    while(true) {
      // first we will to if we are walking too many steps.
      if(current_walking_step > max_walking_step_warning) {
        if(Loci::MPI_rank == 0) {
          std::cout << "*****ParticleSpace Warning*****" << std::endl ;
          std::cout << "* Current walking step > "
                    << max_walking_step_warning << std::endl ;
        }
        int ls = to_locate.size() ;
        int gls = 0 ;
        MPI_Reduce(&ls, &gls, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD) ;
        if(Loci::MPI_rank == 0) {
          std::cout << "* Code may have problems!" << std::endl ;
          std::cout << "* Removing the remaining " << gls
                    << " particles" << std::endl ;
          std::cout << "*******************************" << std::endl ;
        }
        typename std::list<ParticleType>::const_iterator ip ;
        for(ip=to_locate.begin();ip != to_locate.end();++ip) {
          Loci::debugout << "P:pos = " << ip->pos[0] << "::"
                         << ip->pos[1] <<std::endl ;
        }

        to_locate.clear() ;
        
      }

      // check whether walking is completed
      if(Loci::GLOBAL_AND(to_locate.empty()))
        break ;

      // first we need to cache the requested cells
      geom.cache(request_cells) ;

      std::pair<LocateListInfo,LocateListInfo> walk_result ;

      if(current_walking_step < face_history_tracing_start) {
        // walk particles without face history checking
        walk_result = move_particles
          (to_locate, located,
           &ParticleSpace::walk_particle<MeshInteriorFaceTopo,
                                         ParticleType>,
           &ParticleSpace::walk_particle<MeshBoundaryFaceTopo,
                                         ParticleType>,
           &ParticleSpace::bounce_particle<ParticleType>) ;
        
      } else {
        // walk particles with face history checking
        if(!face_history)
          face_history = new std::map<int, std::set<int> > ;

        walk_result = move_particles
          (to_locate, located,
           &ParticleSpace::walk_particle_wht<MeshInteriorFaceTopo,
                                             ParticleType>,
           &ParticleSpace::walk_particle<MeshBoundaryFaceTopo,
                                         ParticleType>,
           &ParticleSpace::bounce_particle_wht<ParticleType>) ;
      }

      // finished walking one step
      ++current_walking_step ;
      // modifying the necessary info according the walk_result
      new_cells += walk_result.second.ref_cells ;
      request_cells = walk_result.first.ref_cells ;
      result_located_len += walk_result.second.len ;
    } // end while(true)

    if(face_history) {
      delete face_history ;
      face_history = 0 ;
    }
    particle_number = result_located_len ;
    last_location_steps = current_walking_step ;
    // location completed, we need to see whether or not
    // to perform a particle distribution
    ++distribution_tick ;

    bool need_distribution = false ;
    std::string distribution_type ;

    // in case of a parallel ParticleSpace, we need to
    // check to see whether we need to perform a particle distribution
    if(parallel_run) {
      // first we will see if the mean max criterion is met
      // divide first so we don't overflow
      int pmeanp = (particle_number+Loci::MPI_processes/2)/Loci::MPI_processes ;
      // sum of pmeanp over all processors approximates the mean size
      mpi_2int pair_l(pmeanp, particle_number), pair_g ;
      MPI_Op mpi_2int_op ;
      MPI_Op_create( (MPI_User_function*)SUM_MAX_2INT, 1, &mpi_2int_op) ;
      MPI_Allreduce(&pair_l, &pair_g, 1,
                    MPI_2INT, mpi_2int_op, MPI_COMM_WORLD) ;
      MPI_Op_free(&mpi_2int_op) ;
      double mean_size =  static_cast<double>(pair_g.first)  ;
      int max_size = pair_g.second ;
      
      if(max_size > 100 && static_cast<double>(max_size) >
         particle_redistribution_mean_max_threshold*mean_size) {
        need_distribution = true ;
        distribution_type = "auto" ;
      }
      // if the mean max criterion is not met, then we see
      // if the pre-set distribution interval is met
      if(!need_distribution) {
        if(distribution_tick >= particle_redistribution_freq) {
          need_distribution = true ;
          distribution_type = "pre" ;
        }
      }
    } // end if(parallel_run)

    if(need_distribution) {
      distribution_tick = 0 ;

      if(Loci::MPI_rank == 0) {
        std::cout << "*****ParticleSpace Message*****" << std::endl ;
        std::cout << "* Distribute particle list ["
                  << distribution_type << "]" ;
      }
      // since the list to_locate points to the
      // passed in list and that it is now empty,
      // we will put the distribution result in it
      std::back_insert_iterator<std::list<ParticleType> >
        result(to_locate) ;

      if(particle_distribution_method == "orb") {
        if(Loci::MPI_rank == 0) {
          std::cout << " (orb) ... " ;
          std::cout.flush() ;
        }
        orb_partition_sequence(located.begin(), located.end(),
                               result, MPI_COMM_WORLD) ;
      } else {
        if(Loci::MPI_rank == 0) {
          std::cout << " (hilbert) ... " ;
          std::cout.flush() ;
        }
        hilbert_partition_sequence(located.begin(), located.end(),
                                   result, MPI_COMM_WORLD);
      }
      
      if(Loci::MPI_rank ==0) {
        std::cout << "Done" << std::endl ;
        std::cout << "*******************************" << std::endl ;
      }
      // we can now get rid of the located list
      std::list<ParticleType>().swap(located) ;
      // we will need to update the mapping info
      new_cells = collect_cells(to_locate.begin(), to_locate.end()) ;
      cell_mapping.update_mapping_info(new_cells) ;
      // update the pull and push comm
      update_comm(new_cells, pull_comm, push_comm) ;
      particle_number = to_locate.size() ;
      
    } else {
      // copy the located list to the to_locate list
      to_locate.splice(to_locate.end(), located) ;
      // don't distribute particles, just update the mapping info
      cell_mapping.update_mapping_info(new_cells) ;
      // updating the pull and push comm
      update_comm(new_cells, pull_comm, push_comm) ;
    }

    // we then proceed to freeze the particle to cell map again
    freeze_particles(to_locate) ;

    // finally we perform a sorting of the particle list
    // based on their reference cells
    to_locate.sort() ;
  }
  
  template<class ParticleType> void
  ParticleSpace::locate_particles(std::list<ParticleType>& lp,
                                  std::vector<ParticleType>& bounce_particles,
                                  std::vector<int>& bounce_faces,
                                  std::vector<ParticleType>&
                                    iblanked_particles) {
    ////////////////////////////////////////////////////
    // the local geometry cache is flushed and recached
    // if the mesh data is somehow changing between
    // particle location steps.
    ////////////////////////////////////////////////////
    if(geom.do_recache) {
      if(Loci::MPI_rank == 0) {
        std::cout << "*****ParticleSpace Message*****" << std::endl ;
        std::cout << " Re-caching all geometric data" << std::endl ;
        std::cout << "*******************************" << std::endl ;
      }
      geom.recache_local() ;
    }
    
    // this variable is used to track how long since the
    // last particle distribution happened.
    static int distribution_tick = 0 ;
    
    std::list<ParticleType>& to_locate = lp ;
    std::list<ParticleType> located ;

    // the first thing we need to do is to thaw the particle to cell map
    thaw_particles(lp) ;

#ifdef DONOTREMOVE
    int current_walking_step = 0 ;
#endif

    // this records all the new reference cells
    // when the location is complete
    entitySet new_cells ;

#ifdef DONOTREMOVE
    // initial cell requests are those in the particle refs
    entitySet request_cells = cell_mapping.cells ;
#endif

    // this records the final located particle number
    size_t result_located_len = 0 ;

    int walked_steps =
      locate_particles_core(to_locate,located,bounce_particles,
                            bounce_faces,iblanked_particles,
                            new_cells,result_located_len);

#ifdef DONOTREMOVE
    // start the location process
    while(true) {
      // first we will to if we are walking too many steps.
      if(current_walking_step > max_walking_step_warning) {
        if(Loci::MPI_rank == 0) {
          std::cout << "*****ParticleSpace Warning*****" << std::endl ;
          std::cout << "* Current walking step > "
                    << max_walking_step_warning << std::endl ;
        }
        int ls = to_locate.size() ;
        int gls = 0 ;
        MPI_Reduce(&ls, &gls, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD) ;
        if(Loci::MPI_rank == 0) {
          std::cout << "* Code may have problems!" << std::endl ;
          std::cout << "* Removing the remaining " << gls
                    << " particles" << std::endl ;
          std::cout << "*******************************" << std::endl ;
        }
        //        typename std::list<ParticleType>::const_iterator ip ;
        //        for(ip=to_locate.begin();ip != to_locate.end();++ip) {
        //          Loci::debugout << "P:pos = " << ip->pos[0] << "::"
        //                         << ip->pos[1] <<std::endl ;
        //        }
        to_locate.clear() ;
      }

      // check whether walking is completed
      if(Loci::GLOBAL_AND(to_locate.empty()))
        break ;

      // first we need to cache the requested cells
      geom.cache(request_cells) ;

      std::pair<LocateListInfo,LocateListInfo> walk_result ;

      if(current_walking_step < face_history_tracing_start) {
        // walk particles without face history checking
        walk_result = move_particles
          (to_locate, located, 
           bounce_particles, bounce_faces, iblanked_particles,
           &ParticleSpace::walk_particle<MeshInteriorFaceTopo,
                                         ParticleType>,
           &ParticleSpace::walk_particle<MeshBoundaryFaceTopo,
                                         ParticleType>) ;
        
      } else {
        // walk particles with face history checking
        if(!face_history)
          face_history = new std::map<int, std::set<int> > ;

        walk_result = move_particles
          (to_locate, located, 
           bounce_particles, bounce_faces, iblanked_particles,
           &ParticleSpace::walk_particle_wht<MeshInteriorFaceTopo,
                                             ParticleType>,
           &ParticleSpace::walk_particle<MeshBoundaryFaceTopo,
                                         ParticleType>) ;
      }

      // finished walking one step
      ++current_walking_step ;
      // modifying the necessary info according the walk_result
      new_cells += walk_result.second.ref_cells ;
      request_cells = walk_result.first.ref_cells ;
      result_located_len += walk_result.second.len ;
    } // end while(true)

    if(face_history) {
      delete face_history ;
      face_history = 0 ;
    }
#endif
    
    particle_number = result_located_len ;
    last_location_steps = walked_steps;
    // location completed, we need to see whether or not
    // to perform a particle distribution
    ++distribution_tick ;

    bool need_distribution = false ;
    std::string distribution_type ;

    // in case of a parallel ParticleSpace, we need to
    // check to see whether we need to perform a particle distribution
    if(parallel_run) {
      // first we will see if the mean max criterion is met
      // divide first so we don't overflow
      int pmeanp = (particle_number+Loci::MPI_processes/2)/Loci::MPI_processes ;
      // sum of pmeanp over all processors approximates the mean size
      mpi_2int pair_l(pmeanp, particle_number), pair_g ;
      MPI_Op mpi_2int_op ;
      MPI_Op_create( (MPI_User_function*)SUM_MAX_2INT, 1, &mpi_2int_op) ;
      MPI_Allreduce(&pair_l, &pair_g, 1,
                    MPI_2INT, mpi_2int_op, MPI_COMM_WORLD) ;
      MPI_Op_free(&mpi_2int_op) ;
      double mean_size = static_cast<double>(pair_g.first) ;
      int max_size = pair_g.second ;
      if(max_size > 100 && static_cast<double>(max_size) >
         particle_redistribution_mean_max_threshold*mean_size) {
        need_distribution = true ;
        distribution_type = "auto" ;
      }
      // if the mean max criterion is not met, then we see
      // if the pre-set distribution interval is met
      if(!need_distribution) {
        if(distribution_tick >= particle_redistribution_freq) {
          need_distribution = true ;
          distribution_type = "pre" ;
        }
      }
    } // end if(parallel_run)

    if(need_distribution) {
      distribution_tick = 0 ;

      if(Loci::MPI_rank == 0) {
        std::cout << "*****ParticleSpace Message*****" << std::endl ;
        std::cout << "* Distribute particle list ["
                  << distribution_type << "]" ;
        std::cout.flush() ;
      }
      // since the list to_locate points to the
      // passed in list and that it is now empty,
      // we will put the distribution result in it
      std::back_insert_iterator<std::list<ParticleType> >
        result(to_locate) ;

      if(particle_distribution_method == "orb") {
        if(Loci::MPI_rank == 0) {
          std::cout << " (orb) ... " ;
          std::cout.flush() ;
        }
        orb_partition_sequence(located.begin(), located.end(),
                               result, MPI_COMM_WORLD) ;
      } else {
        if(Loci::MPI_rank == 0) {
          std::cout << " (hilbert) ... " ;
          std::cout.flush() ;
        }
        hilbert_partition_sequence(located.begin(), located.end(),
                                   result, MPI_COMM_WORLD);
      }
      
      if(Loci::MPI_rank ==0) {
        std::cout << "Done" << std::endl ;
        std::cout << "*******************************" << std::endl ;
      }
      // we can now get rid of the located list
      std::list<ParticleType>().swap(located) ;
      // we will need to update the mapping info
      new_cells = collect_cells(to_locate.begin(), to_locate.end()) ;
      cell_mapping.update_mapping_info(new_cells) ;
      // update the pull and push comm
      update_comm(new_cells, pull_comm, push_comm) ;
      particle_number = to_locate.size() ;
      
    } else {
      // copy the located list to the to_locate list
      to_locate.splice(to_locate.end(), located) ;
      // don't distribute particles, just update the mapping info
      cell_mapping.update_mapping_info(new_cells) ;
      // updating the pull and push comm
      //update_comm(new_cells, pull_comm, push_comm) ;
      update_comm(new_cells, pull_comm, push_comm) ;
    }

    // we then proceed to freeze the particle to cell map again
    freeze_particles(to_locate) ;

    // we then perform a sorting of the particle list
    // based on their reference cells
    to_locate.sort() ;

    if(!parallel_run)
      return ;
    // finally we need to redistribute bounced_particles and bounce_faces
    std::vector<int> pmap ;
    generate_process_map(factsP->get_init_ptn(), bounce_faces, pmap) ;
    
    distribute_vector(pmap, bounce_particles, MPI_COMM_WORLD) ;
    distribute_vector(pmap, bounce_faces, MPI_COMM_WORLD) ;
    // and remap the bounce_faces to the local numbering
    for(size_t i=0;i!=bounce_faces.size();++i) {
      int& f = bounce_faces[i] ;
      f = cell_mapping.loci_g2l[f] ;
    }
    // also redistribute the iblanked_particles
    std::vector<int> iblanked_cells;
    for(size_t i=0;i!=iblanked_particles.size();++i)
      iblanked_cells.push_back(iblanked_particles[i].get_cell());
    pmap.clear();
    generate_process_map(factsP->get_init_ptn(), iblanked_cells, pmap);
    distribute_vector(pmap, iblanked_particles, MPI_COMM_WORLD);
    // and remap the cell to local number
    for(size_t i=0;i!=iblanked_particles.size();++i) {
      int c = iblanked_particles[i].get_cell();
      c = cell_mapping.loci_g2l[c];
      iblanked_particles[i].set_cell(c);
    }
  }
  
  template<class ParticleType> int ParticleSpace::
  locate_particles_core(std::list<ParticleType>& to_locate,
                        std::list<ParticleType>& located,
                        std::vector<ParticleType>& bounce_particles,
                        std::vector<int>& bounce_faces,
                        std::vector<ParticleType>& iblanked_particles,
                        entitySet& new_cells,
                        size_t& result_located_len)
  {
    entitySet request_cells = cell_mapping.cells;
    int current_walking_step = 0;
    result_located_len = 0;
    // start the location process
    while(true) {
      // first we will to if we are walking too many steps.
      if(current_walking_step > max_walking_step_warning) {
        if(Loci::MPI_rank == 0) {
          std::cout << "*****ParticleSpace Warning*****" << std::endl ;
          std::cout << "* Current walking step > "
                    << max_walking_step_warning << std::endl ;
        }
        int ls = to_locate.size() ;
        int gls = 0 ;
        MPI_Reduce(&ls, &gls, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD) ;
        if(Loci::MPI_rank == 0) {
          std::cout << "* Code may have problems!" << std::endl ;
          std::cout << "* Removing the remaining " << gls
                    << " particles" << std::endl ;
          std::cout << "*******************************" << std::endl ;
        }
        //        typename std::list<ParticleType>::const_iterator ip ;
        //        for(ip=to_locate.begin();ip != to_locate.end();++ip) {
        //          Loci::debugout << "P:pos = " << ip->pos[0] << "::"
        //                         << ip->pos[1] <<std::endl ;
        //        }
        to_locate.clear() ;
      }

      // check whether walking is completed
      if(Loci::GLOBAL_AND(to_locate.empty()))
        break ;

      // first we need to cache the requested cells
      geom.cache(request_cells) ;

      std::pair<LocateListInfo,LocateListInfo> walk_result ;

      if(current_walking_step < face_history_tracing_start) {
        // walk particles without face history checking
        walk_result = move_particles
          (to_locate, located, 
           bounce_particles, bounce_faces, iblanked_particles,
           &ParticleSpace::walk_particle<MeshInteriorFaceTopo,
                                         ParticleType>,
           &ParticleSpace::walk_particle<MeshBoundaryFaceTopo,
                                         ParticleType>) ;
        
      } else {
        // walk particles with face history checking
        if(!face_history)
          face_history = new std::map<int, std::set<int> > ;

        walk_result = move_particles
          (to_locate, located, 
           bounce_particles, bounce_faces, iblanked_particles,
           &ParticleSpace::walk_particle_wht<MeshInteriorFaceTopo,
                                             ParticleType>,
           &ParticleSpace::walk_particle<MeshBoundaryFaceTopo,
                                         ParticleType>) ;
      }

      // finished walking one step
      ++current_walking_step ;
      // modifying the necessary info according the walk_result
      new_cells += walk_result.second.ref_cells ;
      request_cells = walk_result.first.ref_cells ;
      result_located_len += walk_result.second.len ;
    } // end while(true)

    if(face_history) {
      delete face_history ;
      face_history = 0 ;
    }
    return current_walking_step;
  }

  template<class ParticleType, class ForwardIterator> void
  ParticleSpace::register_particles(ForwardIterator first,
                                    ForwardIterator last,
                                    std::list<ParticleType>& rp) {
    // in this function, what we will mainly do is
    // to relabel the particle ids and
    // to add the new mapping information and also
    // freeze the particle to cell mapping in the
    // new list of particles.

    // we will need to do a parallel scan to compute
    // the particle id allocation on each process so
    // that the new ids don't collapse
    long local_id = 0 ;
    long len = std::distance(first, last) ;
    // compute a local id alloc
    MPI_Scan(&len, &local_id, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD) ;
    local_id -= len ;
    local_id += particle_id_alloc ;
    // then set the new alloc for next time
    long tot_len = 0 ;
    MPI_Allreduce(&len, &tot_len, 1,
                  MPI_LONG, MPI_SUM, MPI_COMM_WORLD) ;
    particle_id_alloc += tot_len ;
    // then we begin to renumber the particles
    // However we will need to record the range
    // that we will add particles in list "rp"
    std::list<ParticleType> tmp ;
    for(;first!=last;++first) {
      tmp.push_back(*first) ;
      tmp.back().set_id(local_id++) ;
    }

    // remember the old mapping cells since we want to
    // add our new mappings on top of this.
    // this is a somewhat hack for performance, ideally we
    // would like to do this through a supported API instead
    // of directly using the internal components of the
    // mapping object.
    entitySet current_mapping_cells = cell_mapping.cells;
    // then add the new mapping information
    cell_mapping.add_mapping_info(tmp.begin(), tmp.end()) ;

    // before merging the injected particles, we need to walk them
    // first.  this is because usually the containing cells for the
    // newly injected particles are only an approximate location.
    // we wanted to set the accurate cell locations for the injected
    // particles.  also we wanted to check to make sure that no
    // injected particles are outside of the geometry boundaries.
    // by checking the number of bounced particles, we can report
    // to the user if there are any "out of domain" particles in
    // the injected ones.
    std::list<ParticleType> tmp_located;
    std::vector<ParticleType> bounce_p;
    std::vector<int> bounce_f;
    std::vector<ParticleType> iblanked_p;
    entitySet new_cells;
    size_t result_len;

    locate_particles_core(tmp,tmp_located,
                          bounce_p,bounce_f,iblanked_p,
                          new_cells,result_len);

    // give warning if there are any bounced particles
    // (i.e., particles that are currently outside of domain)
    unsigned int local_b = bounce_p.size();
    unsigned int total_b = 0;
    MPI_Reduce(&local_b, &total_b, 1,
               MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
    if(Loci::MPI_rank == 0) {
      if(total_b > 0) {
        std::cout << "*****ParticleSpace Message*****" << std::endl ;
        std::cout << "* " << total_b << " newly registered particles"
                  << std::endl;
        std::cout << "* are outside of geometry domain and are REMOVED!"
                  << std::endl ;
        std::cout << "*******************************" << std::endl ;
      }
    }
    // give warning also if there are particles inside any iblanked
    // cells (particles that are currently inside a hole in the mesh)
    local_b = iblanked_p.size();
    total_b = 0;
    MPI_Reduce(&local_b, &total_b, 1,
        MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
    if(Loci::MPI_rank == 0) {
      if(total_b > 0) {
        std::cout << "*****ParticleSpace Message*****" << std::endl;
        std::cout << "* " << total_b << " newly registered particles"
          << std::endl;
        std::cout << "* are in holes of the mesh and are REMOVED!"
          << std::endl;
        std::cout << "*******************************" << std::endl;
      }
    }

    // re-update the mapping info.
    cell_mapping.update_mapping_info(new_cells+current_mapping_cells);
    // then update the pull and push communicator
    add_comm(tmp_located.begin(), tmp_located.end(), pull_comm, push_comm) ;
    // finally freeze the particle to cell map,
    // and then modify the particle number accordingly
    particle_number += freeze_particles(tmp_located.begin(),
                                        tmp_located.end()) ;
    // then append tmp to rp
    rp.splice(rp.end(), tmp_located) ;
  }
  
  // sp must be store<list<Particle> >
  template<class ParticleType> void
  ParticleSpace::register_particles(Loci::storeRepP sp,
                                    std::list<ParticleType>& rp) {
    // first reconstruct the store
    store<std::list<ParticleType> > pstore(sp) ;
    // we will then get all the particles inside
    // into a temporary vector
    entitySet dom = pstore.domain() ;
    if(parallel_run) {
      // in case of a parallel run, shrink the domain
      // to exclude those clone regions added by Loci
      Loci::fact_db::distribute_infoP df = factsP->get_distribute_info() ;
      dom &= df->my_entities ;
    }
    std::vector<ParticleType> vp ;
    for(entitySet::const_iterator ei=dom.begin();ei!=dom.end();++ei) {
      vp.insert(vp.end(), pstore[*ei].begin(), pstore[*ei].end()) ;
    }
    // if in a parallel run, we will convert all the cell
    // number to the Loci::fact_db's global numbering
    if(parallel_run) {
      Loci::fact_db::distribute_infoP df = factsP->get_distribute_info() ;
      typename std::vector<ParticleType>::iterator vi ;
      for(vi=vp.begin();vi!=vp.end();++vi) {
        Entity c = vi->get_cell() ;
        vi->set_cell(df->l2g[c]) ;
      }
    }
    register_particles(vp.begin(), vp.end(), rp) ;
  }
  
  
  template<class ParticleType> void
  ParticleSpace::register_particles(std::string restart_file,
                                    std::list<ParticleType>& rp) {
    // first we will check to see if the file exists
    int file_exists = 1 ;
    if(Loci::MPI_rank == 0) {
      struct stat buf ;
      if(stat(restart_file.c_str(), &buf) == -1
         || !S_ISREG(buf.st_mode)) file_exists = 0 ;
    }
    MPI_Bcast(&file_exists, 1, MPI_INT, 0, MPI_COMM_WORLD) ;

    if(file_exists != 1) {
      // return without action, maybe print out a message
      if(Loci::MPI_rank == 0) {
        std::cout << "*****ParticleSpace Message*****" << std::endl ;
        std::cout << "* Error: restart file: "
                  << restart_file << " not found!" << std::endl ;
        std::cout << "* No particles reconstructed." << std::endl ;
        std::cout << "*******************************" << std::endl ;
      }
      return ;
    }
    
    // then read the file
    hid_t file_id = Loci::hdf5OpenFile(restart_file.c_str(),
                                       H5F_ACC_RDONLY, H5P_DEFAULT) ;
    int read_success = 1 ;
    if(Loci::MPI_rank == 0) {
      if(file_id < 0)
        read_success = 0 ;
    }
    MPI_Bcast(&read_success, 1, MPI_INT, 0, MPI_COMM_WORLD) ;
    if(read_success == 0) {
      if(Loci::MPI_rank == 0) {
        std::cout << "*****ParticleSpace Message*****" << std::endl ;
        std::cout << "* Error: restart file: "
                  << restart_file << " read error!" << std::endl ;
        std::cout << "* No particles reconstructed." << std::endl ;
        std::cout << "*******************************" << std::endl ;
      }
      return ;
    }
    
    std::vector<ParticleType> buf ;
    readUnorderedVector(file_id, "particle restart", buf) ;
    
    Loci::hdf5CloseFile(file_id) ;

    // an important thing to do is to set the global particle
    // id allocation correctly, otherwise newly generated particles
    // will have ids that clash with older ones
    typename std::vector<ParticleType>::iterator vpi ;
    long local_max_id = 0 ;
    for(vpi=buf.begin();vpi!=buf.end();++vpi) {
      if(vpi->get_id() > local_max_id)
        local_max_id = vpi->get_id() ;
    }
    long global_max_id = 0 ;
    MPI_Allreduce(&local_max_id, &global_max_id,
                  1, MPI_LONG, MPI_MAX, MPI_COMM_WORLD) ;
    particle_id_alloc = global_max_id + 1 ;

    // we will first redistribute particles according to
    // the ORB partition, this may help to reduce the memory
    // consumption in the parallel nearest neighbor search
    // performed in the following step
    std::vector<ParticleType> orb_buf ;
    std::back_insert_iterator<std::vector<ParticleType> >
      result(orb_buf) ;
    
    if(particle_distribution_method == "orb") {
      orb_partition_sequence(buf.begin(),
                             buf.end(), result, MPI_COMM_WORLD) ;
    } else {
      hilbert_partition_sequence(buf.begin(), buf.end(),
                                 result, MPI_COMM_WORLD);
    }
    
    // destroy buf now
    std::vector<ParticleType>().swap(buf) ;
    // get all the cells
    constraint geom_cells = factsP->get_variable("geom_cells") ;
    entitySet cells = *geom_cells ;
    store<vec3d> cellcenter(factsP->get_variable("cellcenter")) ;
    
    std::vector<Loci::kdTree::coord3d> cell_pts(cells.size()) ;
    std::vector<int> cell_ids(cells.size()) ;
    int count = 0 ;

    if(parallel_run) {
      Loci::fact_db::distribute_infoP df = factsP->get_distribute_info() ;
      // in parallel run, we need to put the global cell number in
      for(entitySet::const_iterator ei=cells.begin();
          ei!=cells.end();++ei,++count) {
        const vec3d& cc = cellcenter[*ei] ;
        cell_pts[count][0] = cc.x ;
        cell_pts[count][1] = cc.y ;
        cell_pts[count][2] = cc.z ;
        cell_ids[count]    = df->l2g[*ei] ;
      }
    } else {
      for(entitySet::const_iterator ei=cells.begin();
          ei!=cells.end();++ei,++count) {
        const vec3d& cc = cellcenter[*ei] ;
        cell_pts[count][0] = cc.x ;
        cell_pts[count][1] = cc.y ;
        cell_pts[count][2] = cc.z ;
        cell_ids[count]    = *ei ;
      }
    } // end if(parallel_run)

    std::vector<Loci::kdTree::coord3d> search_pts(orb_buf.size()) ;
    
    for(size_t i=0;i<orb_buf.size();++i) {
      vec3d pp = orb_buf[i].get_position() ;
      search_pts[i][0] = pp.x ;
      search_pts[i][1] = pp.y ;
      search_pts[i][2] = pp.z ;
    }
    
    std::vector<int> closest(orb_buf.size(), -1) ;
//#define USE_APPROX_NN
#ifdef USE_APPROX_NN
    approx_NN(cell_pts, cell_ids, search_pts, closest, MPI_COMM_WORLD) ;
#else
    Loci::parallelNearestNeighbors(cell_pts, cell_ids,
                                   search_pts, closest, MPI_COMM_WORLD) ;
#endif

    // reset the reference cell id
    for(size_t i=0;i<orb_buf.size();++i) {
      orb_buf[i].set_cell(closest[i]) ;
    }
    // update the mapping info
    cell_mapping.update_mapping_info(orb_buf.begin(), orb_buf.end()) ;
    particle_number += freeze_particles(orb_buf.begin(), orb_buf.end()) ;
    
    // create a list from orb_buf
    std::list<ParticleType> tmp(orb_buf.begin(),orb_buf.end()) ;
    // then destroy orb_buf since it is no longer needed
    std::vector<ParticleType>().swap(orb_buf) ;
    
    // then we just locate the particles
    locate_particles(tmp) ;
    
    // then append tmp to the list "rp"
    rp.splice(rp.end(), tmp) ;
  }

  template<class ParticleType> void
  ParticleSpace::insert_new_particles(std::vector<ParticleType> &buf,
                                    std::list<ParticleType>& rp) {
    
    
    // an important thing to do is to set the global particle
    // id allocation correctly, otherwise newly generated particles
    // will have ids that clash with older ones
    typename std::vector<ParticleType>::iterator vpi ;
    long local_max_id = 0 ;
    for(vpi=buf.begin();vpi!=buf.end();++vpi) {
      if(vpi->get_id() > local_max_id)
        local_max_id = vpi->get_id() ;
    }
    long global_max_id = 0 ;
    MPI_Allreduce(&local_max_id, &global_max_id,
                  1, MPI_LONG, MPI_MAX, MPI_COMM_WORLD) ;
    particle_id_alloc = global_max_id + 1 ;

    // we will first redistribute particles according to
    // the ORB partition, this may help to reduce the memory
    // consumption in the parallel nearest neighbor search
    // performed in the following step
    std::vector<ParticleType> orb_buf ;
    std::back_insert_iterator<std::vector<ParticleType> >
      result(orb_buf) ;

    if(particle_distribution_method == "orb") {
      orb_partition_sequence(buf.begin(),
                             buf.end(), result, MPI_COMM_WORLD) ;
    } else {
      hilbert_partition_sequence(buf.begin(), buf.end(),
                                 result, MPI_COMM_WORLD);
    }

    // destroy buf now
    std::vector<ParticleType>().swap(buf) ;
    // get all the cells
    constraint geom_cells = factsP->get_variable("geom_cells") ;
    entitySet cells = *geom_cells ;
    store<vec3d> cellcenter(factsP->get_variable("cellcenter")) ;
    
    std::vector<Loci::kdTree::coord3d> cell_pts(cells.size()) ;
    std::vector<int> cell_ids(cells.size()) ;
    int count = 0 ;

    if(parallel_run) {
      Loci::fact_db::distribute_infoP df = factsP->get_distribute_info() ;
      // in parallel run, we need to put the global cell number in
      for(entitySet::const_iterator ei=cells.begin();
          ei!=cells.end();++ei,++count) {
        const vec3d& cc = cellcenter[*ei] ;
        cell_pts[count][0] = cc.x ;
        cell_pts[count][1] = cc.y ;
        cell_pts[count][2] = cc.z ;
        cell_ids[count]    = df->l2g[*ei] ;
      }
    } else {
      for(entitySet::const_iterator ei=cells.begin();
          ei!=cells.end();++ei,++count) {
        const vec3d& cc = cellcenter[*ei] ;
        cell_pts[count][0] = cc.x ;
        cell_pts[count][1] = cc.y ;
        cell_pts[count][2] = cc.z ;
        cell_ids[count]    = *ei ;
      }
    } // end if(parallel_run)

    std::vector<Loci::kdTree::coord3d> search_pts(orb_buf.size()) ;
    for(size_t i=0;i<orb_buf.size();++i) {
      vec3d pp = orb_buf[i].get_position() ;
      search_pts[i][0] = pp.x ;
      search_pts[i][1] = pp.y ;
      search_pts[i][2] = pp.z ;
    }

    std::vector<int> closest(orb_buf.size(), -1) ;
#ifdef USE_APPROX_NN
    approx_NN(cell_pts, cell_ids, search_pts, closest, MPI_COMM_WORLD) ;
#else
    Loci::parallelNearestNeighbors(cell_pts, cell_ids,
                                   search_pts, closest, MPI_COMM_WORLD) ;
#endif
    
    // reset the reference cell id
    for(size_t i=0;i<orb_buf.size();++i) {
      orb_buf[i].set_cell(closest[i]) ;
    }

    // update the mapping info
    cell_mapping.update_mapping_info(orb_buf.begin(), orb_buf.end()) ;
    particle_number += freeze_particles(orb_buf.begin(), orb_buf.end()) ;
    
    // create a list from orb_buf
    std::list<ParticleType> tmp(orb_buf.begin(),orb_buf.end()) ;
    // then destroy orb_buf since it is no longer needed
    std::vector<ParticleType>().swap(orb_buf) ;

    // then we just locate the particles
    locate_particles(tmp) ;

    // then append tmp to the list "rp"
    rp.splice(rp.end(), tmp) ;
  }

  //////////////////////////////////////////////////////////////////
  // This section provides the prototypes for communication       //
  // routines that are based on the ParticleSpace type            //
  //////////////////////////////////////////////////////////////////
  typedef ParticleSpace::Communicator::P2pCommInfo P2pCommInfo ;

#ifdef TO_BE_REVISED
  // These two methods are to be revised
  
  // this is another version of distributed_map_image that uses
  // point 2 point commnunication instead of the collective ones
  entitySet
  distributed_map_image_p2p(const entitySet& domain,
                            Loci::MapRepP m, Loci::fact_db* factsP) ;

  // this is a point 2 point version of "distributed_map_image"
  // that uses a supplied communication structure
  // NOTE: this version assumes that all the necessary communication
  // links are presented in the send/recv records and will only
  // compute the results based on the send/recv records. i.e.,
  // if "domain" includes a process's local entitySet, then
  // the local partition needs to be put in the send/recv records also
  entitySet
  distributed_map_image_p2p(const entitySet& domain,
                            Loci::MapRepP m, Loci::fact_db* factsP,
                            const std::vector<P2pCommInfo>& send,
                            const std::vector<P2pCommInfo>& recv) ;
#endif
  
  // this is an optimized version that is used to expand
  // the domain of store and storeVec, and all other stores
  // (but excluding Map, multiMap and multiStore)
  // it is similar to the above versions except that it
  // uses a supplied communication structure and performs
  // point 2 point communication instead of collective ones.
  // the communication structure also contains the domain
  // to send and receive.
  // this version also expands a vector of stores that
  // share the same domain for efficiency. it returns a
  // vector of storeRepP for the expanded domain.
  // The pack_remap is used to remap the domain entities
  // to a new numbering when packing for communication.
  // typically, it is a local -> global number remapping.
  // The unpack_remap is used to remap the domain entities
  // when receiving the messages. it is usually a
  // global -> local number remapping
  // Since the returned stores only have
  // the expanded domains contained, we would name the function
  // as "get_remote_stores" for clarity.

  // NOTE:
  // Ideally, we would also like to include Map being able
  // to be expanded using this function. this way, if Maps
  // share the same expansion domain as any stores, then
  // we can save some messages (both number of messages
  // and the size of the messages). The current interfaces
  // are enought for doing this. But sadly, the pack/unpack
  // methods in the Map class only blindly pack the image,
  // what we really want is that they be able to use
  // different numberings when packing/unpacking.
  // We may want to come back in the future to revise the
  // Loci storage containers implementation. Since the present
  // Loci only practically use store and storeVec, only methods
  // in those two classes are well designed and tested.
  // Map and multiMap are never computed and communicated after
  // the program starts, so some of the routines inside may
  // not be sufficient and are not well thought of. multiStore
  // seems never get used and I suspect that its implementation
  // is not robust also.
  void
  get_remote_stores(std::vector<Loci::storeRepP>& in,
                    const std::vector<P2pCommInfo>& send, // the comm
                    const std::vector<P2pCommInfo>& recv, // struct
                    const Map& pack_remap,
                    const dMap& unpack_remap,
                    std::vector<Loci::storeRepP>& out) ;

  // this is the new interface that uses the new ParticleSpace type
  void
  get_remote_stores(std::vector<Loci::storeRepP>& in,
                    std::vector<Loci::storeRepP>& out,
                    const ParticleSpace::Communicator& comm) ;

  // these are versions that are for a single rep
  Loci::storeRepP
  get_remote_stores(Loci::storeRepP in,
                    const std::vector<P2pCommInfo>& send,
                    const std::vector<P2pCommInfo>& recv,
                    const Map& pack_remap, const dMap& unpack_remap) ;

  // the new interface that uses ParticleSpace type
  Loci::storeRepP
  get_remote_stores(Loci::storeRepP in,
                    const ParticleSpace::Communicator& comm) ;

  // this is a version specifically used in the incremental
  // caching code, the storeRepPs in the out vector are
  // assumed to have been allocated with enough space to
  // receive the remote contents
  void
  get_remote_stores_inc(std::vector<Loci::storeRepP>& in,
                        const std::vector<P2pCommInfo>& send,
                        const std::vector<P2pCommInfo>& recv,
                        const Map& pack_remap,
                        const dMap& unpack_remap,
                        std::vector<Loci::storeRepP>& out) ;

  // the new interface that uses ParticleSpace type
  void
  get_remote_stores_inc(std::vector<Loci::storeRepP>& in,
                        std::vector<Loci::storeRepP>& out,
                        const ParticleSpace::Communicator& comm) ;

  // this method is mostly the same as the "get_remote_stores" one,
  // just that this is a specific version applied to maps in Loci.
  // In commnicuating maps, we will need to have more remapping
  // options, specifically the domain the map can be remapped before/afer
  // the communication, also the image of the map can be remapped as well.
  // In addition, this method also works for dynamic maps (dMap) where
  // such remapping might not be desirable.  In this case, if the
  // remapping info are NULL pointers, then we just don't do remapping
  void
  get_remote_maps(std::vector<Loci::storeRepP>& in,
                  const std::vector<P2pCommInfo>& send, // the comm
                  const std::vector<P2pCommInfo>& recv, // struct
                  Map* dom_pack, // remapping info for packing domains
                  dMap* dom_unpack, // remapping info for unpacking domains
                  Map* img_pack, // remapping info for packing images
                  dMap* img_unpack, // remapping infor for unpacking images
                  std::vector<Loci::storeRepP>& out);

#ifdef TO_BE_REVISED
  // These methods are to be revised
  // this version is for Map expansion
  // only the domain is remapped, images are NOT
  Loci::storeRepP
  get_remote_Map(const Map& in,
                 const std::vector<P2pCommInfo>& send,
                 const std::vector<P2pCommInfo>& recv,
                 const dMap& remap, Loci::fact_db* factsP) ;

  // interface that uses ParticleSpace type
  Loci::storeRepP
  get_remote_Map(const Map& in,
                 const ParticleSpace::Communicator& comm,
                 Loci::fact_db* factsP) ;
  
  Loci::storeRepP
  get_remote_Map(const Map& in,
                 const entitySet& domain,
                 const dMap& remap, Loci::fact_db* factsP) ;

  // this version is for multiMap expansion
  // only the domain is remapped, images are NOT
  Loci::storeRepP
  get_remote_multiMap(const Loci::multiMap& in,
                      const std::vector<P2pCommInfo>& send,
                      const std::vector<P2pCommInfo>& recv,
                      const dMap& remap, Loci::fact_db* factsP) ;
  
  Loci::storeRepP
  get_remote_multiMap(const Loci::multiMap& in,
                      const entitySet& domain,
                      const dMap& remap, Loci::fact_db* factsP) ;
#endif
  
  // this function will perform a reduction on the passed in
  // vector of stores, the "in_rep" vector.
  // The passed in vector of stores are assumed to be in a
  // process's local numbering (could be static or dynamic
  // local numbering). they typically contain domains that
  // belong to other processes' inititial domain partition.
  // this function will send/recv the contents of stores
  // based on the supplied send/recv communication structure,
  // it will construct the received stores on each process
  // in the "out_rep" vector. corresponding to the order
  // specified in "in_rep". since the contents of a domain
  // element in the receiving stores may have contributions
  // from multiple processes, this function performs a
  // reduction on such domain entities based on the supplied
  // reduction operation defined in the "join_op" vector.
  // the "unit_op" vector is used for setting the desired
  // unit value for the received stores before reduction.
  // the "pack_map" is used to remap the domain of the
  // sending stores in "in_rep" before sending. it is usually
  // a local->global number remap. the "unpack_map" is used
  // to remap the domain of the received stores in "out_rep"
  // after receiving. it is usually a global->local number
  // remap.
  void
  reduce_remote_stores(std::vector<Loci::storeRepP>& in_rep,
                       std::vector<store_traverserP>& unit_op,
                       std::vector<loci_joinerP>& join_op,
                       const std::vector<P2pCommInfo>& send,
                       const std::vector<P2pCommInfo>& recv,
                       const Map& pack_remap,
                       const dMap& unpack_remap,
                       std::vector<Loci::storeRepP>& out_rep) ;

  // an interface that uses the ParticleSpace::Communicator
  void
  reduce_remote_stores(std::vector<Loci::storeRepP>& in_rep,
                       std::vector<store_traverserP>& unit_op,
                       std::vector<loci_joinerP>& join_op,
                       std::vector<Loci::storeRepP>& out_rep,
                       const ParticleSpace::Communicator& comm) ;

  // different version of reduce_remote_stores
  Loci::storeRepP
  reduce_remote_stores(Loci::storeRepP in_rep,
                       store_traverserP unit_op,
                       loci_joinerP join_op,
                       const std::vector<P2pCommInfo>& send,
                       const std::vector<P2pCommInfo>& recv,
                       const Map& pack_remap,
                       const dMap& unpack_remap) ;

  Loci::storeRepP
  reduce_remote_stores(Loci::storeRepP in_rep,
                       store_traverserP unit_op,
                       loci_joinerP join_op,
                       const ParticleSpace::Communicator& comm) ;

  // a little structure used in the restart code
  struct F2G {
    int f, g;                   // f is file number, g is global number
    F2G(int a=0, int b=0):f(a),g(b) {}
    // mandatory methods for parallel sorting
    size_t
    pack_size() const { return 2*sizeof(int); }

    void
    pack(unsigned char* buffer, size_t& position) const
    {
      int* b = reinterpret_cast<int*>(buffer+position);
      *b = f; ++b; *b=g;
      position += pack_size();
    }

    void
    unpack(unsigned char* buffer, size_t& position)
    {
      int* b = reinterpret_cast<int*>(buffer+position);
      f=*b; ++b; g=*b;
      position += pack_size();
    }
  };
  inline bool
  operator<(const F2G& l, const F2G& r)
  { return l.f < r.f; }
  inline bool
  lessF2G(const F2G& l, const F2G& r)
  { return l.f < r.f; }
  
  // this is a new restart procedure that is based on file numbers.
  // all the particles read in have their cells in the file number.
  // we therefore need to convert all these file numbers into the
  // corresponding global numbers.  one benefit of this approach is
  // that there is no need to relocate the particles, in particular,
  // there is no need to use nearest neighbor prior to the relocation,
  // which saves time and space costs.
  // to convert the file number to the global number, we collect all
  // file number distributions on all processes and then build a new
  // communication structure and then communicate the necessary
  // file number to global number map to each process.
  template<class ParticleType> void
  ParticleSpace::restart_particles(std::string restart_file,
                                   std::list<ParticleType>& rp) {
    // first we will check to see if the file exists
    int file_exists = 1 ;
    if(Loci::MPI_rank == 0) {
      struct stat buf ;
      if(stat(restart_file.c_str(), &buf) == -1
         || !S_ISREG(buf.st_mode)) file_exists = 0 ;
    }
    MPI_Bcast(&file_exists, 1, MPI_INT, 0, MPI_COMM_WORLD) ;

    if(file_exists != 1) {
      // return without action, maybe print out a message
      if(Loci::MPI_rank == 0) {
        std::cout << "*****ParticleSpace Message*****" << std::endl ;
        std::cout << "* Error: restart file: "
                  << restart_file << " not found!" << std::endl ;
        std::cout << "* No particles reconstructed." << std::endl ;
        std::cout << "*******************************" << std::endl ;
      }
      return ;
    }
    
    // then read the file
    hid_t file_id = Loci::hdf5OpenFile(restart_file.c_str(),
                                       H5F_ACC_RDONLY, H5P_DEFAULT) ;
    int read_success = 1 ;
    if(Loci::MPI_rank == 0) {
      if(file_id < 0)
        read_success = 0 ;
    }
    MPI_Bcast(&read_success, 1, MPI_INT, 0, MPI_COMM_WORLD) ;
    if(read_success == 0) {
      if(Loci::MPI_rank == 0) {
        std::cout << "*****ParticleSpace Message*****" << std::endl ;
        std::cout << "* Error: restart file: "
                  << restart_file << " read error!" << std::endl ;
        std::cout << "* No particles reconstructed." << std::endl ;
        std::cout << "*******************************" << std::endl ;
      }
      return ;
    }
    
    std::vector<ParticleType> buf ;
    readUnorderedVector(file_id, "particle restart", buf) ;
    
    Loci::hdf5CloseFile(file_id) ;

    // an important thing to do is to set the global particle
    // id allocation correctly, otherwise newly generated particles
    // will have ids that clash with older ones
    typename std::vector<ParticleType>::iterator vpi ;
    long local_max_id = 0 ;
    for(vpi=buf.begin();vpi!=buf.end();++vpi) {
      if(vpi->get_id() > local_max_id)
        local_max_id = vpi->get_id() ;
    }
    long global_max_id = 0 ;
    MPI_Allreduce(&local_max_id, &global_max_id,
                  1, MPI_LONG, MPI_MAX, MPI_COMM_WORLD) ;
    particle_id_alloc = global_max_id + 1 ;

    // we will first redistribute particles
    std::vector<ParticleType> orb_buf ;
    std::back_insert_iterator<std::vector<ParticleType> >
      result(orb_buf) ;

    if(particle_distribution_method == "orb") {
      orb_partition_sequence(buf.begin(),
                             buf.end(), result, MPI_COMM_WORLD) ;
    } else {
      hilbert_partition_sequence(buf.begin(), buf.end(),
                                 result, MPI_COMM_WORLD);
    }

    // destroy buf now
    std::vector<ParticleType>().swap(buf) ;
    
    // now we need to process the file number and convert it to
    // global numbers, we only do this for a parallel run.
    if(parallel_run) {
      // get all the cells
      Loci::fact_db::distribute_infoP df = factsP->get_distribute_info() ;
      constraint geom_cells = factsP->get_variable("geom_cells") ;
      entitySet cells = *geom_cells ;
      // file number does not exist for clone region entities
      cells &= df->my_entities ;
      // first of all, we only have a global to file number mapping,
      // we need to build a file to global number mapping
      // we build a vector of pairs instead of a Map so that we can sort them
      Map l2g ;  l2g = df->l2g.Rep();
      Loci::dMap g2f;  g2f = df->g2f.Rep();
      std::vector<F2G> f2g_local;
      for(entitySet::const_iterator ei=cells.begin();ei!=cells.end();++ei) {
        Loci::Entity g = l2g[*ei];
        Loci::Entity f = g2f[g];
        f2g_local.push_back(F2G(f,g));
      }
      // then sort and balance the f2g vector
      Loci::parSampleSort(f2g_local, lessF2G, MPI_COMM_WORLD);
      std::vector<F2G> f2gv;
      std::back_insert_iterator<std::vector<F2G> > bii(f2gv);
      balance_sequence(f2g_local.begin(), f2g_local.end(),
                       bii, MPI_COMM_WORLD);
      std::vector<F2G>().swap(f2g_local);
      // then communicate the range of file numbers
      int frange_local[2];
      frange_local[0] = f2gv[0].f;
      frange_local[1] = f2gv[f2gv.size()-1].f;
      std::vector<int> frange(2*Loci::MPI_processes);
      MPI_Allgather(frange_local, 2, MPI_INT,
                    &frange[0], 2, MPI_INT, MPI_COMM_WORLD);
      // gather all cells from local particles
      entitySet cells_f;
      for(size_t i=0;i!=orb_buf.size();++i)
        cells_f += orb_buf[i].get_cell();
      // compute the distribution
      std::vector<Loci::entitySet> dist(Loci::MPI_processes);
      for(int p=0;p!=Loci::MPI_processes;++p)
        dist[p] = cells_f & interval(frange[2*p],frange[2*p+1]);
      // transpose the distribution
      std::vector<Loci::entitySet> dist_t(Loci::MPI_processes);
      transpose_vector_entitySet_opt(dist, dist_t);
      // build the communication structure
      std::vector<P2pCommInfo> send, recv;
      for(size_t i=0;i!=dist.size();++i) {
        const entitySet& es = dist[i];
        if(es != EMPTY) {
          // we don't need local number in this case since we
          // aere just communicating the global number and
          // file numbers in a dynamic map
          recv.push_back(P2pCommInfo(i,es,es));
        }
      }
      for(size_t i=0;i!=dist_t.size();++i) {
        const entitySet& es = dist_t[i];
        if(es != EMPTY)
          send.push_back(P2pCommInfo(i,es,es));
      }
      // to communicate, convert the local vector of pair to a dmap
      Loci::dMap f2g;
      for(size_t i=0;i!=f2gv.size();++i)
        f2g[f2gv[i].f] = f2gv[i].g;
      // we don't need the f2gv vector any more
      std::vector<F2G>().swap(f2gv);
      
      // now we are ready to send/recv the f2g maps
      std::vector<Loci::storeRepP> in(1), out(1);
      in[0] = f2g.Rep();
      // no need to remap any of the domain and images
      // just do a straight send/recv
      get_remote_maps(in, send, recv, 0, 0, 0, 0, out);
      Loci::dMap cells_f2g(out[0]);
      // some checks
      Loci::entitySet recv_f = cells_f2g.domain();
      if(Loci::GLOBAL_OR(cells_f != recv_f)) {
        if(Loci::MPI_rank == 0) {
          std::cout << "*****ParticleSpace Message*****" << std::endl;
          std::cout << "* Error: restart particles failed! " << std::endl;
          std::cout << "* File number communication failed!" << std::endl;
          std::cout << "* No particles reconstructed." << std::endl ;
          std::cout << "*******************************" << std::endl;
        }
        return;
      }
      // renumber the cells to global number in all particles
      for(size_t i=0;i!=orb_buf.size();++i) {
        Loci::Entity c = orb_buf[i].get_cell();
        orb_buf[i].set_cell(cells_f2g[c]);
      }
    } // end of if(parallel_run)

    // update various particle space info
    entitySet new_cells = collect_cells(orb_buf.begin(), orb_buf.end());
    cell_mapping.update_mapping_info(new_cells) ;
    update_comm(new_cells, pull_comm, push_comm);
    
    particle_number += freeze_particles(orb_buf.begin(), orb_buf.end()) ;
    // then append orb_buf to the list "rp"
    rp.insert(rp.end(), orb_buf.begin(), orb_buf.end()) ;
  }

  template<class ParticleType> void
  ParticleSpace::restart_particles_slow(std::string restart_file,
                                        std::list<ParticleType>& rp) {
    // first we will check to see if the file exists
    int file_exists = 1 ;
    if(Loci::MPI_rank == 0) {
      struct stat buf ;
      if(stat(restart_file.c_str(), &buf) == -1
         || !S_ISREG(buf.st_mode)) file_exists = 0 ;
    }
    MPI_Bcast(&file_exists, 1, MPI_INT, 0, MPI_COMM_WORLD) ;

    if(file_exists != 1) {
      // return without action, maybe print out a message
      if(Loci::MPI_rank == 0) {
        std::cout << "*****ParticleSpace Message*****" << std::endl ;
        std::cout << "* Error: restart file: "
                  << restart_file << " not found!" << std::endl ;
        std::cout << "* No particles reconstructed." << std::endl ;
        std::cout << "*******************************" << std::endl ;
      }
      return ;
    }
    
    // then read the file
    hid_t file_id = Loci::hdf5OpenFile(restart_file.c_str(),
                                       H5F_ACC_RDONLY, H5P_DEFAULT) ;
    int read_success = 1 ;
    if(Loci::MPI_rank == 0) {
      if(file_id < 0)
        read_success = 0 ;
    }
    MPI_Bcast(&read_success, 1, MPI_INT, 0, MPI_COMM_WORLD) ;
    if(read_success == 0) {
      if(Loci::MPI_rank == 0) {
        std::cout << "*****ParticleSpace Message*****" << std::endl ;
        std::cout << "* Error: restart file: "
                  << restart_file << " read error!" << std::endl ;
        std::cout << "* No particles reconstructed." << std::endl ;
        std::cout << "*******************************" << std::endl ;
      }
      return ;
    }
    
    std::vector<ParticleType> buf ;
    readUnorderedVector(file_id, "particle restart", buf) ;
    
    Loci::hdf5CloseFile(file_id) ;

    // an important thing to do is to set the global particle
    // id allocation correctly, otherwise newly generated particles
    // will have ids that clash with older ones
    typename std::vector<ParticleType>::iterator vpi ;
    long local_max_id = 0 ;
    for(vpi=buf.begin();vpi!=buf.end();++vpi) {
      if(vpi->get_id() > local_max_id)
        local_max_id = vpi->get_id() ;
    }
    long global_max_id = 0 ;
    MPI_Allreduce(&local_max_id, &global_max_id,
                  1, MPI_LONG, MPI_MAX, MPI_COMM_WORLD) ;
    particle_id_alloc = global_max_id + 1 ;

    // we will first redistribute particles according to the ORB partition
    std::vector<ParticleType> orb_buf ;
    std::back_insert_iterator<std::vector<ParticleType> >
      result(orb_buf) ;
    
    if(particle_distribution_method == "orb") {
      orb_partition_sequence(buf.begin(),
                             buf.end(), result, MPI_COMM_WORLD) ;
    } else {
      hilbert_partition_sequence(buf.begin(), buf.end(),
                                 result, MPI_COMM_WORLD);
    }

    // destroy buf now
    std::vector<ParticleType>().swap(buf) ;
    // now we need to process the file number and convert it to
    // global numbers, we only do this for a parallel run.
    if(parallel_run) {
      // get all the cells
      Loci::fact_db::distribute_infoP df = factsP->get_distribute_info() ;
      constraint geom_cells = factsP->get_variable("geom_cells") ;
      entitySet cells = *geom_cells ;
      cells &= df->my_entities ;
      // first of all, we only have a global to file number mapping,
      // we need to build a file to global number mapping
      Map l2g ;  l2g = df->l2g.Rep();
      Loci::dMap g2f;  g2f = df->g2f.Rep();
      Loci::dMap f2g;
      Loci::entitySet local_f;  // all file numbers on this process
      for(entitySet::const_iterator ei=cells.begin();ei!=cells.end();++ei) {
        Loci::Entity g = l2g[*ei];
        Loci::Entity f = g2f[g];
        f2g[f] = g;
        local_f += f;
      }
      // then we need to gather the file number distribution
      std::vector<entitySet> f_ptn = gather_all_entitySet(local_f);
      // gather all cells from local particles
      entitySet cells_f;
      for(size_t i=0;i!=orb_buf.size();++i)
        cells_f += orb_buf[i].get_cell();
      // compute the distribution
      std::vector<Loci::entitySet> dist(f_ptn.size());
      for(size_t i=0;i!=f_ptn.size();++i)
        dist[i] = cells_f & f_ptn[i];
      // transpose the distribution
      std::vector<Loci::entitySet> dist_t(f_ptn.size());
      transpose_vector_entitySet_opt(dist, dist_t);
      // build the communication structure
      std::vector<P2pCommInfo> send, recv;
      for(size_t i=0;i!=dist.size();++i) {
        const entitySet& es = dist[i];
        if(es != EMPTY) {
          // we don't need local number in this case since we
          // aere just communicating the global number and
          // file numbers in a dynamic map
          recv.push_back(P2pCommInfo(i,es,es));
        }
      }
      for(size_t i=0;i!=dist_t.size();++i) {
        const entitySet& es = dist_t[i];
        if(es != EMPTY)
          send.push_back(P2pCommInfo(i,es,es));
      }
      // now we are ready to send/recv the f2g maps
      std::vector<Loci::storeRepP> in(1), out(1);
      in[0] = f2g.Rep();
      // no need to remap any of the domain and images
      // just do a straight send/recv
      get_remote_maps(in, send, recv, 0, 0, 0, 0, out);
      Loci::dMap cells_f2g(out[0]);
      // some checks
      Loci::entitySet recv_f = cells_f2g.domain();
      if(Loci::GLOBAL_OR(cells_f != recv_f)) {
        if(Loci::MPI_rank == 0) {
          std::cout << "*****ParticleSpace Message*****" << std::endl;
          std::cout << "* Error: restart particles failed! " << std::endl;
          std::cout << "* File number communication failed!" << std::endl;
          std::cout << "* No particles reconstructed." << std::endl ;
          std::cout << "*******************************" << std::endl;
        }
        return;
      }
      // renumber the cells to global number in all particles
      for(size_t i=0;i!=orb_buf.size();++i) {
        Loci::Entity c = orb_buf[i].get_cell();
        orb_buf[i].set_cell(cells_f2g[c]);
      }
    } // end of if(parallel_run)

    // update various particle space info
    entitySet new_cells = collect_cells(orb_buf.begin(), orb_buf.end());
    cell_mapping.update_mapping_info(new_cells) ;
    update_comm(new_cells, pull_comm, push_comm);
    
    particle_number += freeze_particles(orb_buf.begin(), orb_buf.end()) ;
    // then append orb_buf to the list "rp"
    rp.insert(rp.end(), orb_buf.begin(), orb_buf.end()) ;
  }

} // end of namespace lagrangianP

#endif

// ... and end of the big mess ...
