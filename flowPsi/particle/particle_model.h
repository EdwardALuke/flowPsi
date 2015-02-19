#ifndef PARTICLE_MODEL_H
#define PARTICLE_MODEL_H
#include <Loci>

#include "particle_config.h"
#include "particle_base.h"
#include <vector>
#include <map>
#include <string>

namespace lagrangianP {

  class momentumIntegrator ;
  class dropletBreakupMethod ;
  using Loci::CPTR ;
  
  enum BreakUpModel {NO_BREAKUP, SIMPLE_WEBER} ;
  enum DragModel {CLIFT_GAUVIN,DYNAMIC_DROPLET,LOTH} ;

  void extract_species(const Loci::options_list &ol,
		       std::vector<std::string> &species_names,
		       std::vector<double> &species_mf) ;

  struct particleModel {
    std::vector<Loci::options_list> binDB ;
    std::map<std::string,Loci::options_list> speciesDB ;
    struct binInfo {
      string name ;
    } ;

    std::vector<binInfo> bins ;

    // Coupling Model
    std::string momentumCouplingModel ;
    DragModel momentumCouplingMethod ;
  public:
    particleModel() {}

    friend std::ostream &operator<<(std::ostream &s, const particleModel &pm) ;
    friend std::istream &operator>>(std::istream &s, particleModel &pm) ;
  } ;
  std::ostream &operator<<(std::ostream &s, const particleModel &pm) ;
  std::istream &operator>>(std::istream &s, particleModel &pm) ;


  class densityFunction {
    double specificVolume ;
  public:
    densityFunction() {
      specificVolume = -1 ;
    }
    densityFunction(double density) {
      specificVolume = 1./density ;
    }
    densityFunction(const Loci::options_list &ol,
		    const std::map<std::string,Loci::options_list> &speciesDB) ;
    double get_density() const {
      return 1./specificVolume;
    }
  } ;

  class tensionFunction {
    bool is_defined ;
    double tension ;
  public:
    tensionFunction() { tension = 1e30 ; is_defined = false; }
    tensionFunction(double tin) : tension(tin) {}
    tensionFunction(const Loci::options_list &ol,
		    const std::map<std::string,Loci::options_list> &speciesDB) ;
    bool valid() const { return is_defined; }
    double getTension() const {
      return tension ;
    }
  } ;
    
  struct particleBinEoS {
    // mixture fractions
    int ns ;
    string name ;
    Loci::Array<double,2> massFractions ;

    densityFunction rhop ; // particle density calculator
    tensionFunction tensionp ; // particle surface tension
    double COR ; // Coeff of Restitution

    CPTR<momentumIntegrator> momentumUpdater ;
    CPTR<dropletBreakupMethod> breakupUpdater ;

    void initialize(const particleModel &pm, string name,
		    const Loci::options_list &ol,
		    const std::map<std::string,Loci::options_list> &speciesDB) ;
  } ;

  typedef std::vector<lagrangianP::particleBinEoS> particleEoSList;

  class AuxiliaryInfoDatabase {
    struct data_info {
      int size ;
      int id ; 
    } ;
    std::map<string,data_info> DataMap ;
    int AuxSize ;
  public:
    AuxiliaryInfoDatabase() {
      AuxSize = 0 ;
    }
    int insertItem(string name, int size) {
      std::map<string,data_info>::const_iterator ii ;
      ii = DataMap.find(name) ;
      if(ii == DataMap.end()) {
        // not in database so add it
        data_info di ;
        di.size = size ;
        di.id = AuxSize ;
        AuxSize += size ;
        DataMap[name] = di ;
	return di.id ;
      } else {
        if(size != ii->second.size) {
          std::cerr << "incompatible sizes in duplicate requests for '"
               << name << "' in AuxiliaryInfoDatabase::insertItem()"
               << std::endl ;
        }
	return ii->second.id ;
      }
    }
    int getAuxSize() const { return AuxSize ; }
    
    int getItemLoc(string name) const {
      std::map<string,data_info>::const_iterator ii ;
      ii = DataMap.find(name) ;
      if(ii == DataMap.end()) {
        return -1 ;
      }
      return ii->second.id ;
    }

    int getItemSize(string name) const {
      std::map<string,data_info>::const_iterator ii ;
      ii = DataMap.find(name) ;
      if(ii == DataMap.end()) {
        return -1 ;
      }
      return ii->second.size ;
    }
  } ;
      
      
      
}



namespace Loci {
  template<> struct data_schema_traits<lagrangianP::particleModel> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;

    typedef char Converter_Base_Type ;
    typedef StringStreamConverter<lagrangianP::particleModel> Converter_Type ;
  } ;
}


#endif
