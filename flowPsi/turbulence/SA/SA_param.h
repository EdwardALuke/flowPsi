#ifndef SA_PARAM_H
#define SA_PARAM_H

namespace flowPsi {

  struct Spa_All_param {
    real cb1,k_coeff,cv1,cb2,cw2,cw3,sigma ; //coefficients in Spalart_Allmaras
    //                                         //model

    Spa_All_param()   {
      cb1=0.1355; //production coeffient
      k_coeff=0.41; // coeffient for constructing vorticity \tilde{S}
      cv1=7.1;  //coeffient for constructing f_v1 which is the ratio of
      //        //eddy viscoisty to molecular viscosity
      cb2=0.622; // coeffient in diffusion term
      cw2=0.3; // coeffient in destruction term
      cw3=2.0; //coeffienct in destruction term
      sigma=2./3.; //coeffienct in diffusion term
    }
  } ;
}

namespace Loci {
  template<> struct data_schema_traits<flowPsi::Spa_All_param> {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      CompoundDatatypeP ct = CompoundFactory(flowPsi::Spa_All_param()) ;
      LOCI_INSERT_TYPE(ct,flowPsi::Spa_All_param,cb1) ;
      LOCI_INSERT_TYPE(ct,flowPsi::Spa_All_param,k_coeff) ;
      LOCI_INSERT_TYPE(ct,flowPsi::Spa_All_param,cv1) ;
      LOCI_INSERT_TYPE(ct,flowPsi::Spa_All_param,cb2) ;
      LOCI_INSERT_TYPE(ct,flowPsi::Spa_All_param,cw2) ;
      LOCI_INSERT_TYPE(ct,flowPsi::Spa_All_param,cw3) ;
      LOCI_INSERT_TYPE(ct,flowPsi::Spa_All_param,sigma) ;
      return DatatypeP(ct) ;
    }
  } ;
}

#endif
