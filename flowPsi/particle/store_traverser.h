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

#ifndef STORE_TRAVERSER_H
#define STORE_TRAVERSER_H

#include <Loci.h>

namespace lagrangianP {

  // this class is used to traverse the
  // Loci store container class hierarchy and perform
  // user defined actions
  class store_traverser ;
  typedef Loci::CPTR<store_traverser> store_traverserP ;
  
  class store_traverser: public Loci::CPTR_type {
  public:
    virtual store_traverserP clone() = 0 ;
    virtual Loci::storeRepP get_storeRep() = 0 ;
    virtual void set_args(Loci::storeRepP& srp) = 0 ;
    virtual void traverse(const sequence& seq) = 0 ;
  } ;

  // these are designed to perform transformations of
  // the store contents based on user defined unary functions,
  // it relies on using the store_traverser to traverse
  // the store contents
  template <class T, class UnaryFunction>
  class store_transformer: public store_traverser {
    T s ;
    UnaryFunction transform ;
  public:
    virtual store_traverserP
    clone() {
      return store_traverserP(new store_transformer<T,UnaryFunction>) ;
    }

    virtual Loci::storeRepP
    get_storeRep() {
      return s.Rep() ;
    }

    virtual void
    set_args(Loci::storeRepP& srp) {
      s.setRep(srp) ;
    }

    virtual void
    traverse(const sequence& seq) {
      for(sequence::const_iterator si=seq.begin();
          si!=seq.end();++si)
        transform(s[*si]) ;
    }
    
  } ;

  // provide specialized store_transformers for various
  // Loci stores types
  template <class T, class UnaryFunction>
  class store_transformer<param<T>,
                          UnaryFunction>: public store_traverser {
    param<T> p ;
    UnaryFunction transform ;
  public:
    virtual store_traverserP
    clone() {
      return
        store_traverserP(new store_transformer<param<T>,
                                               UnaryFunction>) ;
    }

    virtual Loci::storeRepP
    get_storeRep() {
      return p.Rep() ;
    }

    virtual void
    set_args(Loci::storeRepP& srp) {
      p.setRep(srp) ;
    }

    virtual void
    traverse(const sequence& seq) {
      transform(*p) ;
    }    
  } ;

  template<class T, class UnaryFunction>
  class store_transformer<storeVec<T>,
                          UnaryFunction>: public store_traverser {
    storeVec<T> s ;
    UnaryFunction transform ;
  public:    
    virtual store_traverserP
    clone() {
      return
        store_traverserP(new store_transformer<storeVec<T>,
                                               UnaryFunction>) ;
    }

    virtual Loci::storeRepP
    get_storeRep() {
      return s.Rep() ;
    }

    virtual void
    set_args(Loci::storeRepP& srp) {
      s.setRep(srp) ;
    }

    virtual void
    traverse(const sequence& seq) {
      int size = s.vecSize() ;
      for(sequence::const_iterator si=seq.begin();
          si!=seq.end();++si) {
        Loci::Vect<T> v = s[*si] ;
        for(int i=0;i<size;++i)
          transform(v[i]) ;
      }
    }
  } ;

  template<class T, class UnaryFunction>
  class store_transformer<storeMat<T>,
                          UnaryFunction>: public store_traverser {
    storeMat<T> s ;
    UnaryFunction transform ;
  public:
    virtual store_traverserP
    clone() {
      return
        store_traverserP(new store_transformer<storeMat<T>,
                                               UnaryFunction>) ;
    }

    virtual Loci::storeRepP
    get_storeRep() {
      return s.Rep() ;
    }

    virtual void
    set_args(Loci::storeRepP& srp) {
      s.setRep(srp) ;
    }

    virtual void
    traverse(const sequence& seq) {
      int size = s.vecSize() ;
      size*=size ;
      for(sequence::const_iterator si=seq.begin();
          si!=seq.end();++si) {
        Loci::Mat<T> m = s[*si] ;
        for(int i=0;i<size;++i)
          transform(m[i]) ;
      }
    }
  } ;

  template<class T, class UnaryFunction>
  class store_transformer<multiStore<T>,
                          UnaryFunction>: public store_traverser {
    multiStore<T> s ;
    UnaryFunction transform ;
  public:
    virtual store_traverserP
    clone() {
      return
        store_traverserP(new store_transformer<multiStore<T>,
                                               UnaryFunction>) ;
    }

    virtual Loci::storeRepP
    get_storeRep() {
      return s.Rep() ;
    }

    virtual void
    set_args(Loci::storeRepP& srp) {
      s.setRep(srp) ;
    }

    virtual void
    traverse(const sequence& seq) {
      for(sequence::const_iterator si=seq.begin();
          si!=seq.end();++si) {
        Loci::Vect<T> v = s[*si] ;
        int size = v.getSize() ;
        for(int i=0;i<size;++i)
          transform(v[i]) ;
      }
    }
  } ;

  // a function to create store_transformer objects
  // return a pointer to the base class
  template<class T, class UnaryFunction>
  store_traverserP create_store_transformer() {
    return store_traverserP(new store_transformer<T,UnaryFunction>) ;
  }

  // we predefine one UnaryFunction for transformation.
  // this one will set default value (T()) for the passed in value
  template <typename T>
  struct set_default_value {
    void
    operator()(T& v) {
      v = T() ;
    }
  } ;

  // specialized for Loci::vector3d<T> and Loci::vector2d<T>
  // class since their default constructors do not have proper
  // initialization implemented (which they really should)
  template <typename T>
  struct set_default_value<vector2d<T> > {
    void operator()(vector2d<T>& v) {
      v = vector2d<T>(T(), T()) ;
    }
  } ;
  template <typename T>
  struct set_default_value<vector3d<T> > {
    void operator()(vector3d<T>& v) {
      v = vector3d<T>(T(), T(), T()) ;
    }
  } ;

  // here is a function for creating Loci joiners
  typedef Loci::CPTR<Loci::joiner> loci_joinerP ;

  template<typename T, class Op>
  loci_joinerP create_loci_joiner() {
    return loci_joinerP(new Loci::joinOp<T,Op>) ;
  }
  
} // end of namespace lagrangianP

#endif
