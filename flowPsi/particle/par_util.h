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

// this file contains parallel utility functions
#ifndef PARTICLE_PAR_UTIL_H
#define PARTICLE_PAR_UTIL_H

#include <Loci.h>
#include <iostream>
#include <vector>
#include <list>
#include <iterator>
#include <utility>
#include <algorithm>
#include <mpi.h>

#ifdef NO_CSTDLIB
#include <stdlib.h>
#include <math.h>
#else
#include <cstdlib>
#include <cmath>
#endif

#include "util.h"

namespace Loci {
  extern void ORBPartition(const std::vector<vector3d<float> >&,
                           std::vector<int>&, MPI_Comm) ;
}

namespace lagrangianP {

  // this function transposes the passed in vector<entitySet>
  // by an all to all personalized communication
  std::vector<entitySet>
  transpose_vector_entitySet(const std::vector<entitySet>& in) ;
  // this is an optimized version that communicates the intervals
  // in an entitySet if that is less data to communicate, if the
  // entitySet is sparse, then it will choose to send/recv the elements
  // directly. this one also avoids copy on return, instead using
  // a reference passed in
  void
  transpose_vector_entitySet_opt(const std::vector<entitySet>& in,
                                 std::vector<entitySet>& out) ;

  // transpose a vector<sequence>
  std::vector<sequence>
  transpose_vector_sequence(const std::vector<sequence>& in) ;
  // an optimized version similar to the entitySet version
  void
  transpose_vector_sequence_opt(const std::vector<sequence>& in,
                                std::vector<entitySet>& out) ;
  
  // a utility that returns the global sum
  inline int
  global_sum(int l) {
    int g ;
    MPI_Allreduce(&l, &g, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD) ;
    return g ;
  }

  // a utility function that takes an entitySet from a processor
  // and returns a vector of entitySet gathered from all processors
  std::vector<entitySet>
  gather_all_entitySet(const entitySet& eset) ;
  // an optimized version
  void
  gather_all_entitySet_opt(const entitySet& eset,
                           std::vector<entitySet>& out) ;

  // this function takes an entitySet (in the fact_db's initial
  // global numbering) and a Map, it returns the image of the Map
  // based on the passed in entitySet as the domain. The returned
  // entitySet is also in fact_db's initial global numbering.
  // Note: the passed in entitySet may not reside on a single
  // process, hence this is a distributed Map image function.
  // "facts" is the pointer to the current fact_db in use.
  entitySet
  distributed_map_image(const entitySet& domain,
                        Loci::MapRepP m, Loci::fact_db* facts) ;

#ifdef TO_BE_REMOVED
  // the functions in this section are obsoleted, we should remove
  // them in the future.

  
  // template functions need to stay in the header file
  
  // this function expands a store<T> for the passed in domain (in
  // global numbering), it returns a new store that includes the
  // expanded domain (with the original domain).
  // NOTE: everything in this function is in global numbering
  template<typename T> Loci::storeRepP
  expand_store(const store<T>& s,
               const entitySet& edom, Loci::fact_db& facts) {
    // first get fact_db's init partition in global numbering
    const std::vector<entitySet>& ptn = facts.get_init_ptn() ;
    // figure out domain distributions
    std::vector<entitySet> edom_dist(Loci::MPI_processes) ;
    for(int i=0;i<Loci::MPI_processes;++i) {
      edom_dist[i] = edom & ptn[i] ;
    }
    // transpose "edom_dist"
    std::vector<entitySet> edom_dist_t =
      transpose_vector_entitySet(edom_dist) ;
    // get the storeRepP first for the passed in store<T>
    Loci::storeRepP srp = s.Rep() ;

    // we will need to pack the contents for these domains and
    // communicate them back to the requesting processes.
    // we will first need to prepare the sending buffer size
    int* send_counts = new int[Loci::MPI_processes] ;
    for(int i=0;i<Loci::MPI_processes;++i)
      send_counts[i] = srp->pack_size(edom_dist_t[i]) ;

    int* send_displs = new int[Loci::MPI_processes] ;
    send_displs[0] = 0 ;
    for(int i=1;i<Loci::MPI_processes;++i)
      send_displs[i] = send_displs[i-1] + send_counts[i-1] ;

    // then communicate the size first
    int* recv_counts = new int[Loci::MPI_processes] ;
    MPI_Alltoall(send_counts, 1, MPI_INT,
                 recv_counts, 1, MPI_INT, MPI_COMM_WORLD) ;

    int* recv_displs = new int[Loci::MPI_processes] ;
    recv_displs[0] = 0 ;
    for(int i=1;i<Loci::MPI_processes;++i)
      recv_displs[i] = recv_displs[i-1] + recv_counts[i-1] ;

    // we have computed the size, now we need to pack the buffer
    int buf_size = 0 ;
    for(int i=0;i<Loci::MPI_processes;++i)
      buf_size += send_counts[i] ;

    unsigned char* send_buf = new unsigned char[buf_size] ;
    unsigned char** send_ptr = new unsigned char*[Loci::MPI_processes] ;
    send_ptr[0] = send_buf ;
    for(int i=1;i<Loci::MPI_processes;++i)
      send_ptr[i] = send_ptr[i-1] + send_counts[i-1] ;

    // pack the buffer
    for(int i=0;i<Loci::MPI_processes;++i) {
      int position = 0 ;
      srp->pack(send_ptr[i], position, send_counts[i], edom_dist_t[i]) ;
    }

    delete[] send_ptr ;

    // compute the receive sequence first
    std::vector<sequence> recv_seq(Loci::MPI_processes) ;
    for(int i=0;i<Loci::MPI_processes;++i) {
      recv_seq[i] = sequence(edom_dist_t[i]) ;
    }
    recv_seq = transpose_vector_sequence(recv_seq) ;

    // then allocate receive buffer
    int recv_size = 0 ;
    for(int i=0;i<Loci::MPI_processes;++i)
      recv_size += recv_counts[i] ;

    unsigned char* recv_buf = new unsigned char[recv_size] ;

    // communicate
    MPI_Alltoallv(send_buf, send_counts,
                  send_displs, MPI_PACKED,
                  recv_buf, recv_counts,
                  recv_displs, MPI_PACKED, MPI_COMM_WORLD) ;

    delete[] send_counts ;
    delete[] send_displs ;
    delete[] send_buf ;

    // unpack recv buffer
    // we will need to create a new store first
    store<T> ns ;
    // the domain is the union of the origiinal and the "edom" as passed in.
    entitySet ns_dom = s.domain() + edom ;
    ns.allocate(ns_dom) ;
    Loci::storeRepP ns_srp = ns.Rep() ;

    unsigned char** recv_ptr = new unsigned char*[Loci::MPI_processes] ;
    recv_ptr[0] = recv_buf ;
    for(int i=1;i<Loci::MPI_processes;++i)
      recv_ptr[i] = recv_ptr[i-1] + recv_counts[i-1] ;

    for(int i=0;i<Loci::MPI_processes;++i) {
      // unpack
      int position = 0 ;
      ns_srp->unpack(recv_ptr[i], position, recv_counts[i], recv_seq[i]) ;
    }

    delete[] recv_buf ;
    delete[] recv_ptr ;
    delete[] recv_counts ;
    delete[] recv_displs ;

    // fill in the original stuff
    entitySet original_dom = s.domain() ;
    for(entitySet::const_iterator ei=original_dom.begin();
        ei!=original_dom.end();++ei)
      ns[*ei] = s[*ei] ;

    return ns_srp ;
  }

  // this function expands a store<T> for the passed in domain (in
  // global numbering), however it returns a new store that has the
  // expanded domain only, and the expanded domain is also renumbered
  // based on the passed in renumber map. NOTE: the local original
  // store "s" is in local numbering.
  template<typename T> Loci::storeRepP
  expand_store(const store<T>& s, const entitySet& edom,
               const dMap& remap, Loci::fact_db* facts) {
    // if the global edom is empty, then we would return
    if(Loci::GLOBAL_AND(edom==EMPTY)) {
      store<T> ns ;
      return ns.Rep() ;
    }
    // first get fact_db's init partition in global numbering
    const std::vector<entitySet>& ptn = facts->get_init_ptn() ;
    // then get the distribute_info structure
    Loci::fact_db::distribute_infoP df = facts->get_distribute_info() ;
    // figure out domain distributions
    std::vector<entitySet> edom_dist(Loci::MPI_processes) ;
    for(int i=0;i<Loci::MPI_processes;++i) {
      edom_dist[i] = edom & ptn[i] ;
    }

    // transpose "edom_dist"
    std::vector<entitySet> edom_dist_t =
      transpose_vector_entitySet(edom_dist) ;

    // convert to local numbering
    for(int i=0;i<Loci::MPI_processes;++i)
      edom_dist_t[i] = remap_entitySet(edom_dist_t[i], df->g2l) ;

    // get the storeRepP first for the passed in store<T>
    Loci::storeRepP srp = s.Rep() ;

    // we will need to pack the contents for these domains and
    // communicate them back to the requesting processes.
    // we will first need to prepare the sending buffer size
    int* send_counts = new int[Loci::MPI_processes] ;
    for(int i=0;i<Loci::MPI_processes;++i)
      send_counts[i] = srp->pack_size(edom_dist_t[i]) ;

    int* send_displs = new int[Loci::MPI_processes] ;
    send_displs[0] = 0 ;
    for(int i=1;i<Loci::MPI_processes;++i)
      send_displs[i] = send_displs[i-1] + send_counts[i-1] ;

    // then communicate the size first
    int* recv_counts = new int[Loci::MPI_processes] ;
    MPI_Alltoall(send_counts, 1, MPI_INT,
                 recv_counts, 1, MPI_INT, MPI_COMM_WORLD) ;

    int* recv_displs = new int[Loci::MPI_processes] ;
    recv_displs[0] = 0 ;
    for(int i=1;i<Loci::MPI_processes;++i)
      recv_displs[i] = recv_displs[i-1] + recv_counts[i-1] ;

    // we have computed the size, now we need to pack the buffer
    int buf_size = 0 ;
    for(int i=0;i<Loci::MPI_processes;++i)
      buf_size += send_counts[i] ;

    unsigned char* send_buf = new unsigned char[buf_size] ;
    unsigned char** send_ptr = new unsigned char*[Loci::MPI_processes] ;
    send_ptr[0] = send_buf ;
    for(int i=1;i<Loci::MPI_processes;++i)
      send_ptr[i] = send_ptr[i-1] + send_counts[i-1] ;

    // pack the buffer
    for(int i=0;i<Loci::MPI_processes;++i) {
      int position = 0 ;
      srp->pack(send_ptr[i], position, send_counts[i], edom_dist_t[i]) ;
    }

    delete[] send_ptr ;

    // compute the receive sequence first
    std::vector<sequence> recv_seq(Loci::MPI_processes) ;
    for(int i=0;i<Loci::MPI_processes;++i) {
      recv_seq[i] = sequence(edom_dist_t[i]) ;
      recv_seq[i] = remap_sequence(recv_seq[i], df->l2g) ;
    }
    recv_seq = transpose_vector_sequence(recv_seq) ;

    // then allocate receive buffer
    int recv_size = 0 ;
    for(int i=0;i<Loci::MPI_processes;++i)
      recv_size += recv_counts[i] ;

    unsigned char* recv_buf = new unsigned char[recv_size] ;

    // communicate
    MPI_Alltoallv(send_buf, send_counts,
                  send_displs, MPI_PACKED,
                  recv_buf, recv_counts,
                  recv_displs, MPI_PACKED, MPI_COMM_WORLD) ;

    delete[] send_counts ;
    delete[] send_displs ;
    delete[] send_buf ;

    // unpack recv buffer
    // we will need to create a new store first
    store<T> ns ;
    // the domain is "edom" as passed in, but we need to renumber it
    entitySet ns_dom = remap_entitySet(edom, remap) ;
    ns.allocate(ns_dom) ;
    Loci::storeRepP ns_srp = ns.Rep() ;

    unsigned char** recv_ptr = new unsigned char*[Loci::MPI_processes] ;
    recv_ptr[0] = recv_buf ;
    for(int i=1;i<Loci::MPI_processes;++i)
      recv_ptr[i] = recv_ptr[i-1] + recv_counts[i-1] ;

    for(int i=0;i<Loci::MPI_processes;++i) {
      // compute the domain in the new local numbering
      sequence seq = remap_sequence(recv_seq[i], remap) ;
      int position = 0 ;
      ns_srp->unpack(recv_ptr[i], position, recv_counts[i], seq) ;
    }

    delete[] recv_buf ;
    delete[] recv_ptr ;
    delete[] recv_counts ;
    delete[] recv_displs ;

    return ns_srp ;
  }

  // this function expands a storeVec<T> for the passed in domain (in
  // global numbering), however it returns a new storeVec that has the
  // expanded domain only, and the expanded domain is also renumbered
  // based on the passed in renumber map. NOTE: the local original
  // storeVec "s" is in local numbering.
  template<typename T> Loci::storeRepP
  expand_storeVec(const storeVec<T>& s, const entitySet& edom,
                  const dMap& remap, Loci::fact_db* facts) {
    // if the global edom is empty, then we would return
    if(Loci::GLOBAL_AND(edom==EMPTY)) {
      storeVec<T> ns ;
      return ns.Rep() ;
    }
    // first get fact_db's init partition in global numbering
    const std::vector<entitySet>& ptn = facts->get_init_ptn() ;
    // then get the distribute_info structure
    Loci::fact_db::distribute_infoP df = facts->get_distribute_info() ;
    // figure out domain distributions
    std::vector<entitySet> edom_dist(Loci::MPI_processes) ;
    for(int i=0;i<Loci::MPI_processes;++i) {
      edom_dist[i] = edom & ptn[i] ;
    }
    // transpose "edom_dist"
    std::vector<entitySet> edom_dist_t =
      transpose_vector_entitySet(edom_dist) ;
    // convert to local numbering
    for(int i=0;i<Loci::MPI_processes;++i)
      edom_dist_t[i] = remap_entitySet(edom_dist_t[i], df->g2l) ;

    // get the storeRepP first for the passed in storeVec<T>
    Loci::storeRepP srp = s.Rep() ;

    // we will need to pack the contents for these domains and
    // communicate them back to the requesting processes.
    // we will first need to prepare the sending buffer size
    int* send_counts = new int[Loci::MPI_processes] ;
    for(int i=0;i<Loci::MPI_processes;++i)
      send_counts[i] = srp->pack_size(edom_dist_t[i]) ;

    int* send_displs = new int[Loci::MPI_processes] ;
    send_displs[0] = 0 ;
    for(int i=1;i<Loci::MPI_processes;++i)
      send_displs[i] = send_displs[i-1] + send_counts[i-1] ;

    // then communicate the size first
    int* recv_counts = new int[Loci::MPI_processes] ;
    MPI_Alltoall(send_counts, 1, MPI_INT,
                 recv_counts, 1, MPI_INT, MPI_COMM_WORLD) ;

    int* recv_displs = new int[Loci::MPI_processes] ;
    recv_displs[0] = 0 ;
    for(int i=1;i<Loci::MPI_processes;++i)
      recv_displs[i] = recv_displs[i-1] + recv_counts[i-1] ;

    // we have computed the size, now we need to pack the buffer
    int buf_size = 0 ;
    for(int i=0;i<Loci::MPI_processes;++i)
      buf_size += send_counts[i] ;

    unsigned char* send_buf = new unsigned char[buf_size] ;
    unsigned char** send_ptr = new unsigned char*[Loci::MPI_processes] ;
    send_ptr[0] = send_buf ;
    for(int i=1;i<Loci::MPI_processes;++i)
      send_ptr[i] = send_ptr[i-1] + send_counts[i-1] ;

    // pack the buffer
    for(int i=0;i<Loci::MPI_processes;++i) {
      int position = 0 ;
      srp->pack(send_ptr[i], position, send_counts[i], edom_dist_t[i]) ;
    }

    delete[] send_ptr ;

    // compute the receive sequence first
    std::vector<sequence> recv_seq(Loci::MPI_processes) ;
    for(int i=0;i<Loci::MPI_processes;++i) {
      recv_seq[i] = sequence(edom_dist_t[i]) ;
      recv_seq[i] = remap_sequence(recv_seq[i], df->l2g) ;
    }
    recv_seq = transpose_vector_sequence(recv_seq) ;

    // then allocate receive buffer
    int recv_size = 0 ;
    for(int i=0;i<Loci::MPI_processes;++i)
      recv_size += recv_counts[i] ;

    unsigned char* recv_buf = new unsigned char[recv_size] ;

    // communicate
    MPI_Alltoallv(send_buf, send_counts,
                  send_displs, MPI_PACKED,
                  recv_buf, recv_counts,
                  recv_displs, MPI_PACKED, MPI_COMM_WORLD) ;

    delete[] send_counts ;
    delete[] send_displs ;
    delete[] send_buf ;

    // unpack recv buffer
    // we will need to create a new store first
    storeVec<T> ns ;
    ns.setVecSize(s.vecSize()) ;
    // the domain is "edom" as passed in, but we need to renumber it
    entitySet ns_dom = remap_entitySet(edom, remap) ;
    ns.allocate(ns_dom) ;
    Loci::storeRepP ns_srp = ns.Rep() ;

    unsigned char** recv_ptr = new unsigned char*[Loci::MPI_processes] ;
    recv_ptr[0] = recv_buf ;
    for(int i=1;i<Loci::MPI_processes;++i)
      recv_ptr[i] = recv_ptr[i-1] + recv_counts[i-1] ;

    for(int i=0;i<Loci::MPI_processes;++i) {
      // compute the domain in the new local numbering
      sequence seq = remap_sequence(recv_seq[i], remap) ;
      int position = 0 ;
      ns_srp->unpack(recv_ptr[i], position, recv_counts[i], seq) ;
    }

    delete[] recv_buf ;
    delete[] recv_ptr ;
    delete[] recv_counts ;
    delete[] recv_displs ;

    return ns_srp ;
  }

  // this function expands a Map for the passed in domain (in
  // global numbering), however it returns a new Map that has the
  // expanded domain only, and the expanded domain is also renumbered
  // based on the passed in renumber map
  Loci::storeRepP
  expand_Map(const Map& m, const entitySet& edom,
             const dMap& remap, Loci::fact_db* facts) ;


  // This is a special expand_Map function, everything is
  // the same as the expand_Map function above, except for
  // that in the newly expanded Map, only its domain is
  // renumbered, the image of the Map is left unchanged (not remapped).
  Loci::storeRepP
  expand_Map2(const Map& m, const entitySet& edom,
              const dMap& remap, Loci::fact_db* facts) ;

  // this function expands a multiMap for the passed in domain (in
  // global numbering), however it returns a new multiMap that has the
  // expanded domain only, and the expanded domain is also renumbered
  // based on the passed in renumber map
  Loci::storeRepP
  expand_multiMap(const multiMap& m, const entitySet& edom,
                  const dMap& remap, Loci::fact_db* facts) ;
#endif

  ////////////////////////////////////////////////////////////
  // these are functions that perform particle redistribution

  // this function generates a process id map based on a
  // supplied process entitySet partition.
  void
  generate_process_map(const std::vector<entitySet>& ptn,
                       const std::vector<int>& items,
                       std::vector<int>& pmap) ;
  
  // this function will distribute an std::vector<T> according
  // to a supplied index array that maps the ith element in
  // the vector to a process.
  // NOTE: it is assumed that the size of an object of T == sizeof(T)
  template<typename T> void
  distribute_vector(const std::vector<int>& pmap,
                    std::vector<T>& v, MPI_Comm comm) {
    if(pmap.size() != v.size()) {
      throw Loci::StringError
        ("distribute_vector failed! pmap.size != v.size") ;
    }
    int np = 0 ;
    MPI_Comm_size(comm, &np) ;
    // first communicate the send/recv size
    std::vector<int> send_items(np, 0) ;
    for(size_t i=0;i!=pmap.size();++i)
      send_items[pmap[i]]++ ;

    std::vector<int> recv_items(np, 0) ;
    MPI_Alltoall(&send_items[0], 1, MPI_INT,
                 &recv_items[0], 1, MPI_INT, comm) ;

    // now obtain the true send/recv size (in bytes)
    std::vector<int> send_sizes(np, 0) ;
    for(int i=0;i!=np;++i)
      send_sizes[i] = send_items[i] * sizeof(T) ;

    std::vector<int> recv_sizes(np, 0) ;
    for(int i=0;i!=np;++i)
      recv_sizes[i] = recv_items[i] * sizeof(T) ;

    std::vector<int> send_displs(np, 0) ;
    std::vector<int> recv_displs(np, 0) ;
    for(int i=1;i<np;++i) {
      send_displs[i] = send_displs[i-1] + send_items[i-1] ;
      recv_displs[i] = recv_displs[i-1] + recv_items[i-1] ;
    }

    int tot_send_items = 0 ;
    int tot_recv_items = 0 ;
    for(int i=0;i<np;++i) {
      tot_send_items += send_items[i] ;
      tot_recv_items += recv_items[i] ;
    }
    // allocate the send buffer
    std::vector<T> send_buffer(tot_send_items) ;
    // pack the send buffer
    std::vector<int> buffer_idx = send_displs ;
    for(size_t i=0;i!=pmap.size();++i) {
      int pid = pmap[i] ;
      int& po = buffer_idx[pid] ;
      send_buffer[po] = v[i] ;
      po++ ;
    }

    // rescale displs in bytes
    for(int i=0;i!=np;++i) {
      send_displs[i] = send_displs[i] * sizeof(T) ;
      recv_displs[i] = recv_displs[i] * sizeof(T) ;
    }
    
    // allocate recv buffer
    std::vector<T> recv_buffer(tot_recv_items) ;
    MPI_Alltoallv(&send_buffer[0], &send_sizes[0],
                  &send_displs[0], MPI_BYTE,
                  &recv_buffer[0], &recv_sizes[0],
                  &recv_displs[0], MPI_BYTE, comm) ;

    v.swap(recv_buffer) ;
  }

  // this is function that balances the items in the passed
  // in sequence among the participating processes in the
  // MPI communicator. the main difference with the current
  // implementation in the Loci code base (the function
  // "balanceDistribution" in "parSampleSort.h") is that
  // this one does not assume the objects in the sequence
  // are a primitive type (i.e., they could have pointers or
  // other complex data-structures). this function therefore
  // rely on the objects in the sequence have the "pack"
  // and "unpack" methods defined. besides, this version
  // also generalized the interface so that it is more
  // flexible. after this function, each process in the
  // MPI communicator will have approximately the same
  // number of objects while also maintaining the original
  // global ordering.
  template <class ForwardIterator, class OutputIterator> OutputIterator
  balance_sequence(ForwardIterator first, ForwardIterator last,
                   OutputIterator result, MPI_Comm comm) {
    // get total number of processes
    int p = 0 ;
    MPI_Comm_size(comm,&p) ;

    if(p == 1) {
      // if serial run, just copy the contents
      for(;first!=last;++first,++result)
        *result = *first ;
      return result ;
    }

    // get process rank
    int r=0 ;
    MPI_Comm_rank(comm,&r) ;

    // compute the sequence size
    int sz = std::distance(first, last) ;
    // buffer to store the size of sequence on all processes
    std::vector<int> sizes(p) ;
    MPI_Allgather(&sz,1,MPI_INT,&sizes[0],1,MPI_INT,comm) ;
    long long tot_size = 0  ;
    for(int i=0;i<p;++i) {
      tot_size += sizes[i] ;
    }
    // bsizes holds the balanced size on each process
    std::vector<long long> bsizes(p) ;
    int psz = tot_size/p ;
    int rem = tot_size%p ;
    for(int i=0;i<p;++i)
      bsizes[i] = psz + (i<rem?1:0) ;

    // recv_size is the buffer to hold size of
    // message from other processes
    std::vector<int> recv_size(p,0) ;
    // the send size
    std::vector<int> send_size(p,0) ;
    // sdist is the buffer of the current object distribution
    // in the sequence on each process. rdist is the buffer for
    // the balanced object distribution (global index range)
    // in the global sequence on each process.
    std::vector<long long> sdist(p+1),rdist(p+1) ;
    sdist[0] = 0 ;
    rdist[0] = 0 ;
    for(int i=0;i<p;++i) {
      sdist[i+1] = sizes[i]+sdist[i] ;
      rdist[i+1] = bsizes[i]+rdist[i] ;
    }

    // we will first need to communicate the size of the message
    
    // Perform Irecvs first
    std::vector<MPI_Request> requests(p) ;
    int req = 0 ;

    // [i1,i2] defines the global index interval of the
    // balanced sequence on process "r"
    long long i1 = rdist[r] ;
    long long i2 = rdist[r+1]-1 ;

    // to compute the receive range is to intersect [i1,i2]
    // with the current global index distribution on all
    // processes
    for(int i=0;i<p;++i) {
      if(sdist[i]<=i2 && sdist[i+1]-1>=i1) {
        // this means that process "r" has an intersection
        // with process "i" and is therefore required to
        // obtain receive size from process i.
        if(i == r) { // local copy
          // [li,ri] defines the intersected interval with
          // process "i"
          long long li = std::max(i1,sdist[i]) ;
          long long ri = std::min(i2,sdist[i+1]-1) ;
          int len = ri-li+1 ;

          // s1 is the offset to the send buffer,
          // in other words, iterator "first+s1"
          // points to the first object to be sent
          int s1 = li-sdist[r] ;

          ForwardIterator aux = first ;
          std::advance(aux, s1) ;
          
          for(int j=0;j<len;++j,++aux)
            recv_size[i] += aux->pack_size() ;
          
        } else {
          MPI_Irecv(&recv_size[i],1,
                    MPI_INT,i,55,comm,&requests[req++]) ;
        }
      }
    }

    // Perform sends of message size

    // [i1,i2] defines the global index interval of the
    // current object distribution on process "r"
    i1 = sdist[r] ;
    i2 = sdist[r+1]-1 ;
    // then this interval is intersected with the
    // balanced object global index interval on
    // each process to figure out the sub-interval
    // to be sent to each process
    for(int i=0;i<p;++i) {
      if(i != r && rdist[i]<=i2 && rdist[i+1]-1>=i1) {
        // intersected, [li,ri] defines the
        // intersected interval with process "i"'s
        // balanced index interval
        long long li = std::max(i1,rdist[i]) ;
        long long ri = std::min(i2,rdist[i+1]-1) ;
        int len = ri-li+1 ;
        // s1 is the offset to the send buffer,
        // i.e., iterator "first+s1" points to
        // the first object to be sent
        int s1 = li-sdist[r] ;
        // compute the sent size
        ForwardIterator aux = first ;
        std::advance(aux, s1) ;

        for(int j=0;j<len;++j,++aux)
          send_size[i] += aux->pack_size() ;

        MPI_Send(&send_size[i],1,MPI_INT,i,55,comm) ;
      }
    }

    if(req > 0) {
      std::vector<MPI_Status> status(p) ;
      MPI_Waitall(req,&requests[0],&status[0]) ;
    }
    
    // we are now ready to pack the message and send/recv
    // first we'll allocate the recv buffer
    int total_recv_size = 0 ;
    for(int i=0;i<p;++i)
      total_recv_size += recv_size[i] ;
    std::vector<int> recv_displs(p,0) ;
    for(int i=1;i<p;++i)
      recv_displs[i] = recv_displs[i-1] + recv_size[i-1] ;

    unsigned char* recv_buffer = new unsigned char[total_recv_size] ;

    req = 0 ;
    // perform Irecvs first
    i1 = rdist[r] ;
    i2 = rdist[r+1]-1 ;
    for(int i=0;i<p;++i) {
      if(sdist[i]<=i2 && sdist[i+1]-1>=i1) {
        int li = std::max(i1,sdist[i]) ;
        int ri = std::min(i2,sdist[i+1]-1) ;
        int len = ri-li+1 ;

        unsigned char* recv_start = recv_buffer + recv_displs[i] ;
        if(i == r) { // local copy
          int s1 = li-sdist[r] ;

          ForwardIterator aux = first ;
          std::advance(aux, s1) ;

          size_t position = 0 ;
          for(int j=0;j<len;++j,++aux)
            aux->pack(recv_start, position) ;

        } else {
          MPI_Irecv(recv_start,recv_size[i],
                    MPI_BYTE,i,66,comm,&requests[req++]) ;
        }
      }
    }

    // pack and send
    i1 = sdist[r] ;
    i2 = sdist[r+1]-1 ;
    for(int i=0;i<p;++i) {
      if(i != r && rdist[i]<=i2 && rdist[i+1]-1>=i1) {
        long long li = std::max(i1,rdist[i]) ;
        long long ri = std::min(i2,rdist[i+1]-1) ;
        int len = ri-li+1 ;
        int s1 = li-sdist[r] ;

        ForwardIterator aux = first ;
        std::advance(aux, s1) ;

        unsigned char* send_buffer = new unsigned char[send_size[i]] ;

        size_t position = 0 ;
        for(int j=0;j<len;++j,++aux)
          aux->pack(send_buffer, position) ;
        
        MPI_Send(send_buffer,send_size[i],MPI_BYTE,i,66,comm) ;

        delete[] send_buffer ;
      }
    }

    if(req > 0) {
      std::vector<MPI_Status> status(p) ;
      MPI_Waitall(req,&requests[0],&status[0]) ;
    }

    // finished, now we need to extract the recv buffer

    size_t position = 0 ;
    typename std::iterator_traits<ForwardIterator>::value_type tmp ;
    for(int i=0;i<bsizes[r];++i,++result) {
      tmp.unpack(recv_buffer, position) ;
      *result = tmp ;
    }

    delete[] recv_buffer ;

    return result ;
  }

  // this is a simplified version that just returns the
  // length of the balanced sequence, but doesn't actually
  // do the balance
  template <class InputIterator, class OutputIterator> int
  balance_sequence(InputIterator first,
                   InputIterator last, MPI_Comm comm) {
    // get total number of processes
    int p = 0 ;
    MPI_Comm_size(comm,&p) ;
    
    // get process rank
    int r=0 ;
    MPI_Comm_rank(comm,&r) ;
    
    // compute the sequence size
    long long int sz = std::distance(first, last) ;
    
    // buffer to store the size of sequence on all processes
    long long int tot_size = 0 ;
    MPI_Allreduce(&sz, &tot_size, 1, MPI_LONG_LONG_INT, MPI_SUM, comm) ;
    
    int psz = tot_size/p ;
    int rem = tot_size%p ;
    
    int bsize = psz + ( (r<rem) ? 1 : 0) ;

    return bsize ;
  }

  // a helper function used for the "orb_partition_sequence"
  // it implements the main communication routines, but
  // assumes that the sequence passed in is indeed valid
  // for performing the orb partition. the addtional parameter
  // sequence_len indicates the sequence length (to avoid
  // recompute it again.
  template <class ForwardIterator, class OutputIterator> OutputIterator
  orb_partition_sequence_helper
  (ForwardIterator first, ForwardIterator last,
   size_t len, OutputIterator result, MPI_Comm comm) {
    // number of processes in the comm
    int np = 0 ;
    MPI_Comm_size(comm, &np) ;
    // first we will pack the position into an vector
    std::vector<Loci::vector3d<float> > coords(len) ;
    ForwardIterator aux = first ;
    for(size_t i=0;i<len;++i,++aux) {
      const Loci::vector3d<double>& pos = aux->get_position() ;
      coords[i].x = static_cast<float>(pos.x) ;
      coords[i].y = static_cast<float>(pos.y) ;
      coords[i].z = static_cast<float>(pos.z) ;
    }
    // this is used to hold item to process mapping
    std::vector<int> i2p(len) ;
    // perform an orb partition on the coordinates
    Loci::ORBPartition(coords, i2p, comm) ;
    // then we have the item to process mapping, we
    // need to package and communicate them
    // first we need to communicate the send/recv size

    // send/recv items are the items to be send/recv
    // they are not the actual bytes to send/recv,
    // they are useful for unpacking in the final stage
    std::vector<int> send_items(np, 0) ;
    for(size_t i=0;i<len;++i)
      send_items[i2p[i]]++ ;

    std::vector<int> recv_items(np, 0) ;
    MPI_Alltoall(&send_items[0], 1, MPI_INT,
                 &recv_items[0], 1, MPI_INT, comm) ;

    // now we need to obtain the send/recv size (the real byte size)
    std::vector<int> send_sizes(np, 0) ;
    aux = first ;
    for(size_t i=0;i<len;++i,++aux)
      send_sizes[i2p[i]] += aux->pack_size() ;

    std::vector<int> recv_sizes(np, 0) ;
    MPI_Alltoall(&send_sizes[0], 1, MPI_INT,
                 &recv_sizes[0], 1, MPI_INT, comm) ;

    std::vector<int> send_displs(np, 0) ;
    std::vector<int> recv_displs(np, 0) ;
    for(int i=1;i<np;++i) {
      send_displs[i] = send_displs[i-1] + send_sizes[i-1] ;
      recv_displs[i] = recv_displs[i-1] + recv_sizes[i-1] ;
    }

    int tot_send_size = 0 ;
    int tot_recv_size = 0 ;
    int tot_recv_items = 0 ;
    for(int i=0;i<np;++i) {
      tot_send_size += send_sizes[i] ;
      tot_recv_size += recv_sizes[i] ;
      tot_recv_items += recv_items[i] ;
    }
    // allocate the send buffer
    unsigned char* send_buffer = new unsigned char[tot_send_size] ;
    // pack the send buffer
    std::vector<int> buffer_idx = send_displs ;
    aux = first ;
    for(size_t i=0;i<len;++i,++aux) {
      int pid = i2p[i] ;
      size_t po = static_cast<size_t>(buffer_idx[pid]) ;
      aux->pack(send_buffer, po) ;
      buffer_idx[pid] = po ;
    }
    // allocate the recv buffer
    unsigned char* recv_buffer = new unsigned char[tot_recv_size] ;
    // communicate the buffer
    MPI_Alltoallv(&send_buffer[0], &send_sizes[0],
                  &send_displs[0], MPI_BYTE,
                  &recv_buffer[0], &recv_sizes[0],
                  &recv_displs[0], MPI_BYTE, comm) ;

    delete[] send_buffer ;
    // then we need to extract the results
    size_t position = 0 ;
    typename std::iterator_traits<ForwardIterator>::value_type tmp ;
    for(int i=0;i<tot_recv_items;++i,++result) {
      tmp.unpack(recv_buffer, position) ;
      *result = tmp ;
    }
    
    delete[] recv_buffer ;

    return result ;
  }
  
  // this function partitions a sequence using the ORB
  // (orthogonal recursive bisection) method. the requirements
  // for the object in the sequence is that they provide
  // "pack", "unpack", "pack_size", and "get_position"
  // methods.
  template <class ForwardIterator, class OutputIterator> OutputIterator
  orb_partition_sequence(ForwardIterator first, ForwardIterator last,
                         OutputIterator result, MPI_Comm comm) {
    // get the sequence length
    unsigned int len = std::distance(first, last) ;
    // get the max sequence size in the comm
    unsigned int max_len = 0 ;
    MPI_Allreduce(&len, &max_len, 1, MPI_UNSIGNED, MPI_MAX, comm) ;

    if(max_len == 0)
      return result ;           // no elements at all

    // number of processes in the comm
    int np = 0 ;
    MPI_Comm_size(comm, &np) ;

    // get the minimum sequence length in the comm
    unsigned int min_len = 0 ;
    MPI_Allreduce(&len, &min_len, 1, MPI_UNSIGNED, MPI_MIN, comm) ;

    if(min_len < static_cast<unsigned int>(np)) {
      // if some processes have less than "np" items, perform a
      // balance first
      typedef typename std::iterator_traits<ForwardIterator>
        ::value_type value_type ;
      std::vector<value_type> new_sequence ;

      std::back_insert_iterator<std::vector<value_type> >
        bii(new_sequence) ;
      balance_sequence(first, last, bii, comm) ;
      
      len = new_sequence.size() ;
      // check to see if the minimum sequence length on
      // any process is still less than "np" or not.
      // in case of yes, just return without doing
      // anything, since we cannot perform the orb
      // partition because it requires that the items
      // on any process must be >= "np"
      MPI_Allreduce(&len, &min_len, 1, MPI_UNSIGNED, MPI_MIN, comm) ;

      if(min_len < static_cast<unsigned int>(np)) {
        // since we cannot perform orb partition, we
        // just copy the original sequence to the result
        for(typename std::vector<value_type>::const_iterator
              b=new_sequence.begin(),e=new_sequence.end();
            b!=e;++b,++result)
          *result = *b ;
        return result ;
      } else                      // else we orb partition the new_sequence
        return orb_partition_sequence_helper
          (new_sequence.begin(),
           new_sequence.end(), new_sequence.size(), result, comm) ;
    }
    // otherwise just perform the orb partition on the original sequence
    return orb_partition_sequence_helper
      (first, last, len, result, comm) ;
  }

#ifdef TO_BE_REMOVED
  // the functions in this section are obsoleted, we should remove
  // them in the future revision

  
  // reduction on distributed store and storeVec
  
  // first we define a routine that does reduction on stores
  // given a store and the initial entity domain partition,
  // this routine redistributes the stores on each process
  // such that the stores match the processes' local entitySet
  // partition. after the redistribution, each process then
  // combines the store together to produce a final store.
  // NOTE: that the passed in store "s" is assumed in the
  // dynamic temporary local numbering
  // "result_dom" is the domain for the resulting store
  // (it is in every process's local numbering)
  template<typename T> Loci::storeRepP
  redistribute_reduce_store(const store<T>& s,
                            const Map& dyn_l2g,
                            const dMap& dyn_g2l,
                            const entitySet& result_dom,
                            Loci::fact_db* factsP) {
    // we will not need to do anything for the sequential case
    if(Loci::MPI_processes == 1)
      return s.Rep() ;
        
    // get some distribute information
    Loci::fact_db::distribute_infoP df = factsP->get_distribute_info() ;
    // we will need to create a vector of domain partition based
    // on the "result_dom" supplied so that everyone knows how to
    // split their domain.
    // first we convert "result_dom" to the global numbering
    // so that everyone understands it
    entitySet result_dom_g = remap_entitySet(result_dom, df->l2g) ;
    // then gather them in a vector
    std::vector<entitySet> dom_ptn = gather_all_entitySet(result_dom_g) ;
    // partition domain
    entitySet dom_l = s.domain() ;
    // remap to global numbering
    entitySet dom = remap_entitySet(dom_l, dyn_l2g) ;
    std::vector<entitySet> dom_dist(Loci::MPI_processes) ;
    for(int i=0;i<Loci::MPI_processes;++i)
      dom_dist[i] = dom & dom_ptn[i] ;
    // then at this point, we can start to pack contents
    // for communication, the dom_dist vector is the thing
    // to send to others, but for packing purpose, we need
    // to obtain a dynamic local numbering for dom_dist
    std::vector<entitySet> dom_dist_dl(Loci::MPI_processes) ;
    for(int i=0;i<Loci::MPI_processes;++i)
      dom_dist_dl[i] = remap_entitySet(dom_dist[i], dyn_g2l) ;

    // then we proceed to compute the proper size parameters
    Loci::storeRepP srp = s.Rep() ;
    std::vector<int> send_counts(Loci::MPI_processes) ;
    for(int i=0;i<Loci::MPI_processes;++i)
      send_counts[i] = srp->pack_size(dom_dist_dl[i]) ;

    std::vector<int> send_displs(Loci::MPI_processes) ;
    send_displs[0] = 0 ;
    for(int i=1;i<Loci::MPI_processes;++i)
      send_displs[i] = send_displs[i-1] + send_counts[i-1] ;

    std::vector<int> recv_counts(Loci::MPI_processes) ;
    MPI_Alltoall(&send_counts[0], 1, MPI_INT,
                 &recv_counts[0], 1, MPI_INT, MPI_COMM_WORLD) ;

    std::vector<int> recv_displs(Loci::MPI_processes) ;
    recv_displs[0] = 0 ;
    for(int i=1;i<Loci::MPI_processes;++i)
      recv_displs[i] = recv_displs[i-1] + recv_counts[i-1] ;

    // we now have the info. to allocate send buffer
    int send_buf_size = send_displs[Loci::MPI_processes-1]
      + send_counts[Loci::MPI_processes-1] ;

    // allocate the buffer
    std::vector<unsigned char> send_buf(send_buf_size) ;
    std::vector<unsigned char*> send_ptr(Loci::MPI_processes) ;
    send_ptr[0] = &send_buf[0] ;
    for(int i=1;i<Loci::MPI_processes;++i)
      send_ptr[i] = send_ptr[i-1] + send_counts[i-1] ;

    // pack the buffer
    for(int i=0;i<Loci::MPI_processes;++i) {
      int position = 0 ;
      srp->pack(send_ptr[i], position, send_counts[i], dom_dist_dl[i]) ;
    }

    // before receiving, we need to communicate the packing sequence
    // first so that the unpack can be done in the right order
    std::vector<sequence> recv_seq(Loci::MPI_processes) ;
    for(int i=0;i<Loci::MPI_processes;++i) {
      recv_seq[i] = sequence(dom_dist_dl[i]) ;
      // renumber it to global numbering so that everyone understands
      recv_seq[i] = remap_sequence(recv_seq[i], dyn_l2g) ;
    }
    // transpose it so that everyone knows what to receive
    recv_seq = transpose_vector_sequence(recv_seq) ;

    // finally allocate the receive buffer
    int recv_buf_size = recv_displs[Loci::MPI_processes-1]
      + recv_counts[Loci::MPI_processes-1] ;

    std::vector<unsigned char> recv_buf(recv_buf_size) ;

    // final communicate step to put everything in the recv_buf
    MPI_Alltoallv(&send_buf[0], &send_counts[0],
                  &send_displs[0], MPI_PACKED,
                  &recv_buf[0], &recv_counts[0],
                  &recv_displs[0], MPI_PACKED, MPI_COMM_WORLD) ;
    // these buffer can now go away
    std::vector<unsigned char*>().swap(send_ptr) ;
    std::vector<int>().swap(send_counts) ;
    std::vector<int>().swap(send_displs) ;
    std::vector<unsigned char>().swap(send_buf) ;

    // unpack the receive buffer and perform the reduction
    // operation to combine all the stores together (currently
    // a summation operation)
    // first allocate a resulting store
    store<T> ns ;
    ns.allocate(result_dom) ;

    // we also need a temporary store for extraction
    // of the receive buffer use
    store<T> tmp ;
    tmp.allocate(result_dom) ;
    Loci::storeRepP tmp_rep = tmp.Rep() ;

    // set the new and tmp store to the default value
    for(entitySet::const_iterator ei=result_dom.begin();
        ei!=result_dom.end();++ei) {
      ns[*ei] = T() ;
      tmp[*ei] = T() ;
    }

    // extract and combine the receive buffer
    std::vector<unsigned char*> recv_ptr(Loci::MPI_processes) ;
    recv_ptr[0] = &recv_buf[0] ;
    for(int i=1;i<Loci::MPI_processes;++i)
      recv_ptr[i] = recv_ptr[i-1] + recv_counts[i-1] ;

    for(int i=0;i<Loci::MPI_processes;++i) {
      // get the unpack sequence in local numbering
      sequence seq = remap_sequence(recv_seq[i], df->g2l) ;
      int position = 0 ;
      // unpack to tmp store
      tmp_rep->unpack(recv_ptr[i], position, recv_counts[i], seq) ;
      // then combine tmp store to the resulting new store
      for(sequence::const_iterator si=seq.begin();si!=seq.end();++si)
        ns[*si] += tmp[*si] ;
    }

    // return the new combined store
    return ns.Rep() ;  
  }

  // this is the storeVec version
  template<typename T> Loci::storeRepP
  redistribute_reduce_storeVec(const storeVec<T>& s,
                               const Map& dyn_l2g,
                               const dMap& dyn_g2l,
                               const entitySet& result_dom,
                               Loci::fact_db* factsP) {
    // we will not need to do anything for the sequential case
    if(Loci::MPI_processes == 1)
      return s.Rep() ;
        
    // get some distribute information
    Loci::fact_db::distribute_infoP df = factsP->get_distribute_info() ;
    // we will need to create a vector of domain partition based
    // on the "result_dom" supplied so that everyone knows how to
    // split their domain.
    // first we convert "result_dom" to the global numbering
    // so that everyone understands it
    entitySet result_dom_g = remap_entitySet(result_dom, df->l2g) ;
    // then gather them in a vector
    std::vector<entitySet> dom_ptn = gather_all_entitySet(result_dom_g) ;
    // partition domain
    entitySet dom_l = s.domain() ;
    // remap to global numbering
    entitySet dom = remap_entitySet(dom_l, dyn_l2g) ;
    std::vector<entitySet> dom_dist(Loci::MPI_processes) ;
    for(int i=0;i<Loci::MPI_processes;++i)
      dom_dist[i] = dom & dom_ptn[i] ;
    // then at this point, we can start to pack contents
    // for communication, the dom_dist vector is the thing
    // to send to others, but for packing purpose, we need
    // to obtain a dynamic local numbering for dom_dist
    std::vector<entitySet> dom_dist_dl(Loci::MPI_processes) ;
    for(int i=0;i<Loci::MPI_processes;++i)
      dom_dist_dl[i] = remap_entitySet(dom_dist[i], dyn_g2l) ;

    // then we proceed to compute the proper size parameters
    Loci::storeRepP srp = s.Rep() ;
    std::vector<int> send_counts(Loci::MPI_processes) ;
    for(int i=0;i<Loci::MPI_processes;++i)
      send_counts[i] = srp->pack_size(dom_dist_dl[i]) ;

    std::vector<int> send_displs(Loci::MPI_processes) ;
    send_displs[0] = 0 ;
    for(int i=1;i<Loci::MPI_processes;++i)
      send_displs[i] = send_displs[i-1] + send_counts[i-1] ;

    std::vector<int> recv_counts(Loci::MPI_processes) ;
    MPI_Alltoall(&send_counts[0], 1, MPI_INT,
                 &recv_counts[0], 1, MPI_INT, MPI_COMM_WORLD) ;

    std::vector<int> recv_displs(Loci::MPI_processes) ;
    recv_displs[0] = 0 ;
    for(int i=1;i<Loci::MPI_processes;++i)
      recv_displs[i] = recv_displs[i-1] + recv_counts[i-1] ;

    // we now have the info. to allocate send buffer
    int send_buf_size = send_displs[Loci::MPI_processes-1]
      + send_counts[Loci::MPI_processes-1] ;

    // allocate the buffer
    std::vector<unsigned char> send_buf(send_buf_size) ;
    std::vector<unsigned char*> send_ptr(Loci::MPI_processes) ;
    send_ptr[0] = &send_buf[0] ;
    for(int i=1;i<Loci::MPI_processes;++i)
      send_ptr[i] = send_ptr[i-1] + send_counts[i-1] ;

    // pack the buffer
    for(int i=0;i<Loci::MPI_processes;++i) {
      int position = 0 ;
      srp->pack(send_ptr[i], position, send_counts[i], dom_dist_dl[i]) ;
    }

    // before receiving, we need to communicate the packing sequence
    // first so that the unpack can be done in the right order
    std::vector<sequence> recv_seq(Loci::MPI_processes) ;
    for(int i=0;i<Loci::MPI_processes;++i) {
      recv_seq[i] = sequence(dom_dist_dl[i]) ;
      // renumber it to global numbering so that everyone understands
      recv_seq[i] = remap_sequence(recv_seq[i], dyn_l2g) ;
    }
    // transpose it so that everyone knows what to receive
    recv_seq = transpose_vector_sequence(recv_seq) ;

    // finally allocate the receive buffer
    int recv_buf_size = recv_displs[Loci::MPI_processes-1]
      + recv_counts[Loci::MPI_processes-1] ;

    std::vector<unsigned char> recv_buf(recv_buf_size) ;

    // final communicate step to put everything in the recv_buf
    MPI_Alltoallv(&send_buf[0], &send_counts[0],
                  &send_displs[0], MPI_PACKED,
                  &recv_buf[0], &recv_counts[0],
                  &recv_displs[0], MPI_PACKED, MPI_COMM_WORLD) ;
    // these buffer can now go away
    std::vector<unsigned char*>().swap(send_ptr) ;
    std::vector<int>().swap(send_counts) ;
    std::vector<int>().swap(send_displs) ;
    std::vector<unsigned char>().swap(send_buf) ;

    // unpack the receive buffer and perform the reduction
    // operation to combine all the stores together (currently
    // a summation operation)
    // first allocate a resulting store
    int vec_size = s.vecSize() ;
    storeVec<T> ns ;
    ns.setVecSize(vec_size) ;
    ns.allocate(result_dom) ;

    // we also need a temporary store for extraction
    // of the receive buffer use
    storeVec<T> tmp ;
    tmp.setVecSize(vec_size) ;
    tmp.allocate(result_dom) ;
    Loci::storeRepP tmp_rep = tmp.Rep() ;

    // set the new and tmp store to the default value
    for(entitySet::const_iterator ei=result_dom.begin();
        ei!=result_dom.end();++ei) {
      for(int i=0;i<vec_size;++i) {
        ns[*ei][i] = T() ;
        tmp[*ei][i] = T() ;
      }
    }

    // extract and combine the receive buffer
    std::vector<unsigned char*> recv_ptr(Loci::MPI_processes) ;
    recv_ptr[0] = &recv_buf[0] ;
    for(int i=1;i<Loci::MPI_processes;++i)
      recv_ptr[i] = recv_ptr[i-1] + recv_counts[i-1] ;

    for(int i=0;i<Loci::MPI_processes;++i) {
      // get the unpack sequence in local numbering
      sequence seq = remap_sequence(recv_seq[i], df->g2l) ;
      int position = 0 ;
      // unpack to tmp store
      tmp_rep->unpack(recv_ptr[i], position, recv_counts[i], seq) ;
      // then combine tmp store to the resulting new store
      for(sequence::const_iterator si=seq.begin();si!=seq.end();++si) {
        for(int i=0;i<vec_size;++i)
          ns[*si][i] += tmp[*si][i] ;
      }
    }

    // return the new combined store
    return ns.Rep() ;  
  }
#endif

  // this function reads in data from an hdf5 file into a vector
  template<typename T> void
  readUnorderedVector(hid_t group_id,
                      const char* element_name, std::vector<T>& v) {
    if(Loci::MPI_rank == 0) {
      // first figure out the size that each process will receive
      hid_t dataset = H5Dopen(group_id, element_name) ;
      hid_t dataspace = H5Dget_space(dataset);
      
      typedef data_schema_traits<T> traits_type ;
      Loci::DatatypeP dp = traits_type::get_type() ;

      hsize_t rank = H5Sget_simple_extent_ndims(dataspace);
      hsize_t len = 0 ;
      H5Sget_simple_extent_dims(dataspace, &len, NULL);

      if(rank != 1) {
        std::cerr << "Error: readUnorderedVector hdf5 dataspace rank != 1"
             << ", Aborting..." << std::endl ;
        Loci::Abort() ;
      }

      // broadcast the total length of the data
      // Note: this only works on homogeneous systems
      // because the byte-order has to be the same
      int hsize_t_size = sizeof(hsize_t) ;
      MPI_Bcast(&len,hsize_t_size,MPI_BYTE,0,MPI_COMM_WORLD) ;

      hsize_t my_len = len / Loci::MPI_processes ;
      if( (hsize_t)Loci::MPI_rank < len % Loci::MPI_processes)
        my_len += 1 ;

      // rank 0 reads in data consecutively and sends it
      // to other processes.

#ifdef H5_INTERFACE_1_6_4
      hsize_t offset = 0 ;
#else
      hssize_t offset = 0 ;
#endif
      for(int i=1;i<Loci::MPI_processes;++i) {
        hsize_t ilen = len / Loci::MPI_processes ;
        if( (hsize_t)i < len % Loci::MPI_processes)
          ilen += 1 ;
        if(ilen > 0) {
          // create a hyperslab in the hdf5 dataset
          H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                              &offset, NULL, &ilen, NULL) ;
          // create a memory space
          hid_t memspace = H5Screate_simple(rank, &ilen, NULL) ;
          // create a memory hyperslab
#ifdef H5_INTERFACE_1_6_4
          hsize_t mem_offset = 0 ;
#else
          hssize_t mem_offset = 0 ;
#endif
          H5Sselect_hyperslab(memspace, H5S_SELECT_SET,
                              &mem_offset, NULL, &ilen, NULL);
          // define a buffer to read data from the hdf5 dataspace
          std::vector<T> buffer(ilen) ;
          // read the hyperslab
          H5Dread(dataset, dp->get_hdf5_type(),
                  memspace, dataspace, H5P_DEFAULT, &buffer[0]) ;
          
          H5Sclose(memspace);
          // then send it to process i
          MPI_Send(&buffer[0], sizeof(T)*ilen,
                   MPI_BYTE, i, 0, MPI_COMM_WORLD) ;
        }
        offset += ilen ;
      }

      // reads the hyperslab data for rank 0
      if(my_len > 0) {
        H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                            &offset, NULL, &my_len, NULL) ;
        hid_t memspace = H5Screate_simple(rank, &my_len, NULL) ;
#ifdef H5_INTERFACE_1_6_4
        hsize_t mem_offset = 0 ;
#else
        hssize_t mem_offset = 0 ;
#endif
        H5Sselect_hyperslab(memspace, H5S_SELECT_SET,
                            &mem_offset, NULL, &my_len, NULL);
        v.resize(my_len) ;
        H5Dread(dataset, dp->get_hdf5_type(),
                memspace, dataspace, H5P_DEFAULT, &v[0]) ;
        H5Sclose(memspace);
      }

      H5Dclose(dataset);
      H5Sclose(dataspace);
      
    } else {
      // first wait to receive the total data size
      hsize_t len = 0 ;
      int hsize_t_size = sizeof(hsize_t) ;
      MPI_Bcast(&len,hsize_t_size,MPI_BYTE,0,MPI_COMM_WORLD) ;
      int my_len = len / Loci::MPI_processes ;
      if( (unsigned int)Loci::MPI_rank < len % Loci::MPI_processes)
        my_len += 1 ;

      // receive messages from rank 0
      if(my_len > 0) {
        v.resize(my_len) ;
        MPI_Recv(&v[0], sizeof(T)*my_len,
                 MPI_BYTE, 0, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE) ;
      }
    }
  }

  // a structure for MPI_2INT
  struct mpi_2int {
    int first, second ;
    mpi_2int(int f=0,int s=std::numeric_limits<int>::min())
      :first(f), second(s) {}
  } ;
  // a user function defined intended for reduction on
  // MPI_2INT (the first being a MPI_SUM, the second MPI_MAX)
  void
  SUM_MAX_2INT(mpi_2int *invec, mpi_2int *inoutvec,
               int *len, MPI_Datatype *dtype) ;

  // functions that computes median value for a group of processes
  int
  mpi_median(int value, MPI_Comm comm) ;
  long
  mpi_median(long value, MPI_Comm comm) ;
  double
  mpi_median(double value, MPI_Comm comm) ;


  ///////////////////////////////////////////////////
  // these are the codes for hilbert curve partition
  ///////////////////////////////////////////////////

  // Hilbert Key
  typedef Loci::Array<unsigned int, 3> HilbertCode ;
  // Integer coordinate
  typedef Loci::Array<unsigned int, 3> IntCoord3 ;

  template<typename T>
  struct HilbertKey {
    HilbertCode key ;
    T obj ;

    size_t
    pack_size() const {
      return sizeof(HilbertKey) ;
    }
    
    void
    pack(unsigned char* buffer, size_t& position) const {
      HilbertKey* b = reinterpret_cast<HilbertKey*>(buffer+position) ;
      *b = *this ;
      position += pack_size() ;
    }

    void
    unpack(unsigned char* buffer, size_t& position) {
      HilbertKey* b = reinterpret_cast<HilbertKey*>(buffer+position) ;
      *this = *b ;
      position += pack_size() ;
    }
  } ;
  template<typename T> inline bool
  operator<(const HilbertKey<T>& k1, const HilbertKey<T>& k2) {
    return
      ( (k1.key[2]<k2.key[2]) ||
        (k1.key[2]==k2.key[2]&&k1.key[1]<k2.key[1]) ||
        (k1.key[2]==k2.key[2]&&k1.key[1]==k2.key[1]&&k1.key[0]<k2.key[0])) ;
  }
  // this function generates a hilbert index for a 3D integer coordinate
  HilbertCode hilbert_encode(const IntCoord3& p) ;

  // a helper function used for the "hilbert_partition_sequence"
  // it implements the main communication routines, but
  // assumes that the sequence passed in is indeed valid
  // for performing the hilbert partition. (len >= np) the addtional
  // parameter sequence_len indicates the sequence length (to avoid
  // recompute it again.
  template <class ForwardIterator, class OutputIterator> OutputIterator
  hilbert_partition_sequence_helper
  (ForwardIterator first, ForwardIterator last,
   size_t len, OutputIterator result, MPI_Comm comm) {
    // number of processes in the comm
    int np = 0 ;
    MPI_Comm_size(comm, &np) ;

    // first we'll compute an integer coordinate for each object
    vector3d<double> max_pos, min_pos ;
    ForwardIterator b = first ;

    max_pos = b->get_position() ;
    min_pos = max_pos ;
    
    for(++b;b!=last;++b) {
      const vector3d<double>& pos = b->get_position() ;
      max_pos = vector3d<double>(std::max(max_pos.x, pos.x),
                                 std::max(max_pos.y, pos.y),
                                 std::max(max_pos.z, pos.z)) ;
      min_pos = vector3d<double>(std::min(min_pos.x, pos.x),
                                 std::min(min_pos.y, pos.y),
                                 std::min(min_pos.z, pos.z)) ;
    }
    // get the global max and min
    vector3d<double> global_max_pos, global_min_pos ;
    MPI_Allreduce(&max_pos.x, &global_max_pos.x,
                  3, MPI_DOUBLE, MPI_MAX, comm) ;
    MPI_Allreduce(&min_pos.x, &global_min_pos.x,
                  3, MPI_DOUBLE, MPI_MIN, comm) ;
    
    typedef typename std::iterator_traits<ForwardIterator>
      ::value_type OBJ ;

    std::vector<HilbertKey<OBJ> > key_list(len) ;

    vector3d<double> s = global_max_pos - global_min_pos ;
    double scale = 4e9 / std::max(s.x, std::max(s.y,s.z)) ;
    
    int idx = 0 ;
    for(b=first;b!=last;++b,++idx) {
      IntCoord3 p ;
      vector3d<double> p_base = scale*(b->get_position() - global_min_pos) ;
      p[0] = (unsigned int)(p_base.x) ;
      p[1] = (unsigned int)(p_base.y) ;
      p[2] = (unsigned int)(p_base.z) ;

      key_list[idx].key = hilbert_encode(p) ;
      key_list[idx].obj = *b ;

    }

    Loci::parSampleSort(key_list, comm) ;



    // finally balance the sorted list of keys
    std::vector<HilbertKey<OBJ> > bal_key_list ;
    std::back_insert_iterator<std::vector<HilbertKey<OBJ> > >
      bii(bal_key_list) ;


    balance_sequence(key_list.begin(), key_list.end(), bii, comm) ;

    std::vector<HilbertKey<OBJ> >().swap(key_list) ;

    // we've done, just push the results to the result iterator
    for(size_t i=0;i<bal_key_list.size();++i,++result) {
      *result = bal_key_list[i].obj ;
    }

    return result ;
  }
  

  // this function partitions a sequence using the Hilbert curve
  // mapping. the requirements for the object in the sequence is
  // that they provide "pack", "unpack", "pack_size", and
  // "get_position" methods.
  template <class ForwardIterator, class OutputIterator> OutputIterator
  hilbert_partition_sequence(ForwardIterator first, ForwardIterator last,
                             OutputIterator result, MPI_Comm comm) {
    // we'll first balance the list
    
    // get the sequence length
    unsigned int len = std::distance(first, last) ;
    // get the max sequence size in the comm
    unsigned int max_len = 0 ;
    MPI_Allreduce(&len, &max_len, 1, MPI_UNSIGNED, MPI_MAX, comm) ;

    if(max_len == 0)
      return result ;           // no elements at all

    // number of processes in the comm
    int np = 0 ;
    MPI_Comm_size(comm, &np) ;

    // get the minimum sequence length in the comm
    unsigned int min_len = 0 ;
    MPI_Allreduce(&len, &min_len, 1, MPI_UNSIGNED, MPI_MIN, comm) ;

    if(min_len < static_cast<unsigned int>(np)) {
      // if some processes have less than "np" items, perform a
      // balance first
      typedef typename std::iterator_traits<ForwardIterator>
        ::value_type value_type ;
      std::vector<value_type> new_sequence ;

      std::back_insert_iterator<std::vector<value_type> >
        bii(new_sequence) ;

      balance_sequence(first, last, bii, comm) ;
      
      len = new_sequence.size() ;
      // check to see if the minimum sequence length on
      // any process is still less than "np" or not.
      // in case of yes, just return without doing
      // anything, since we cannot perform the hilbert
      // partition because it requires that the items
      // on any process must be >= "np" (due to the parallel sample sort)
      MPI_Allreduce(&len, &min_len, 1, MPI_UNSIGNED, MPI_MIN, comm) ;

      if(min_len < static_cast<unsigned int>(np)) {
        // since we cannot perform hilbert partition, we
        // just copy the original sequence to the result
        for(typename std::vector<value_type>::const_iterator
              b=new_sequence.begin(),e=new_sequence.end();
            b!=e;++b,++result)
          *result = *b ;
        return result ;
      } else                      // else we partition the new_sequence
        return hilbert_partition_sequence_helper
          (new_sequence.begin(),
           new_sequence.end(), new_sequence.size(), result, comm) ;
    }
    // otherwise just perform the hilbert partition on the original sequence
    return hilbert_partition_sequence_helper
      (first, last, len, result, comm) ;
  }

  // this function sorts a vector according to hilbert index
  template <class ParticleType> void
  hilbert_sort_vector(std::vector<ParticleType>& vp, MPI_Comm comm) {
    unsigned int len = vp.size() ;
    unsigned int max_len = 0 ;
    MPI_Allreduce(&len, &max_len, 1, MPI_UNSIGNED, MPI_MAX, comm) ;

    if(max_len == 0)
      return ;                  // do not sort

    int np = 0 ;
    MPI_Comm_size(comm, &np) ;

    unsigned int min_len = 0 ;
    MPI_Allreduce(&len, &min_len, 1, MPI_UNSIGNED, MPI_MIN, comm) ;

    if(min_len < static_cast<unsigned int>(np)) {
      // balance the vector first since the min partition
      // has a length less than the number of processes
      Loci::balanceDistribution(vp, comm) ;

      len = vp.size() ;
      MPI_Allreduce(&len, &min_len, 1, MPI_UNSIGNED, MPI_MIN, comm) ;

      if(min_len < static_cast<unsigned int>(np)) {
        // if the minimum partition length is still less than
        // the number of parallel processes, we won't sort, just quit
        return ;
      }
    }
    // go ahead and compute the hilbert index for each item in "vp"
    vector3d<double> max_pos, min_pos ;
    typename std::vector<ParticleType>::iterator b = vp.begin() ;
    max_pos = b->get_position() ;
    min_pos = max_pos ;

    for(++b;b!=vp.end();++b) {
      const vector3d<double>& pos = b->get_position() ;
      max_pos = vector3d<double>(std::max(max_pos.x, pos.x),
                                 std::max(max_pos.y, pos.y),
                                 std::max(max_pos.z, pos.z)) ;
      min_pos = vector3d<double>(std::min(min_pos.x, pos.x),
                                 std::min(min_pos.y, pos.y),
                                 std::min(min_pos.z, pos.z)) ;
    }
    vector3d<double> global_max_pos, global_min_pos ;
    MPI_Allreduce(&max_pos.x, &global_max_pos.x,
                  3, MPI_DOUBLE, MPI_MAX, comm) ;
    MPI_Allreduce(&min_pos.x, &global_min_pos.x,
                  3, MPI_DOUBLE, MPI_MIN, comm) ;

    std::vector<HilbertKey<ParticleType> > key_list(len) ;

    vector3d<double> s = global_max_pos - global_min_pos ;
    double scale = 4e9 / std::max(s.x, std::max(s.y,s.z)) ;

    int idx = 0 ;
    for(b=vp.begin();b!=vp.end();++b,++idx) {
      IntCoord3 p ;
      vector3d<double> p_base = scale*(b->get_position() - global_min_pos) ;
      p[0] = (unsigned int)(p_base.x) ;
      p[1] = (unsigned int)(p_base.y) ;
      p[2] = (unsigned int)(p_base.z) ;
      key_list[idx].key = hilbert_encode(p) ;
      key_list[idx].obj = *b ;
    }

    // sort the vector
    Loci::parSampleSort(key_list, comm) ;

    vp.clear() ;

    for(typename std::vector<HilbertKey<ParticleType> >::const_iterator
          vi=key_list.begin();vi!=key_list.end();++vi)
      vp.push_back(vi->obj) ;
  }


} // end of namespace lagrangianP

#endif
