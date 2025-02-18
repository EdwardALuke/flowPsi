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

// implementation of "par_util.h"
#include "particle_config.h"
#include "par_util.h"
#include <parSampleSort.h>
#include <algorithm>

using std::cout ;
using std::cerr ;
using std::endl ;
using std::vector ;
using std::list ;
#include <Tools/block_hash.h>
using Loci::block_hash ;

namespace Loci {
  extern void ORBPartition(const vector<vector3d<float> >&,
                           vector<int>&, MPI_Comm) ;
}

namespace lagrangianP {
  // this function transposes the passed in vector<entitySet>
  // by an all to all personalized communication
  vector<entitySet>
  transpose_vector_entitySet(const vector<entitySet>& in) {
    // first compute the send count and displacement
    int* send_counts = new int[Loci::MPI_processes] ;
    for(int i=0;i<Loci::MPI_processes;++i)
      send_counts[i] = in[i].size() ;
    int* send_displs = new int[Loci::MPI_processes] ;
    send_displs[0] = 0 ;
    for(int i=1;i<Loci::MPI_processes;++i)
      send_displs[i] = send_displs[i-1] + send_counts[i-1] ;
    // then communicate this get the recv info.
    int* recv_counts = new int[Loci::MPI_processes] ;
    MPI_Alltoall(send_counts, 1, MPI_INT,
                 recv_counts, 1, MPI_INT, MPI_COMM_WORLD) ;
    int* recv_displs = new int[Loci::MPI_processes] ;
    recv_displs[0] = 0 ;
    for(int i=1;i<Loci::MPI_processes;++i)
      recv_displs[i] = recv_displs[i-1] + recv_counts[i-1] ;
    // all info. gathered, ready to do MPI_Alltoallv
    // first pack data into a raw buffer.
    int buf_size = 0 ;
    for(int i=0;i<Loci::MPI_processes;++i)
      buf_size += send_counts[i] ;
    int* send_buf = new int[buf_size] ;
    int buf_idx = 0 ;
    for(int i=0;i<Loci::MPI_processes;++i) {
      const entitySet& eset = in[i] ;
      for(entitySet::const_iterator ei=eset.begin();
          ei!=eset.end();++ei,++buf_idx)
        send_buf[buf_idx] = *ei ;
    }
    // allocate receive buffer
    int recv_size = 0 ;
    for(int i=0;i<Loci::MPI_processes;++i)
      recv_size += recv_counts[i] ;
    int* recv_buf = new int[recv_size] ;
    // communicate
    MPI_Alltoallv(send_buf, send_counts,
                  send_displs, MPI_INT,
                  recv_buf, recv_counts,
                  recv_displs, MPI_INT, MPI_COMM_WORLD) ;
    delete[] send_counts ;
    delete[] send_displs ;
    delete[] recv_counts ;
    delete[] send_buf ;
    // unpack recv buffer into a vector of entitySet
    vector<entitySet> out(Loci::MPI_processes) ;    
    int k = 0 ;
    for(int i=0;i<Loci::MPI_processes;++i) {
      int limit ;
      if(i == Loci::MPI_processes-1)
        limit = recv_size ;
      else
        limit = recv_displs[i+1] ;
      for(;k<limit;++k)
        out[i] += recv_buf[k] ;
    }
    delete[] recv_displs ;
    delete[] recv_buf ;

    return out ;
  }

  void
  transpose_vector_entitySet_opt(const vector<entitySet>& in,
                                 vector<entitySet>& out) {
    vector<bool> pack_interval(Loci::MPI_processes,false) ;
    // first compute the send count and displacement
    vector<int> send_counts(Loci::MPI_processes,0) ;
    for(int i=0;i<Loci::MPI_processes;++i) {
      // we don't pack empty entitySet
      if(in[i] == EMPTY) {
        send_counts[i] = 0 ;
        continue ;
      }
      // since if sending intervals, the send size will
      // be 2 * the number of intervals, then if that is
      // larger than the total element size, then we choose
      // to communicate elements directly, otherwise, we
      // choose to send the intervals
      int num_intervals = in[i].num_intervals() ;
      num_intervals*=2 ;
      int size = in[i].size() ;
      if(num_intervals >= size) {
        pack_interval[i] = false ;
        // +1 since we include a flag in the head to indicate
        // whether this is interval, or elements packed
        send_counts[i] = size + 1 ;
      } else {
        pack_interval[i] = true ;
        send_counts[i] = num_intervals + 1 ;
      }
    }
    vector<int> send_displs(Loci::MPI_processes,0) ;
    send_displs[0] = 0 ;
    for(int i=1;i<Loci::MPI_processes;++i)
      send_displs[i] = send_displs[i-1] + send_counts[i-1] ;
    // then communicate this get the recv info.
    vector<int> recv_counts(Loci::MPI_processes,0) ;
    MPI_Alltoall(&send_counts[0], 1, MPI_INT,
                 &recv_counts[0], 1, MPI_INT, MPI_COMM_WORLD) ;
    vector<int> recv_displs(Loci::MPI_processes,0) ;
    recv_displs[0] = 0 ;
    for(int i=1;i<Loci::MPI_processes;++i)
      recv_displs[i] = recv_displs[i-1] + recv_counts[i-1] ;
    // all info. gathered, ready to do MPI_Alltoallv
    // first pack data into a raw buffer.
    int buf_size = send_counts[Loci::MPI_processes-1] +
      send_displs[Loci::MPI_processes-1] ;
    
    vector<int> send_buf(buf_size) ;
    int buf_idx = 0 ;
    for(int i=0;i<Loci::MPI_processes;++i) {
      // only pack non-empty entitySet
      if(send_counts[i] == 0)
        continue ;
      
      const entitySet& eset = in[i] ;
      if(pack_interval[i]) {
        // packing intervals
        send_buf[buf_idx++] = 1 ; // set flag to indicate
                                  // that this is intervals packed
        for(size_t k=0;k<eset.num_intervals();++k) {
          send_buf[buf_idx++] = eset[k].first ;
          send_buf[buf_idx++] = eset[k].second ;
        }
      } else {
        send_buf[buf_idx++] = 0 ; // elements directly packed
        for(entitySet::const_iterator ei=eset.begin();
            ei!=eset.end();++ei,++buf_idx)
          send_buf[buf_idx] = *ei ;
      }
    }
    // allocate receive buffer
    int recv_size = recv_displs[Loci::MPI_processes-1] +
      recv_counts[Loci::MPI_processes-1] ;

    vector<int> recv_buf(recv_size) ;
    // communicate
    MPI_Alltoallv(&send_buf[0], &send_counts[0],
                  &send_displs[0], MPI_INT,
                  &recv_buf[0], &recv_counts[0],
                  &recv_displs[0], MPI_INT, MPI_COMM_WORLD) ;
    // release buffers that are not needed
    vector<int>().swap(send_counts) ;
    vector<int>().swap(send_displs) ;
    vector<int>().swap(send_buf) ;

    // unpack recv buffer into a vector of entitySet
    if(out.size() != (size_t)Loci::MPI_processes)
      out.resize(Loci::MPI_processes) ;    

    for(int i=0;i<Loci::MPI_processes;++i) {
      int b = recv_displs[i] ;
      int e = b + recv_counts[i] ;
      if(b == e)
        continue ;              // empty buffer

      entitySet& eset = out[i] ;
      int flag = recv_buf[b] ;
      ++b ;
      if(flag == 1) {
        // if packed with interval
        for(;b<e;b+=2) {
          int l = recv_buf[b] ;
          int u = recv_buf[b+1] ;
          eset += interval(l,u) ;
        }
      } else {
        // packed with elements
        for(;b<e;++b) {
          eset += recv_buf[b] ;
        }
      }
    }

    // the end...
  }
  
  // transpose a vector<sequence>
  vector<sequence>
  transpose_vector_sequence(const vector<sequence>& in) {
    // first compute the send count and displacement
    int* send_counts = new int[Loci::MPI_processes] ;
    for(int i=0;i<Loci::MPI_processes;++i)
      send_counts[i] = in[i].size() ;
    int* send_displs = new int[Loci::MPI_processes] ;
    send_displs[0] = 0 ;
    for(int i=1;i<Loci::MPI_processes;++i)
      send_displs[i] = send_displs[i-1] + send_counts[i-1] ;
    // then communicate this get the recv info.
    int* recv_counts = new int[Loci::MPI_processes] ;
    MPI_Alltoall(send_counts, 1, MPI_INT,
                 recv_counts, 1, MPI_INT, MPI_COMM_WORLD) ;
    int* recv_displs = new int[Loci::MPI_processes] ;
    recv_displs[0] = 0 ;
    for(int i=1;i<Loci::MPI_processes;++i)
      recv_displs[i] = recv_displs[i-1] + recv_counts[i-1] ;
    // all info. gathered, ready to do MPI_Alltoallv
    // first pack data into a raw buffer.
    int buf_size = 0 ;
    for(int i=0;i<Loci::MPI_processes;++i)
      buf_size += send_counts[i] ;
    int* send_buf = new int[buf_size] ;
    int buf_idx = 0 ;
    for(int i=0;i<Loci::MPI_processes;++i) {
      const sequence& seq = in[i] ;
      for(sequence::const_iterator si=seq.begin();
          si!=seq.end();++si,++buf_idx)
        send_buf[buf_idx] = *si ;
    }
    // allocate receive buffer
    int recv_size = 0 ;
    for(int i=0;i<Loci::MPI_processes;++i)
      recv_size += recv_counts[i] ;
    int* recv_buf = new int[recv_size] ;
    // communicate
    MPI_Alltoallv(send_buf, send_counts,
                  send_displs, MPI_INT,
                  recv_buf, recv_counts,
                  recv_displs, MPI_INT, MPI_COMM_WORLD) ;
    delete[] send_counts ;
    delete[] send_displs ;
    delete[] recv_counts ;
    delete[] send_buf ;
    // unpack recv buffer into a vector of entitySet
    vector<sequence> out(Loci::MPI_processes) ;    
    int k = 0 ;
    for(int i=0;i<Loci::MPI_processes;++i) {
      int limit ;
      if(i == Loci::MPI_processes-1)
        limit = recv_size ;
      else
        limit = recv_displs[i+1] ;
      for(;k<limit;++k)
        out[i] += recv_buf[k] ;
    }
    delete[] recv_displs ;
    delete[] recv_buf ;

    return out ;
  }

  void
  transpose_vector_sequence_opt(const vector<sequence>& in,
                                 vector<sequence>& out) {
    vector<bool> pack_interval(Loci::MPI_processes,false) ;
    // first compute the send count and displacement
    vector<int> send_counts(Loci::MPI_processes,0) ;
    for(int i=0;i<Loci::MPI_processes;++i) {
      if(in[i] == EMPTY) {
        send_counts[i] = 0 ;
        continue ;
      }
      // since if sending intervals, the send size will
      // be 2 * the number of intervals, then if that is
      // larger than the total element size, then we choose
      // to communicate elements directly, otherwise, we
      // choose to send the intervals
      int num_intervals = in[i].num_intervals() ;
      num_intervals*=2 ;
      int size = in[i].size() ;
      if(num_intervals >= size) {
        pack_interval[i] = false ;
        // +1 since we include a flag in the head to indicate
        // whether this is interval, or elements packed
        send_counts[i] = size + 1 ;
      } else {
        pack_interval[i] = true ;
        send_counts[i] = num_intervals + 1 ;
      }
    }
    vector<int> send_displs(Loci::MPI_processes,0) ;
    send_displs[0] = 0 ;
    for(int i=1;i<Loci::MPI_processes;++i)
      send_displs[i] = send_displs[i-1] + send_counts[i-1] ;
    // then communicate this get the recv info.
    vector<int> recv_counts(Loci::MPI_processes,0) ;
    MPI_Alltoall(&send_counts[0], 1, MPI_INT,
                 &recv_counts[0], 1, MPI_INT, MPI_COMM_WORLD) ;
    vector<int> recv_displs(Loci::MPI_processes,0) ;
    recv_displs[0] = 0 ;
    for(int i=1;i<Loci::MPI_processes;++i)
      recv_displs[i] = recv_displs[i-1] + recv_counts[i-1] ;
    // all info. gathered, ready to do MPI_Alltoallv
    // first pack data into a raw buffer.
    int buf_size = send_counts[Loci::MPI_processes-1] +
      send_displs[Loci::MPI_processes-1] ;
    
    vector<int> send_buf(buf_size) ;
    int buf_idx = 0 ;
    for(int i=0;i<Loci::MPI_processes;++i) {
      if(send_counts[i] == 0)
        continue ;              // do not pack for empty message
      
      const sequence& seq = in[i] ;
      if(pack_interval[i]) {
        // packing intervals
        send_buf[buf_idx++] = 1 ; // set flag to indicate
                                  // that this is intervals packed
        for(size_t k=0;k<seq.num_intervals();++k) {
          send_buf[buf_idx++] = seq[k].first ;
          send_buf[buf_idx++] = seq[k].second ;
        }
      } else {
        send_buf[buf_idx++] = 0 ; // elements directly packed
        for(sequence::const_iterator si=seq.begin();
            si!=seq.end();++si,++buf_idx)
          send_buf[buf_idx] = *si ;
      }
    }
    // allocate receive buffer
    int recv_size = recv_displs[Loci::MPI_processes-1] +
      recv_counts[Loci::MPI_processes-1] ;

    vector<int> recv_buf(recv_size) ;
    // communicate
    MPI_Alltoallv(&send_buf[0], &send_counts[0],
                  &send_displs[0], MPI_INT,
                  &recv_buf[0], &recv_counts[0],
                  &recv_displs[0], MPI_INT, MPI_COMM_WORLD) ;
    // release buffers that are not needed
    vector<int>().swap(send_counts) ;
    vector<int>().swap(send_displs) ;
    vector<int>().swap(send_buf) ;

    // unpack recv buffer into a vector of entitySet
    if(out.size() != (size_t)Loci::MPI_processes)
      out.resize(Loci::MPI_processes) ;    

    for(int i=0;i<Loci::MPI_processes;++i) {
      int b = recv_displs[i] ;
      int e = b + recv_counts[i] ;
      if(b==e)
        continue ;
      
      sequence& seq = out[i] ;

      int flag = recv_buf[b] ;
      ++b ;
      if(flag == 1) {
        // if packed with interval
        for(;b<e;b+=2) {
          int f = recv_buf[b] ;
          int s = recv_buf[b+1] ;
          seq += interval(f,s) ;
        }
      } else {
        // packed with elements
        for(;b<e;++b) {
          seq += recv_buf[b] ;
        }
      }
    }

    // the end...
  }

  // a utility function that takes an entitySet from a processor
  // and returns a vector of entitySet gathered from all processors
  vector<entitySet>
  gather_all_entitySet(const entitySet& eset) {
    int local_size = eset.size() ;
    int global_size = global_sum(local_size) ;
    // compute receive counts from all processors
    int* recv_counts = new int[Loci::MPI_processes] ;
    MPI_Allgather(&local_size, 1, MPI_INT,
                  recv_counts, 1, MPI_INT, MPI_COMM_WORLD) ;
    // then compute receive displacement
    int* recv_displs = new int[Loci::MPI_processes] ;
    recv_displs[0] = 0 ;
    for(int i=1;i<Loci::MPI_processes;++i)
      recv_displs[i] = recv_displs[i-1] + recv_counts[i-1] ;
    // pack the local eset into an array
    int* local_eset = new int[local_size] ;
    int count = 0 ;
    for(entitySet::const_iterator ei=eset.begin();
        ei!=eset.end();++ei,++count)
      local_eset[count] = *ei ;
    // allocate the entire array for all data from all processors
    int* global_eset = new int[global_size] ;
    // communicate to obtain all esets from every processors
    MPI_Allgatherv(local_eset, local_size, MPI_INT,
                   global_eset, recv_counts, recv_displs,
                   MPI_INT, MPI_COMM_WORLD) ;
    delete[] local_eset ;
    delete[] recv_counts ;
    // unpack the raw buffer into a vector<entitySet>
    vector<entitySet> ret(Loci::MPI_processes) ;
    int k = 0 ;
    for(int i=0;i<Loci::MPI_processes;++i) {
      int limit ;
      if(i == Loci::MPI_processes-1)
        limit = global_size ;
      else
        limit = recv_displs[i+1] ;
      for(;k<limit;++k)
        ret[i] += global_eset[k] ;
    }
    delete[] recv_displs ;
    delete[] global_eset ;

    return ret ;
  }

  void
  gather_all_entitySet_opt(const entitySet& eset,
                           vector<entitySet>& out) {
    bool pack_interval = false ;
    int elem_size = eset.size() ;
    int interval_size = eset.num_intervals() ;
    pack_interval = (2*interval_size) < elem_size ;

    int local_size = 0 ;
    if(eset != EMPTY) {
      if(pack_interval)
        local_size = 2*interval_size + 1 ;
      else
        local_size = elem_size + 1 ;
    }
    
    // compute receive counts from all processors
    vector<int> recv_counts(Loci::MPI_processes,0) ;
    MPI_Allgather(&local_size, 1, MPI_INT,
                  &recv_counts[0], 1, MPI_INT, MPI_COMM_WORLD) ;
    // then compute receive displacement
    vector<int> recv_displs(Loci::MPI_processes,0) ;
    recv_displs[0] = 0 ;
    for(int i=1;i<Loci::MPI_processes;++i)
      recv_displs[i] = recv_displs[i-1] + recv_counts[i-1] ;

    int global_size = recv_displs[Loci::MPI_processes-1] +
      recv_counts[Loci::MPI_processes-1] ;

    // pack the local eset into an array
    vector<int> local_eset(local_size) ;

    if(local_size > 0) {
      if(pack_interval) {
        local_eset[0] = 1 ;
        int count = 1 ;
        for(size_t i=0;i<eset.num_intervals();++i) {
          local_eset[count++] = eset[i].first ;
          local_eset[count++] = eset[i].second ;
        }
      } else {
        local_eset[0] = 0 ;
        int count = 1 ;
        for(entitySet::const_iterator ei=eset.begin();
            ei!=eset.end();++ei,++count)
          local_eset[count] = *ei ;
      }
    }
    
    // allocate the entire array for all data from all processors
    vector<int> global_eset(global_size) ;
    // communicate to obtain all esets from every processors
    MPI_Allgatherv(&local_eset[0], local_size, MPI_INT,
                   &global_eset[0], &recv_counts[0], &recv_displs[0],
                   MPI_INT, MPI_COMM_WORLD) ;
    vector<int>().swap(local_eset) ;
    // unpack the raw buffer into a vector<entitySet>
    if(out.size() != (size_t)Loci::MPI_processes)
      out.resize(Loci::MPI_processes) ;

    for(int i=0;i<Loci::MPI_processes;++i) {
      int b = recv_displs[i] ;
      int e = b + recv_counts[i] ;
      if(b==e) {
        out[i] = EMPTY ;
        continue ;
      }
      
      entitySet& eset = out[i] ;
      int flag = global_eset[b] ;
      ++b ;
      if(flag == 1) {
        // unpack intervals
        for(;b<e;b+=2) {
          int f = global_eset[b] ;
          int s = global_eset[b+1] ;
          eset += interval(f,s) ;
        }
      } else if(flag ==0) {
        // unpack elements
        for(;b<e;++b)
          eset += global_eset[b] ;
      } else {
        // impossible unless something is wrong
        cerr << "gather_all_entitySet_opt failed on process "
             << Loci::MPI_rank << ", with unpack flag = " << flag << endl ;
        Loci::Abort() ;
      }
    }

  }
  
  // this function takes an entitySet (in the fact_db's initial
  // global numbering) and a Map, it returns the image of the Map
  // based on the passed in entitySet as the domain. The returned
  // entitySet is also in fact_db's initial global numbering.
  // Note: the passed in entitySet may not reside on a single
  // process, hence this is a distributed Map image function.
  // "facts" is the pointer to the current fact_db in use.
  entitySet
  distributed_map_image(const entitySet& domain,
                        Loci::MapRepP m, Loci::fact_db* facts) {
    // first get fact_db's init partition in global numbering
    const vector<entitySet>& ptn = facts->get_init_ptn() ;
    // then get the distribute_info structure
    Loci::fact_db::distribute_infoP df = facts->get_distribute_info() ;  
    // figure out domain distributions
    vector<entitySet> domain_dist(Loci::MPI_processes) ;
    for(int i=0;i<Loci::MPI_processes;++i) {
      domain_dist[i] = domain & ptn[i] ;
    }
    // we will do a transpose on domain_dist to get all domains
    // entitySet that are on each local process
    vector<entitySet> domain_local =
      transpose_vector_entitySet(domain_dist) ;
    // we then proceed to get the images of the map m based
    // on the domain_local vector. but first we'll need to
    // convert the domain_local to the local numbering
    for(int i=0;i<Loci::MPI_processes;++i)
      domain_local[i] = remap_entitySet(domain_local[i], df->g2l) ;
    // then we can get the images
    vector<entitySet> images(Loci::MPI_processes) ;
    for(int i=0;i<Loci::MPI_processes;++i) {
      images[i] = m->image(domain_local[i]) ;
      // and we need to convert it back to global numbering
      images[i] = remap_entitySet(images[i], df->l2g) ;
    }
    // finally we transpose images again to return the image
    // values to corresponding process
    images = transpose_vector_entitySet(images) ;
    // merge the vector and return results
    entitySet im ;
    for(int i=0;i<Loci::MPI_processes;++i)
      im += images[i] ;

    return im ;
  }

#ifdef TO_BE_REMOVED
  // we need to remove functions in this section in the future
  // since they are obseleted.

  // this function expands a Map for the passed in domain (in
  // global numbering), however it returns a new Map that has the
  // expanded domain only, and the expanded domain is also renumbered
  // based on the passed in renumber map
  Loci::storeRepP
  expand_Map(const Map& m, const entitySet& edom,
             const dMap& remap, Loci::fact_db* facts) {
    // if edom is globally empty, we just return
    if(Loci::GLOBAL_AND(edom==EMPTY)) {
      Map nm ;
      return nm.Rep() ;
    }
    // first get fact_db's init partition in global numbering
    const vector<entitySet>& ptn = facts->get_init_ptn() ;
    // then get the distribute_info structure
    Loci::fact_db::distribute_infoP df = facts->get_distribute_info() ;  
    // figure out domain distributions
    vector<entitySet> edom_dist(Loci::MPI_processes) ;
    for(int i=0;i<Loci::MPI_processes;++i) {
      edom_dist[i] = edom & ptn[i] ;
    }
    // transpose "edom_dist"
    vector<entitySet> edom_dist_t = transpose_vector_entitySet(edom_dist) ;
    // convert to local numbering
    for(int i=0;i<Loci::MPI_processes;++i)
      edom_dist_t[i] = remap_entitySet(edom_dist_t[i], df->g2l) ;

    // prepare buffer to send/recv actual Map contents
    int* send_counts = new int[Loci::MPI_processes] ;
    for(int i=0;i<Loci::MPI_processes;++i)
      send_counts[i] = 2 * edom_dist_t[i].size() ;

    int* send_displs = new int[Loci::MPI_processes] ;
    send_displs[0] = 0 ;
    for(int i=1;i<Loci::MPI_processes;++i)
      send_displs[i] = send_displs[i-1] + send_counts[i-1] ;

    // then communicate this to get the recv info.
    int* recv_counts = new int[Loci::MPI_processes] ;
    MPI_Alltoall(send_counts, 1, MPI_INT,
                 recv_counts, 1, MPI_INT, MPI_COMM_WORLD) ;

    int* recv_displs = new int[Loci::MPI_processes] ;
    recv_displs[0] = 0 ;
    for(int i=1;i<Loci::MPI_processes;++i)
      recv_displs[i] = recv_displs[i-1] + recv_counts[i-1] ;

    // all the size info have been gathered, we will need
    // to pack the contents into a raw buffer and communicate
    int buf_size = 0 ;
    for(int i=0;i<Loci::MPI_processes;++i)
      buf_size += send_counts[i] ;

    int* send_buf = new int[buf_size] ;
    int buf_idx = 0 ;
    for(int i=0;i<Loci::MPI_processes;++i) {
      const entitySet& es = edom_dist_t[i] ;
      for(entitySet::const_iterator ei=es.begin();
          ei!=es.end();++ei) {
        // first store the domain entity (in global numbering)
        send_buf[buf_idx++] = df->l2g[*ei] ;
        // then store the content entity (in global numbering)
        send_buf[buf_idx++] = df->l2g[m[*ei]] ;
      }
    }
  
    // allocate receive buffer
    int recv_size = 0 ;
    for(int i=0;i<Loci::MPI_processes;++i)
      recv_size += recv_counts[i] ;

    int* recv_buf = new int[recv_size] ;
    // communicate
    MPI_Alltoallv(send_buf, send_counts,
                  send_displs, MPI_INT,
                  recv_buf, recv_counts,
                  recv_displs, MPI_INT, MPI_COMM_WORLD) ;

    delete[] send_counts ;
    delete[] send_displs ;
    delete[] recv_counts ;
    delete[] send_buf ;
    // extract the buffer into a new Map
    Map nm ;
    entitySet nm_dom = remap_entitySet(edom, remap) ;
    nm.allocate(nm_dom) ;

    int k=0 ;
    for(int i=0;i<Loci::MPI_processes;++i) {
      int limit ;
      if(i == Loci::MPI_processes-1)
        limit = recv_size ;
      else
        limit = recv_displs[i+1] ;

      for(;k<limit;k+=2) {
        int first = remap[recv_buf[k]] ;
        int second = remap[recv_buf[k+1]] ;
        nm[first] = second ;
      }
    }

    delete[] recv_displs ;
    delete[] recv_buf ;

    return nm.Rep() ;
  }

  // This is a special expand_Map function, everything is
  // the same as the expand_Map function above, except for
  // that in the newly expanded Map, only its domain is
  // renumbered, the image of the Map is left unchanged (not remapped).
  Loci::storeRepP
  expand_Map2(const Map& m, const entitySet& edom,
              const dMap& remap, Loci::fact_db* facts) {
    // if edom is globally empty, we just return
    if(Loci::GLOBAL_AND(edom==EMPTY)) {
      Map nm ;
      return nm.Rep() ;
    }
    // first get fact_db's init partition in global numbering
    const vector<entitySet>& ptn = facts->get_init_ptn() ;
    // then get the distribute_info structure
    Loci::fact_db::distribute_infoP df = facts->get_distribute_info() ;  
    // figure out domain distributions
    vector<entitySet> edom_dist(Loci::MPI_processes) ;
    for(int i=0;i<Loci::MPI_processes;++i) {
      edom_dist[i] = edom & ptn[i] ;
    }
    // transpose "edom_dist"
    vector<entitySet> edom_dist_t = transpose_vector_entitySet(edom_dist) ;
    // convert to local numbering
    for(int i=0;i<Loci::MPI_processes;++i)
      edom_dist_t[i] = remap_entitySet(edom_dist_t[i], df->g2l) ;

    // prepare buffer to send/recv actual Map contents
    int* send_counts = new int[Loci::MPI_processes] ;
    for(int i=0;i<Loci::MPI_processes;++i)
      send_counts[i] = 2 * edom_dist_t[i].size() ;

    int* send_displs = new int[Loci::MPI_processes] ;
    send_displs[0] = 0 ;
    for(int i=1;i<Loci::MPI_processes;++i)
      send_displs[i] = send_displs[i-1] + send_counts[i-1] ;

    // then communicate this to get the recv info.
    int* recv_counts = new int[Loci::MPI_processes] ;
    MPI_Alltoall(send_counts, 1, MPI_INT,
                 recv_counts, 1, MPI_INT, MPI_COMM_WORLD) ;

    int* recv_displs = new int[Loci::MPI_processes] ;
    recv_displs[0] = 0 ;
    for(int i=1;i<Loci::MPI_processes;++i)
      recv_displs[i] = recv_displs[i-1] + recv_counts[i-1] ;

    // all the size info have been gathered, we will need
    // to pack the contents into a raw buffer and communicate
    int buf_size = 0 ;
    for(int i=0;i<Loci::MPI_processes;++i)
      buf_size += send_counts[i] ;

    int* send_buf = new int[buf_size] ;
    int buf_idx = 0 ;
    for(int i=0;i<Loci::MPI_processes;++i) {
      const entitySet& es = edom_dist_t[i] ;
      for(entitySet::const_iterator ei=es.begin();
          ei!=es.end();++ei) {
        // first store the domain entity (in global numbering)
        send_buf[buf_idx++] = df->l2g[*ei] ;
        // then store the content entity (in global numbering)
        send_buf[buf_idx++] = df->l2g[m[*ei]] ;
      }
    }
  
    // allocate receive buffer
    int recv_size = 0 ;
    for(int i=0;i<Loci::MPI_processes;++i)
      recv_size += recv_counts[i] ;

    int* recv_buf = new int[recv_size] ;
    // communicate
    MPI_Alltoallv(send_buf, send_counts,
                  send_displs, MPI_INT,
                  recv_buf, recv_counts,
                  recv_displs, MPI_INT, MPI_COMM_WORLD) ;

    delete[] send_counts ;
    delete[] send_displs ;
    delete[] recv_counts ;
    delete[] send_buf ;
    // extract the buffer into a new Map
    Map nm ;
    entitySet nm_dom = remap_entitySet(edom, remap) ;
    nm.allocate(nm_dom) ;

    int k=0 ;
    for(int i=0;i<Loci::MPI_processes;++i) {
      int limit ;
      if(i == Loci::MPI_processes-1)
        limit = recv_size ;
      else
        limit = recv_displs[i+1] ;

      for(;k<limit;k+=2) {
        int first = remap[recv_buf[k]] ;
        // do not renumber the image!
        int second = recv_buf[k+1] ;
        nm[first] = second ;
      }
    }

    delete[] recv_displs ;
    delete[] recv_buf ;

    return nm.Rep() ;
  }

  // this function expands a multiMap for the passed in domain (in
  // global numbering), however it returns a new multiMap that has the
  // expanded domain only, and the expanded domain is also renumbered
  // based on the passed in renumber map
  Loci::storeRepP
  expand_multiMap(const multiMap& m, const entitySet& edom,
                  const dMap& remap, Loci::fact_db* facts) {
    // if edom is globally empty, then just return
    if(Loci::GLOBAL_AND(edom==EMPTY)) {
      multiMap nmm ;
      return nmm.Rep() ;
    }
    // first get fact_db's init partition in global numbering
    const vector<entitySet>& ptn = facts->get_init_ptn() ;
    // then get the distribute_info structure
    Loci::fact_db::distribute_infoP df = facts->get_distribute_info() ;  
    // figure out domain distributions
    vector<entitySet> edom_dist(Loci::MPI_processes) ;
    for(int i=0;i<Loci::MPI_processes;++i) {
      edom_dist[i] = edom & ptn[i] ;
    }
    // transpose "edom_dist"
    vector<entitySet> edom_dist_t = transpose_vector_entitySet(edom_dist) ;
    // convert to local numbering
    for(int i=0;i<Loci::MPI_processes;++i)
      edom_dist_t[i] = remap_entitySet(edom_dist_t[i], df->g2l) ;

    // prepare buffer to send/recv actual Map contents
    int* send_counts = new int[Loci::MPI_processes] ;
    for(int i=0;i<Loci::MPI_processes;++i) {
      send_counts[i] = 2 * edom_dist_t[i].size() ;
      for(entitySet::const_iterator ei=edom_dist_t[i].begin();
          ei!=edom_dist_t[i].end();++ei)
        send_counts[i] += m.num_elems(*ei) ;
    }

    int* send_displs = new int[Loci::MPI_processes] ;
    send_displs[0] = 0 ;
    for(int i=1;i<Loci::MPI_processes;++i)
      send_displs[i] = send_displs[i-1] + send_counts[i-1] ;

    // then communicate this to get the recv info.
    int* recv_counts = new int[Loci::MPI_processes] ;
    MPI_Alltoall(send_counts, 1, MPI_INT,
                 recv_counts, 1, MPI_INT, MPI_COMM_WORLD) ;

    int* recv_displs = new int[Loci::MPI_processes] ;
    recv_displs[0] = 0 ;
    for(int i=1;i<Loci::MPI_processes;++i)
      recv_displs[i] = recv_displs[i-1] + recv_counts[i-1] ;

    // all the size info have been gathered, we will need
    // to pack the contents into a raw buffer and communicate
    int buf_size = 0 ;
    for(int i=0;i<Loci::MPI_processes;++i)
      buf_size += send_counts[i] ;

    int* send_buf = new int[buf_size] ;
    int buf_idx = 0 ;
    for(int i=0;i<Loci::MPI_processes;++i) {
      const entitySet& es = edom_dist_t[i] ;
      for(entitySet::const_iterator ei=es.begin();
          ei!=es.end();++ei) {
        // first store the domain entity (in global numbering)
        send_buf[buf_idx++] = df->l2g[*ei] ;
        int size = m.num_elems(*ei) ;
        // then store the size of the elements
        send_buf[buf_idx++] = size ;
        // finally store the elements (in global numbering)
        for(int k=0;k<size;++k)
          send_buf[buf_idx++] = df->l2g[m[*ei][k]] ;
      }
    }
  
    // allocate receive buffer
    int recv_size = 0 ;
    for(int i=0;i<Loci::MPI_processes;++i)
      recv_size += recv_counts[i] ;

    int* recv_buf = new int[recv_size] ;
    // communicate
    MPI_Alltoallv(send_buf, send_counts,
                  send_displs, MPI_INT,
                  recv_buf, recv_counts,
                  recv_displs, MPI_INT, MPI_COMM_WORLD) ;

    delete[] send_counts ;
    delete[] send_displs ;
    delete[] recv_counts ;
    delete[] send_buf ;
    // extract the buffer into a new Map
    multiMap nm ;
    entitySet nm_dom = remap_entitySet(edom, remap) ;

    store<int> alloc ;
    alloc.allocate(nm_dom) ;
  
    // figure out the allocate size
    int k=0;
    for(int i=0;i<Loci::MPI_processes;++i) {
      int limit ;
      if(i == Loci::MPI_processes-1)
        limit = recv_size ;
      else
        limit = recv_displs[i+1] ;

      for(;k<limit;) {
        int first = remap[recv_buf[k]] ;
        int size = recv_buf[k+1] ;
        alloc[first] = size ;
        k = k + 2 + size ;
      }
    }

    nm.allocate(alloc) ;

    // fill in the actual contents
    k=0 ;
    for(int i=0;i<Loci::MPI_processes;++i) {
      int limit ;
      if(i == Loci::MPI_processes-1)
        limit = recv_size ;
      else
        limit = recv_displs[i+1] ;

      for(;k<limit;) {
        int first = remap[recv_buf[k]] ;
        int size = recv_buf[k+1] ;
        k += 2 ;
        for(int j=0;j<size;++j)
          nm[first][j] = remap[recv_buf[k+j]] ;
        k += size ;
      }
    }

    delete[] recv_displs ;
    delete[] recv_buf ;

    return nm.Rep() ;
  }
#endif

  /////////////////////////////////////////////////////////////
  // user function for mpi_2int reduction (SUM, MAX)
  void
  SUM_MAX_2INT(mpi_2int *invec, mpi_2int *inoutvec,
               int *len, MPI_Datatype *dtype) {
    int i;
    for ( i=0; i<*len; i++ ) {
      inoutvec[i].first += invec[i].first ;
      if(invec[i].second > inoutvec[i].second)
        inoutvec[i].second = invec[i].second ;
    }
  }

  // functions that computes median value for a group of processes
  int
  mpi_median(int value, MPI_Comm comm) {
    // get total number of processes
    int p = 0 ;
    MPI_Comm_size(comm,&p) ;

    vector<int> v(p) ;
    // gather values onto all processes
    MPI_Allgather(&value, 1, MPI_INT,
                  &v[0], 1, MPI_INT, comm) ;

    if(p%2 != 0) {
      // the middle of the element
      int m = p/2 ;
      std::nth_element(v.begin(), v.begin()+m, v.end()) ;
      return v[m] ;
    } else {
      // mean of the two middle elements
      int m = p/2 ;
      std::nth_element(v.begin(), v.begin()+m, v.end()) ;
      int n1 = v[m] ;
      std::nth_element(v.begin(), v.begin()+m-1, v.begin()+m) ;
      int n2 = v[m-1] ;
      return (n1+n2)/2 ;
    }
  }

  long
  mpi_median(long value, MPI_Comm comm) {
    // get total number of processes
    int p = 0 ;
    MPI_Comm_size(comm,&p) ;

    vector<long> v(p) ;
    // gather values onto all processes
    MPI_Allgather(&value, 1, MPI_LONG,
                  &v[0], 1, MPI_LONG, comm) ;

    if(p%2 != 0) {
      // the middle of the element
      int m = p/2 ;
      std::nth_element(v.begin(), v.begin()+m, v.end()) ;
      return v[m] ;
    } else {
      // mean of the two middle elements
      int m = p/2 ;
      std::nth_element(v.begin(), v.begin()+m, v.end()) ;
      long n1 = v[m] ;
      std::nth_element(v.begin(), v.begin()+m-1, v.begin()+m) ;
      long n2 = v[m-1] ;
      return (n1+n2)/2 ;
    }
  }
  
  double
  mpi_median(double value, MPI_Comm comm) {
    // get total number of processes
    int p = 0 ;
    MPI_Comm_size(comm,&p) ;

    vector<double> v(p) ;
    // gather values onto all processes
    MPI_Allgather(&value, 1, MPI_DOUBLE,
                  &v[0], 1, MPI_DOUBLE, comm) ;

    if(p%2 != 0) {
      // the middle of the element
      int m = p/2 ;
      std::nth_element(v.begin(), v.begin()+m, v.end()) ;
      return v[m] ;
    } else {
      // mean of the two middle elements
      int m = p/2 ;
      std::nth_element(v.begin(), v.begin()+m, v.end()) ;
      double n1 = v[m] ;
      std::nth_element(v.begin(), v.begin()+m-1, v.begin()+m) ;
      double n2 = v[m-1] ;
      return (n1+n2)/2 ;
    }
  }
  
  void
  generate_process_map(const std::vector<entitySet>& ptn,
                       const std::vector<int>& items,
                       std::vector<int>& pmap) {
    std::vector<int> distinct_items ;
    std::back_insert_iterator<std::vector<int> > di(distinct_items) ;
    std::unique_copy(items.begin(), items.end(), di) ;
    entitySet all_items =
      Loci::create_intervalSet(distinct_items.begin(), distinct_items.end()) ;
    std::vector<entitySet> dist(ptn.size()) ;
    entitySet tmp ;
    for(size_t i=0;i!=ptn.size();++i) {
      dist[i] = all_items & ptn[i] ;
      tmp += dist[i] ;
    }
    if(tmp != all_items)
      throw Loci::StringError
        ("generate_process_map failed! there are items outside of ptn") ;
    dMap proc_map ;
    for(size_t p=0;p!=dist.size();++p) {
      for(entitySet::const_iterator ei=dist[p].begin();
          ei!=dist[p].end();++ei)
        proc_map[*ei] = p ;
    }
    pmap.resize(items.size(), 0) ;
    for(size_t i=0;i<pmap.size();++i)
      pmap[i] = proc_map[items[i]] ;
  }

  ///////////////////////////////////////////////////
  // these are the codes for hilbert curve partition
  ///////////////////////////////////////////////////
  // mask for 3D
  const unsigned int g_mask[] = {4,2,1} ;
  
  HilbertCode hilbert_encode(const IntCoord3& p) {
    const int DIM = 3 ;
    const int WORDBITS = 32 ;
    const int NUMBITS = 32 ;
    unsigned int mask = (unsigned long)1 << (WORDBITS - 1) ;
    unsigned int element, temp1, temp2, A, W = 0, S, tS, T, tT, J, P = 0, xJ;
    
    HilbertCode	h;
    h[0] = 0;
    h[1] = 0;
    h[2] = 0;
    
    int	i = NUMBITS * DIM - DIM, j;
    
    for (j = A = 0; j < DIM; j++)
      if (p[j] & mask)
        A |= g_mask[j];
    
    S = tS = A;
    
    P |= S & g_mask[0];
    for (j = 1; j < DIM; j++)
      if( (S & g_mask[j]) ^ ((P >> 1) & g_mask[j]))
        P |= g_mask[j];
    
    /* add in DIM bits to hcode */
    element = i / WORDBITS;
    if (i % WORDBITS > WORDBITS - DIM) {
        h[element] |= P << i % WORDBITS;
        h[element + 1] |= P >> (WORDBITS - i % WORDBITS) ;
    } else
      h[element] |= P << (i - element * WORDBITS) ;

    J = DIM;
    for (j = 1; j < DIM; j++)
      if ((P >> j & 1) == (P & 1))
        continue;
      else
        break;
    if (j != DIM)
      J -= j;
    xJ = J - 1;
    
    if (P < 3)
      T = 0;
    else
      if (P % 2)
        T = (P - 1) ^ (P - 1) / 2;
      else
        T = (P - 2) ^ (P - 2) / 2;
    tT = T;
    
    for (i -= DIM, mask >>= 1; i >=0; i -= DIM, mask >>= 1) {
      for (j = A = 0; j < DIM; j++)
        if (p[j] & mask)
          A |= g_mask[j];
      
      W ^= tT;
      tS = A ^ W;
      if (xJ % DIM != 0) {
        temp1 = tS << xJ % DIM;
        temp2 = tS >> (DIM - xJ % DIM) ;
        S = temp1 | temp2;
        S &= ((unsigned int)1 << DIM) - 1;
      } else
        S = tS;

      P = S & g_mask[0];
      for (j = 1; j < DIM; j++)
        if( (S & g_mask[j]) ^ ((P >> 1) & g_mask[j]))
          P |= g_mask[j];
      
      /* add in DIM bits to hcode */
      element = i / WORDBITS;
      if (i % WORDBITS > WORDBITS - DIM) {
        h[element] |= P << (i % WORDBITS) ;
        h[element + 1] |= P >> (WORDBITS - i % WORDBITS) ;
      } else
        h[element] |= P << (i - element * WORDBITS);

      if (i > 0) {
        if (P < 3)
          T = 0;
        else
          if (P % 2)
            T = (P - 1) ^ (P - 1) / 2;
          else
            T = (P - 2) ^ (P - 2) / 2;
        
        if (xJ % DIM != 0) {
          temp1 = T >> (xJ % DIM) ;
          temp2 = T << (DIM - xJ % DIM) ;
          tT = temp1 | temp2;
          tT &= ((unsigned int)1 << DIM) - 1;
        } else
          tT = T;
        
        J = DIM;
        for (j = 1; j < DIM; j++)
          if ((P >> j & 1) == (P & 1))
            continue;
          else
            break;
        if (j != DIM)
          J -= j;
        
        xJ += J - 1;

      }
    }
    return h;
  }
  

} // end of namespace lagrangianP
