//#############################################################################
//#
//# Copyright 2014-2017, Mississippi State University
//#
//# The GridMover is free software: you can redistribute it and/or modify
//# it under the terms of the Lesser GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The GridMover software is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# Lesser GNU General Public License for more details.
//#
//# You should have received a copy of the Lesser GNU General Public License
//# along with the GridMover software. If not, see <http://www.gnu.org/licenses>
//#
//#############################################################################
#include "iws.h"

using std::vector;
using std::cout ;
using std::endl ;
using std::cerr ;



namespace gridMotion {
  
  // Initialize the schedule for first time (note, the size of data on each
  // processor is assumed to remain fixed throught execution)
  void IWSSchedule::initialize(MPI_Comm C, int local_size)
  {
    using std::max ;
    comm = C ;
    
    MPI_Comm_rank(comm, &MPI_rank);
    MPI_Comm_size(comm, &MPI_processes);
    
    sendChunks.clear() ;
    recvChunks.clear() ;
    //Finding the size of the iterate space.
    int lsz = local_size ;

    int gsz = 0;
    MPI_Allreduce(&lsz, &gsz, 1, MPI_INT, MPI_SUM, comm);
    //Chunk size
    int p;
    int MPI_rank;
    MPI_Comm_rank(comm, &MPI_rank);  
    MPI_Comm_size(comm, &p);
    int CHUNKING_FACTOR = 256;
    int csz = max(gsz/(p*CHUNKING_FACTOR), 10);
    int nc = max(lsz/csz-1, 0);
    int first_chunk = lsz - csz*nc;
      
    int start = 0;
    int last = local_size ;
  
    chunkInfo tmp;
    tmp.start = 0;
    tmp.end = start + first_chunk - 1;
    tmp.chunkTime[0] = 0;
    tmp.chunkTime[1] = 0;
    chunkData.push_back(tmp);


    for (start = first_chunk; start < last; start += csz) {
      tmp.start = start;
      tmp.end = start+csz-1;
      chunkData.push_back(tmp);
    }
  
    int nchunks = chunkData.size();
    for(int i = 0; i < nchunks; ++i)
      selfChunks.push_back(i);

    numBalances = 0;
    numRemoteChunks = 0;
  }
  
  // Generate a work migration schedule based on previously measured chunk times
  void IWSSchedule::generateSchedule() {
    using std::vector ;
    using std::abs ;
    using std::max ;
    using std::cout ;
    using std::cerr ;
    using std::endl ;
    using std::pair ;
    const float eff_tol = 0.01 ;
    // Compute balanced schedule, first compute total time
    float time = 0 ;
    for(size_t i=0;i<chunkData.size();++i)
      time += chunkData[i].chunkTime[0]+chunkData[i].chunkTime[1] ;
    float total_time = 0 ;
#ifdef VERBOSE
    cout << "reducing " << time << endl ;
#endif
    MPI_Allreduce(&time,&total_time,1,MPI_FLOAT,MPI_SUM,comm) ;
    float mean_time = total_time/float(MPI_processes) ;
    //    cout << mean_time << " (mean time) = " << total_time << " (total time) / " << MPI_processes << " (MPI_proc)." << endl;
    //    cout << "My time: " << time << endl;
    float diff_time = time-mean_time ;
    //    cout << "Diff time: " << diff_time << endl;
    if(fabs(diff_time) < mean_time*eff_tol) // If close to average time, zero
      diff_time = 0 ;
    
    // Compute chunks that will be sent to other processors
    vector<int> send_list ;
    if(diff_time > 0) { // If this processor is taking too much time, 
      // allocate chunks that get it to the mean_time ;
      vector<pair<float,int> > chunk_times(chunkData.size()-1) ;
      for(size_t i=1;i<chunkData.size();++i) {
        chunk_times[i-1].first = 
          chunkData[i].chunkTime[0]+chunkData[i].chunkTime[1] ;
        chunk_times[i-1].second = i ;
      }
      std::sort(chunk_times.begin(),chunk_times.end()) ;
      
      
#define SELECT_GIVE
#ifdef SELECT_GIVE
      // select chunks to give away
      float time_give = 0 ;
      for(int kk = chunk_times.size()-1;kk>=0;--kk) {
        if(time_give > diff_time)
          break ;
        if(time_give+chunk_times[kk].first < diff_time+mean_time*eff_tol) {
          send_list.push_back(chunk_times[kk].second) ;
          time_give += chunk_times[kk].first ;
        }
      }
#endif
#ifdef SELECT_KEEP
      // This code for this processor to keep what it can to meet mean time
      float time_keep = chunkData[0].chunkTime[0]+chunkData[0].chunkTime[1] ;
      // keep the smallest timed chunks up to keep_percent of mean to reduce
      // number of chunks communicated
      float keep_percent = 0.01 ;
      int stop = 0 ;
      while(time_keep+chunk_times[stop].first < keep_percent*mean_time &&
            stop < int(chunk_times.size())) {
        time_keep += chunk_times[stop].first ;
        stop++ ;
      }
      
      int kk ;
      for(kk = chunk_times.size()-1;kk>=stop;--kk) {
        if(time_keep > mean_time)
          break ;
        if(time_keep+chunk_times[kk].first < mean_time*(1.+eff_tol)) {
          time_keep += chunk_times[kk].first ;
        } else
          send_list.push_back(chunk_times[kk].second) ;
      }
      for(;kk>=0;--kk)
        send_list.push_back(chunk_times[kk].second) ;
#endif
    }
    if(send_list.size() != 0) {
      diff_time = 0 ;
      for(size_t i=0;i<send_list.size();++i) {
        diff_time += (chunkData[send_list[i]].chunkTime[0]+
                      chunkData[send_list[i]].chunkTime[1]) ;
      }
    }
    vector<float> time_xfers(MPI_processes) ;
    MPI_Allgather(&diff_time, 1, MPI_FLOAT, &time_xfers[0],1,MPI_FLOAT,
                  comm) ;
#ifdef VERBOSE
    if(MPI_rank == 0) {
      cout << "mean_time = " << mean_time << endl ;
      cout << "time xfers =";
      for(int i=0;i<MPI_processes;++i)
        cout << " " << time_xfers[i] ;
      cout << endl ;
    }
#endif
    vector<pair<float,int> > send_chunks ;
    vector<pair<float,int> > recv_chunks ;
    for(int i=0;i<MPI_processes;++i) {
      if(time_xfers[i] > 0.0) 
        send_chunks.push_back(pair<float,int>(time_xfers[i],i)) ;
      if(time_xfers[i] < 0.0) 
        recv_chunks.push_back(pair<float,int>(-time_xfers[i],i)) ;
    }
    sort(send_chunks.begin(),send_chunks.end()) ;
    sort(recv_chunks.begin(),recv_chunks.end()) ;
    
    // Compute chunk sendto schedule
    vector<pair<int,vector<int> > > sendto ;
    
    for(int i=recv_chunks.size()-1;i>=0;--i) {
      double ptime = recv_chunks[i].first ;
      for(int j=send_chunks.size()-1;j>=0;--j) {
        if(send_chunks[j].second >=0) {
          if(ptime - send_chunks[j].first > -mean_time*eff_tol) {
            ptime -= send_chunks[j].first ;
            //assign all chunks from this processor
            if(send_chunks[j].second == MPI_rank) {
#ifdef VERBOSE
              cout << "adding " << send_list.size() << "chunks to sendto"
                   << endl ;
#endif
              if(send_list.size() > 0) 
                sendto.push_back(pair<int,vector<int> >(recv_chunks[i].second,
                                                        send_list)) ;
              send_list.clear() ;
            }
            time_xfers[send_chunks[j].second] = 0 ;
            send_chunks[j].second = -1 ;
            if(ptime < 0)
              break ;
          }
        }
      }
      ptime = -max(ptime,0.0) ;
      if(fabs(ptime) < mean_time*eff_tol)
        ptime = 0 ;
      time_xfers[recv_chunks[i].second] = ptime ;
    }
    
#ifdef VERBOSE
    
    if(MPI_rank == 0) {
      cout << "time xfers2 =";
      for(int i=0;i<MPI_processes;++i)
        cout << " " << time_xfers[i] ;
      cout << endl ;
    }
    
#endif
    
    bool rebalance = false ;
    for(int i=0;i<MPI_processes;++i)
      if(time_xfers[i] > 0)
        rebalance = true ;
    bool sources = false ;
    for(int i=0;i<MPI_processes;++i)
      if(time_xfers[i] < 0)
        sources = true ;
    if(rebalance && sources) { // we have residual work to distribute
      int nchunks = send_list.size() ;
      vector<int> chunk_groups(MPI_processes) ;
      MPI_Allgather(&nchunks, 1, MPI_INT, &chunk_groups[0],1,MPI_INT,
                    comm) ;
      vector<float> chunk_times(nchunks) ;
      for(int i=0;i<nchunks;++i) {
        int ch = send_list[i] ;
        chunk_times[i] = (chunkData[ch].chunkTime[0]+
                          chunkData[ch].chunkTime[1]) ;
      }
      vector<int> chunk_displ(MPI_processes) ;
      chunk_displ[0] = 0 ;
      for(int i=1;i<MPI_processes;++i) 
        chunk_displ[i] = chunk_displ[i-1]+chunk_groups[i-1] ;
      
      int ntchunks = (chunk_displ[MPI_processes-1]+
                      chunk_groups[MPI_processes-1]) ;
      
      vector<float> chunk_time_gather(ntchunks) ;
      
      MPI_Allgatherv(&chunk_times[0],nchunks,MPI_FLOAT,
                     &chunk_time_gather[0],
                     &chunk_groups[0],
                     &chunk_displ[0],
                     MPI_FLOAT,
                     comm) ;
      vector<pair<float,pair<int,int> > > chunkQueue(ntchunks) ;
      int cnk = 0 ;
      for(int i=0;i<MPI_processes;++i) {
        for(int j=0;j<chunk_groups[i];++j) {
          pair<int,int> chunk_info(i,j) ;
          float chunk_time = chunk_time_gather[chunk_displ[i]+j] ;
          chunkQueue[cnk] = pair<float,pair<int,int> >(chunk_time,
                                                       chunk_info) ;
          cnk++ ;
        }
      }
      cnk = 0 ;
      for(int i=0;i<MPI_processes;++i) {
        if(chunk_groups[i] > 1) {
          std::sort(&chunkQueue[cnk],&chunkQueue[cnk+chunk_groups[i]]) ;
        }
        cnk += chunk_groups[i] ;
      }
      
      //	sort(chunkQueue.begin(),chunkQueue.end()) ;
      
      int chunkQueueStart = chunkQueue.size()-1 ;
      for(int i=MPI_processes-1;i>=0;--i) 
        if(time_xfers[i] < 0) {
          vector<int> sendto_p ;
          float time_x = -time_xfers[i] ;
          for(int j=chunkQueueStart;j >=0;--j)
            if(time_x > 0 && (chunkQueue[j].first > 0) &&
               (time_x - chunkQueue[j].first > -mean_time*eff_tol)) {
              // assign chunk 
              time_x -= chunkQueue[j].first ;
              const int cp = chunkQueue[j].second.first ;
              if(cp == MPI_rank)
                sendto_p.push_back(send_list[chunkQueue[j].second.second]) ;
              time_xfers[cp] -= chunkQueue[j].first ;
              if(time_xfers[cp] < mean_time*eff_tol)
                time_xfers[cp] = 0 ;
              chunkQueue[j].first = -1.0 ; // remove from consideration
              if(time_x < 0)
                break ;
            }
          //skip deleted entries
          for(;chunkQueueStart>0;--chunkQueueStart)
            if(chunkQueue[chunkQueueStart].first > 0)
              break ;
          if(sendto_p.size() > 0) {
            sendto.push_back(pair<int,vector<int> >(i,sendto_p)) ;
          }
          time_x = -time_x ;
          if(fabs(time_x) < mean_time*eff_tol)
            time_x = 0 ;
          time_xfers[i] = time_x ;
        }
      
      
      if(MPI_rank == 0) {
        bool found = false ;
        for(int j=chunkQueueStart;j >=0;--j) 
          if((chunkQueue[j].first > 0)) {
            found = true ;
          }
        if(found)  {
          cout << "chunks remaining in queue:" << endl ;
          for(int j=chunkQueueStart;j >=0;--j) 
            if((chunkQueue[j].first > 0)) {
              cout << "chunk from p=" << chunkQueue[j].second.first
                   << ", time = " << chunkQueue[j].first << endl ;
            }
        }
      }
      
      
    }
#ifdef VERBOSE
        
    if(MPI_rank == 0) {
      cout << "time xfers3 =";
      for(int i=0;i<MPI_processes;++i)
        cout << " " << time_xfers[i] ;
      cout << endl ;
    }
    
#endif
    
    
    
    //------------------------------------------------------------------
    // convert sendto to execution and communication schedules
    //------------------------------------------------------------------
    
#ifdef VERBOSE
    
    cout << "sendto:" << endl ;
    for(size_t i=0;i<sendto.size();++i)
      cout << sendto[i].second.size() << "chunks to " <<
        sendto[i].first << " from " << MPI_rank << endl ;
    
#endif
    int numChunks = chunkData.size() ;
    // Convert sendto to schedule
    vector<int> chunk_list(numChunks,-1) ;
    for(size_t i=0;i<sendto.size();++i)
      for(size_t j=0;j<sendto[i].second.size();++j)
        chunk_list[sendto[i].second[j]] = sendto[i].first ;
    // Setup local execution list
    
    vector<int> local_list ;
    for(int i=0;i<numChunks;++i)
      if(chunk_list[i] == -1)
        local_list.push_back(i) ;
    selfChunks.swap(local_list) ;
    
    int local_chk_info[2],global_chk_info[2] ;
    local_chk_info[0] = selfChunks.size() ;
    local_chk_info[1] = numChunks ;
    MPI_Allreduce(&local_chk_info[0],&global_chk_info[0],2,MPI_INT,
                  MPI_SUM, comm) ;
#ifndef SILENT    
    if(MPI_rank == 0) {
      cout << "IWS Schedule: Communicating "
           << global_chk_info[1]-global_chk_info[0]
           << " chunks, "
           << 100.0*(1.-double(global_chk_info[0])/double(global_chk_info[1])) << "% of all chunks." << endl ;
    }
#endif
    
    vector<IWSSchedule::chunkCommInfo> local_send_list ;
    for(size_t i=0;i<sendto.size();++i) {
      IWSSchedule::chunkCommInfo tmp ;
      tmp.proc = sendto[i].first ;
      tmp.chunkList = sendto[i].second ;

      int sz = 0 ;
      for(int j=0;j<tmp.chunkList.size();++j)
        sz += chunkData[tmp.chunkList[j]].end -
          chunkData[tmp.chunkList[j]].start+1 ;
      
      tmp.send_size = sz ;
      tmp.recv_size = sz ;
      local_send_list.push_back(tmp);
    }
    sendChunks.swap(local_send_list) ;
    
    // now invert sendChunks
    vector<int> sendSizes(MPI_processes,0) ;
    for(size_t i=0;i<sendChunks.size();++i)
      sendSizes[sendChunks[i].proc] = sendChunks[i].chunkList.size() ;
    
    vector<int> recvSizes(MPI_processes,0) ;
    MPI_Alltoall(&sendSizes[0],1,MPI_INT,
                 &recvSizes[0],1,MPI_INT,
                 comm) ;
    numRemoteChunks = 0 ;
    
    vector<IWSSchedule::chunkCommInfo> local_recv_list ;
    for(int i=0;i<MPI_processes;++i) {
      numRemoteChunks += recvSizes[i] ;
      if(recvSizes[i]!=0) {
        IWSSchedule::chunkCommInfo tmp ;
        tmp.proc = i ;
        for(int k=0;k<recvSizes[i];++k)
          tmp.chunkList.push_back(k) ;
        local_recv_list.push_back(tmp) ;
      }
    }
    if(numRemoteChunks != 0 && sendChunks.size() !=0) {
      cerr << "logic error in iterative weighted static LB method" << endl ;
      //TODO: REPLACE LOCI:ABORT HERE.
      //Loci::Abort() ;
    }
    for(size_t i=0;i<local_recv_list.size();++i) {
      MPI_Status tStatus;
      MPI_Recv(&local_recv_list[i].recv_size,1,MPI_INT,
               local_recv_list[i].proc,TAG_INFO,
               comm,&tStatus) ;
    }
    for(size_t i=0;i<sendChunks.size();++i) {
      MPI_Send(&sendChunks[i].send_size,1,MPI_INT,
               sendChunks[i].proc,TAG_INFO,comm) ;
    }
    for(size_t i=0;i<local_recv_list.size();++i) {
      MPI_Status tStatus;
      MPI_Recv(&local_recv_list[i].send_size,1,MPI_INT,
               local_recv_list[i].proc,TAG_INFO,
               comm,&tStatus) ;
    }
    for(size_t i=0;i<sendChunks.size();++i) {
      MPI_Send(&sendChunks[i].recv_size,1,MPI_INT,
               sendChunks[i].proc,TAG_INFO,comm) ;
    }
    recvChunks.swap(local_recv_list) ;
#ifdef VERBOSE
    
    if(sendChunks.size() != 0) {
      cout << "sendChunks: " << endl ;
      for(size_t i=0;i<sendChunks.size();++i)
        cout << sendChunks[i].proc << ' ' << sendChunks[i].chunkList.size() << ' ' << sendChunks[i].send_size << ' ' << sendChunks[i].recv_size << endl ;
    }
    if(recvChunks.size() != 0) {
      cout << "recvChunks: " << endl ;
      for(size_t i=0;i<recvChunks.size();++i)
        cout << recvChunks[i].proc << ' ' << recvChunks[i].chunkList.size() << ' ' << recvChunks[i].send_size << ' ' << recvChunks[i].recv_size << endl ;
    }
    
#endif
    int nsendrecvl = recvChunks.size()+sendChunks.size() ;
    int nsendrecvg = 0 ;
    int nsendrecvs = 0 ;
    MPI_Allreduce(&nsendrecvl,&nsendrecvg,1,MPI_INT,MPI_MAX, comm) ;
    MPI_Allreduce(&nsendrecvl,&nsendrecvs,1,MPI_INT,MPI_SUM, comm) ;

#ifndef SILENT
    if(MPI_rank == 0) {
      cout << "IWS Schedule: Each processor communicating with an average of "
           << double(nsendrecvs)/double(MPI_processes)
           << " processors." << endl
           << "IWS Schedule: Maximum number of communicating partners: "
           << nsendrecvg << " processors." << endl ;
    }
#endif
  }
}

#ifdef TEST

// Dummy work routine
inline void work(int &output, int input){
  if(input == 0)
    cout << "zeroinput!" << endl ;
  usleep(input * 100);
  output = input *10 ;
}

int main(int argc, char** argv[])
{
  MPI_Init(&argc, argv);
  int p;
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  int MPI_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &MPI_rank);

  // Setup Data
  vector<int> node_data;
  for(int i = MPI_rank*100; i < MPI_rank*100 + 100; ++i)
    node_data.push_back(i+1);
  
  // Create iterative weighting schedule class
  gridMotion::IWSSchedule iwss;
  iwss.initialize(MPI_COMM_WORLD,node_data.size()) ;
  
  int timesteps = 20;
  vector<int > output_data(node_data.size()) ; ;
  for(int i = 0; i < timesteps; ++i) {
    // execute work using IWS load balancing
    iterative_weighted_static(iwss,&output_data[0],&node_data[0],work) ;
    // Check that returned values were consistent
    for(int i=0;i<node_data.size();++i)
      if(output_data[i] != node_data[i]*10) {
        cerr << "comm error" << endl ;
        cout << MPI_rank << "- i=" << i << "out[i]=" << output_data[i]
             << "in[i]" << node_data[i] << endl ;
        exit(-1) ;
      }
    // zero out output data so we can find a mistake next round
    for(int i=0;i<node_data.size();++i)
      output_data[i] = -1 ;
  }
  
  MPI_Finalize() ;
}
#endif
