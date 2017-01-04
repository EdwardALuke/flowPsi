//#############################################################################
//#
//# Copyright 2015, Mississippi State University
//#
//# This file is part of the CHEM solver framework.
//#
//# The CHEM solver framework is free software: you can redistribute it 
//# and/or modify it under the terms of the Lesser GNU General Public License 
//# as published by the Free Software Foundation, either version 3 of the 
//# License, or (at your option) any later version.
//#
//# The CHEM solver is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# Lesser GNU General Public License for more details.
//#
//# You should have received a copy of the Lesser GNU General Public License
//# along with the CHEM solver.  If not, see <http://www.gnu.org/licenses>
//#
//# NOTICE: The GNU license describes the parameters of your redistribution
//# rights granted by Mississippi State University but the redistribution
//# of this software is also constrained by export control laws.  This
//# software may also be considered as covered by EAR or ITAR regulation
//# regimes and so it is your responsibility to ensure that this software
//# is redistributed only to US persons and any redistribution is in
//# conformance with applicable United States export control laws.
//#
//#############################################################################
#include <Loci.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>
#include <algorithm>
#include <string>
#include <vector>
#include <string>
#include <sys/types.h>
#include <sys/stat.h>

using std::string ;

//using namespace std ;
using std::ifstream ;
using std::ofstream ;
using std::ostringstream ;
using std::cout ;
using std::cin ;
using std::endl ;
using std::cerr ;
using std::vector ;
using std::ios ;

typedef vector3d<double> vect3d ;

void dump_scalar(store<float> &c2n, string sname, string iter, string casename) {
  string filename = "output/" + sname + "_sca." + iter + "_" + casename ;
  
  if(Loci::MPI_rank == 0)
    cout << "writing file " << filename << endl ;
  
  fact_db facts ;

  hid_t file_id = Loci::hdf5CreateFile(filename.c_str(),H5F_ACC_TRUNC,
                                       H5P_DEFAULT, H5P_DEFAULT) ;
  
  Loci::writeContainer(file_id,sname,c2n.Rep(),facts) ;

  Loci::hdf5CloseFile(file_id) ;
}

entitySet read_domain(string sname, string iter, string casename) {
  
  string filename = "output/" + sname + "_sca." + iter + "_" + casename ;
  hid_t file_id, group_id ;
  file_id = H5Fopen(filename.c_str(),H5F_ACC_RDONLY, H5P_DEFAULT) ;
#ifdef H5_USE_16_API
  group_id = H5Gopen(file_id,sname.c_str()) ;
#else
  group_id = H5Gopen(file_id,sname.c_str(),H5P_DEFAULT) ;
#endif

  entitySet valdom ;
  Loci::HDF5_ReadDomain(group_id,valdom);

  H5Gclose(group_id) ;
  H5Fclose(file_id) ;
  return valdom ;
}

void read_scalar(store<float> &var, string sname, string iter, string casename) {
  string filename = "output/" + sname + "_sca." + iter + "_" + casename ;
  hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
                                     H5F_ACC_RDONLY,
                                     H5P_DEFAULT) ;
  if(file_id < 0) {
    cerr << "unable to open file '" << filename << "'!" << endl ;
    return ;
  }

  fact_db facts ;
  Loci::readContainer(file_id,sname,var.Rep(),EMPTY,facts) ;
  Loci::hdf5CloseFile(file_id) ;
}

void Usage() {
  cout << "Usage:" << endl
       << "adpt <iter> <sensitivity> <casename> <var1> <var2> ..." << endl
       << " where <iter> is the interation number," << endl
       << " <sensitivity> is the error sensitivity (std deviations from mean)," << endl
       << " <casename> is the name of the case file" << endl
       << " <var1> is the first error variable, there can be 1 or more"<< endl
       << "        error variables listed."<< endl << endl ;
}
string getPosFile(string output_dir,string iteration, string casename) {
  string posname = output_dir+"/grid_pos." + iteration + "_" + casename ;
  struct stat tmpstat ;
  if(stat(posname.c_str(),&tmpstat) != 0) {
    posname = output_dir+"/grid_pos." + casename ;
  } else if(tmpstat.st_size == 0) {
    posname = output_dir+"/grid_pos." + casename ;
  }
  return posname ;
}

int main(int ac, char *av[]) {
  Loci::Init(&ac,&av) ;
  
  if(Loci::MPI_processes != 1) {
    cerr << "This program is serial and should only be run on one processor!"
         << endl ;
    Loci::Abort() ;
    exit(-1) ;
  }
      
  double xvalmax = 1e33 ;
  double xvalmin = -1e33 ;
  double yvalmax = 1e33 ;
  double yvalmin = -1e33 ;
  double zvalmax = 1e33 ;
  double zvalmin = -1e33 ;
  
  bool checkx = false ;
  double mindist = 1e-10 ;
  double maxdist = 1e-10 ;
  bool checkmin = false ;
  bool checkmax = false ;
  while(ac > 6 && av[1][0] == '-') {
    if(!strcmp(av[1],"-xmax")) {
      xvalmax = atof(av[2]) ;
      checkx = true ;
      av+=2 ;
      ac-=2 ;
    } else if(!strcmp(av[1],"-xmin")) {
      xvalmin = atof(av[2]) ;
      checkx = true ;
      av+=2 ;
      ac-=2 ;
    } else if(!strcmp(av[1],"-ymax")) {
      yvalmax = atof(av[2]) ;
      checkx = true ;
      av+=2 ;
      ac-=2 ;
    } else if(!strcmp(av[1],"-ymin")) {
      yvalmin = atof(av[2]) ;
      checkx = true ;
      av+=2 ;
      ac-=2 ;
    } else if(!strcmp(av[1],"-zmax")) {
      zvalmax = atof(av[2]) ;
      checkx = true ;
      av+=2 ;
      ac-=2 ;
    } else if(!strcmp(av[1],"-zmin")) {
      zvalmin = atof(av[2]) ;
      checkx = true ;
      av+=2 ;
      ac-=2 ;
    } else if(!strcmp(av[1],"-dmin")) {
      mindist = atof(av[2]) ;
      checkmin = true ;
      av+=2 ;
      ac-=2 ;
    } else if(!strcmp(av[1],"-dmax")) {
      maxdist = atof(av[2]) ;
      checkmax = true ;
      av += 2 ;
      ac -= 2 ;
    }
  }

  if(ac < 5) {
    cerr << "incorrect number of arguments to adpt" << endl ;
    Usage() ;
    exit(-1) ;
  }

  char *iter = av[1] ;
  double sensitivity = atof(av[2]) ;
  char *problem_name = av[3] ;


  entitySet valdom = read_domain(av[4],iter,problem_name) ;

  store<int> flag ;
  flag.allocate(valdom) ;
  for(entitySet::const_iterator ei=valdom.begin(); ei!=valdom.end();++ei) {
    flag[*ei] = 1 ;
  }


  if(checkx) {
    store<vector3d<float> > pos ;
    string posname = getPosFile("output",string(iter),string(problem_name)) ;
    hid_t file_id = Loci::hdf5OpenFile(posname.c_str(),
                                       H5F_ACC_RDONLY,
                                       H5P_DEFAULT) ;
    if(file_id < 0) {
      cerr << "unable to get grid positions for iteration " << iter
           << endl ;
      cerr << "does file '" << posname << "' exist?" << endl ;
      exit(-1) ;
    }

    fact_db facts ;
    Loci::readContainer(file_id,"pos",pos.Rep(),EMPTY,facts) ;
    Loci::hdf5CloseFile(file_id) ;
    
    for(entitySet::const_iterator ei=valdom.begin(); ei!=valdom.end();++ei) {
      if(pos[*ei].x> xvalmax || pos[*ei].x < xvalmin ||
         pos[*ei].y> yvalmax || pos[*ei].y < yvalmin ||
         pos[*ei].z> zvalmax || pos[*ei].z < zvalmin) 
        flag[*ei] = 0 ;
    }
  }

  if(checkmin) {
    store<float> var ;
    read_scalar(var,"dist_noslip",iter,problem_name) ;

    cout << "threshold for the minimum distance to a viscous wall = " << mindist << endl ;
    cout << "counting number of cells eliminated..." << endl ;

    int count = 0 ;
    for(entitySet::const_iterator ei=valdom.begin(); ei!=valdom.end();++ei) {
      if(var[*ei] < mindist) {
        flag[*ei] = 0 ;
        count++ ;
      }
    }
    cout << "number of boundary layer cells eliminated for refinement = " << count << endl ;
  }

  store<float>  adapt;
  adapt.allocate(valdom) ;

  for(entitySet::const_iterator ei=valdom.begin(); ei!=valdom.end();++ei) {
    adapt[*ei] = 0 ;
  }
  if(checkmax) {
    store<float> var ;
    read_scalar(var,"dist_noslip",iter,problem_name) ;

    cout << "threshold for the minimum distance to a viscous wall = " << maxdist << endl ;
    cout << "counting number of cells added..." << endl ;

    int count = 0 ;
    for(entitySet::const_iterator ei=valdom.begin(); ei!=valdom.end();++ei) {
      if(var[*ei] < maxdist) {
        adapt[*ei] = 1 ;
        count++ ;
      }
    }
    cout << "number of boundary layer cells added for refinement = " << count << endl ;
  }


  entitySet dom ;
  for(int i=4;i<ac;++i) {
    store<float> var ;
    read_scalar(var,av[i],iter,problem_name) ;

    double total = 0.0 ;
    int number =  0 ;

    for(entitySet::const_iterator ei=valdom.begin(); ei!=valdom.end();++ei) {
      if(flag[*ei] != 0) {
        total += var[*ei] ;
        number++ ;
      }
    }
    if(number==0)
      number++ ;
    double average = total/double(number) ;

    double sigmatot = 0.0 ;
    for(entitySet::const_iterator ei=valdom.begin(); ei!=valdom.end();++ei) {
      if(flag[*ei] != 0) {
        double diff = var[*ei]-average ;
        sigmatot += diff*diff ;
      }
    }
    
    double sigma = sqrt(sigmatot/double(number)) ;

    cout << "average = " << average << ", sigma = "<< sigma << endl ;
    double threshold =  average + sensitivity*sigma ;
    cout << "threshold = " << threshold << endl ;

    for(entitySet::const_iterator ei=valdom.begin(); ei!=valdom.end();++ei) {
      if(flag[*ei] != 0 && var[*ei] > threshold) 
        adapt[*ei] += 1. ;
    }
  }


  cout << "counting refined nodes" << endl ;
  int number = valdom.size() ;
  int count = 0;
  for(entitySet::const_iterator ei=valdom.begin(); ei!=valdom.end();++ei) {
    if(adapt[*ei] > 0)
      count++ ;
  }

  cout << "percentage refined nodes = " << double(count)/double(number)*100.0
       << endl ;

  dump_scalar(adapt,"refine",iter,problem_name) ;

  ofstream outfile("refine.dat",ios::out) ;
  for(entitySet::const_iterator ei=valdom.begin(); ei!=valdom.end();++ei) {
    if(adapt[*ei] > 0)
      outfile << '1' << endl ;
    else
      outfile << '0' << endl ;
  }
  

  Loci::Finalize() ;
}
                      
