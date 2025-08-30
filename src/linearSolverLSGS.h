#ifndef LGSI_H
#define LGSI_H

#include <vector>
#include <Loci>
#include "flowTypes.h"

namespace LGSI {
  using flowPsi::real_fj ;
  
  struct lineLayout {
    int eqnstart,eqnend,Lstart,Ustart,OuterStart,OuterEnd ;
  } ;

  struct matrixStructure {
    int access_counter ;
    int nmyeqn ;
    // equation variable map
    std::vector<int> eql2g ; // Local to global map for equation information
    // line information
    std::vector<lineLayout> lineData ;
    // maximum line length
    int mxlinelen ;
    // Global to local map for off diagonal matrix information
    Loci::const_Map globalL2l,globalU2l ;
    // connectivity information
    std::vector<pair<int,int> > outerCLR ; // outer face left/right eqn numbers
    //
    // Communication structure
    std::vector<int> send_processors ;
    std::vector<int> send_offsets ;
    std::vector<int> send_entities ;
    std::vector<int> send_entitiesl ;
    std::vector<int> recv_processors ;
    std::vector<int> recv_offsets ;
    std::vector<int> recv_entities ;
    std::vector<int> recv_entitiesl ;

    std::vector<pair<int,int> > periodicCopy ;
    matrixStructure(const matrixStructure &ms) {
      access_counter = ms.access_counter ;
      nmyeqn = ms.nmyeqn ;
      eql2g = ms.eql2g ;
      lineData = ms.lineData ;
      mxlinelen = ms.mxlinelen ;
      outerCLR = ms.outerCLR ;
      globalL2l.setRep(ms.globalL2l.Rep()) ;
      globalU2l.setRep(ms.globalU2l.Rep()) ;
      send_processors = ms.send_processors ;
      send_offsets = ms.send_offsets ;
      send_entities = ms.send_entities ;
      send_entitiesl = ms.send_entitiesl ;
      recv_processors = ms.recv_processors ;
      recv_offsets = ms.recv_offsets ;
      recv_entities = ms.recv_entities ;
      recv_entitiesl = ms.recv_entitiesl ;
      periodicCopy = ms.periodicCopy ;
    }
    matrixStructure &operator=(const matrixStructure &ms) {
      access_counter = ms.access_counter ;
      nmyeqn = ms.nmyeqn ;
      eql2g = ms.eql2g ;
      lineData = ms.lineData ;
      mxlinelen = ms.mxlinelen ;
      outerCLR = ms.outerCLR ;
      globalL2l.setRep(ms.globalL2l.Rep()) ;
      globalU2l.setRep(ms.globalU2l.Rep()) ;
      send_processors = ms.send_processors ;
      send_offsets = ms.send_offsets ;
      send_entities = ms.send_entities ;
      send_entitiesl = ms.send_entitiesl ;
      recv_processors = ms.recv_processors ;
      recv_offsets = ms.recv_offsets ;
      recv_entities = ms.recv_entities ;
      recv_entitiesl = ms.recv_entitiesl ;
      periodicCopy = ms.periodicCopy ;
      return *this ;
    }
    matrixStructure() { access_counter = -1 ;}
  } ;


  template <class REAL> struct scalarMatrixData {
    std::vector<REAL> D ;
    std::vector<REAL> F ;
  } ;

  template <class REAL> struct scalarRHSData {
    std::vector<REAL> B ;
  } ;

  template <class REAL, class REALFJ> void
  LSGSScalarSolve(std::vector<REAL> &result,
                  const struct matrixStructure & restrict matrix,
                  const struct scalarMatrixData<REALFJ> & restrict matdata,
                  const struct std::vector<REAL> & restrict rhsdata,
                  double LSGSRelaxation,
                  double LSGSAbsTol,
                  double LSGSRelTol,
                  int LSGSMaxIter,
                  int LSGSDiagnosticLevel,
                  std::string LSGSMatrixName) {
    using std::cout ;
    using std::cerr ;
    using std::endl ;
    
    // initialize x to zero
    std::vector<REAL> x1(matrix.eql2g.size(),0.0), xt(matrix.nmyeqn) ;

    REAL * X = &x1[0] ;
    std::vector<REAL> rhs(matrix.mxlinelen) ;

    for(size_t l=0;l<matrix.lineData.size();++l) { // loop over lines
      const REALFJ *L = &matdata.F[matrix.lineData[l].Lstart] ;
      const REALFJ *U = &matdata.F[matrix.lineData[l].Ustart] ;
      const REALFJ *D = &matdata.D[matrix.lineData[l].eqnstart] ;
      const REAL *B = &rhsdata[matrix.lineData[l].eqnstart] ;
      
      const int lsz = matrix.lineData[l].eqnend-matrix.lineData[l].eqnstart ;
      
      const int eqb = matrix.lineData[l].eqnstart ;

      // First compute rhs ;
      for(int i=0;i<lsz;++i)
        rhs[i] = B[i] ;

      REAL *restrict x = &X[eqb] ;
      // now solve line
      // forward solve 
      x[0] = rhs[0] ;
      for(int i=1;i<lsz;++i) {
        x[i] = rhs[i] - L[i-1]*D[i-1]*x[i-1] ;
      }
      // backward solve
      x[lsz-1] = D[lsz-1]*x[lsz-1] ;
      for(int i=lsz-2;i>=0;--i) {
        x[i] = D[i]*(x[i] - x[i+1]*U[i]) ;
      }
    }
    std::vector<REAL> send_data(matrix.send_entitiesl.size()) ;
    std::vector<REAL> recv_data(matrix.recv_entitiesl.size()) ;
    std::vector<MPI_Request> requests(matrix.recv_processors.size()) ;
    std::vector<MPI_Status> status(matrix.recv_processors.size()) ;

    double zero_resid = 1 ;
    REAL w = LSGSRelaxation ;
    int num_iter = LSGSMaxIter ;
    // do num_iter steps
    for(int iter = 0;iter < num_iter;++iter) {
      
      int psz = matrix.periodicCopy.size() ;
      for(int i=0;i<psz;++i)
        X[matrix.periodicCopy[i].first] = X[matrix.periodicCopy[i].second] ;
      if(iter == 0) {
        // copy send data for first iteration, aftwards the unrelaxed
        // values are copied later in the loop.
        for(size_t i=0;i<send_data.size();++i)
          send_data[i] = X[matrix.send_entitiesl[i]] ;
      }
      
      for(size_t i=0;i<matrix.recv_processors.size();++i) {
        int rs = matrix.recv_offsets[i] ;
        int re = matrix.recv_offsets[i+1] ;
        int rsz = re-rs ;
        int p = matrix.recv_processors[i] ;
        int tsz = sizeof(REAL) ;
        MPI_Irecv(&recv_data[rs],rsz*tsz,MPI_BYTE,p,99,
                  MPI_COMM_WORLD,&requests[i]) ;
      }
      for(size_t i=0;i<matrix.send_processors.size();++i) {
        int ss = matrix.send_offsets[i] ;
        int se = matrix.send_offsets[i+1] ;
        int ssz = se-ss ;
        int p = matrix.send_processors[i] ;
        int tsz = sizeof(REAL) ;
        MPI_Send(&send_data[ss],ssz*tsz,MPI_BYTE,p,99,MPI_COMM_WORLD) ;
      }
      MPI_Waitall(requests.size(),&requests[0],&status[0]) ;
      // copy recieved data 
      for(size_t i=0;i<recv_data.size();++i)
        X[matrix.recv_entitiesl[i]] = recv_data[i] ;
      for(int i=0;i<psz;++i)
        X[matrix.periodicCopy[i].first] = X[matrix.periodicCopy[i].second] ;
      


      // Forward sweep
      for(size_t l=0;l<matrix.lineData.size();++l) { // loop over lines
	const REALFJ *L = &matdata.F[matrix.lineData[l].Lstart] ;
	const REALFJ *U = &matdata.F[matrix.lineData[l].Ustart] ;
	const REALFJ *D = &matdata.D[matrix.lineData[l].eqnstart] ;
	const REAL *B = &rhsdata[matrix.lineData[l].eqnstart] ;
	const REAL *F = &matdata.F[0] ;

        const int ls = matrix.lineData[l].eqnstart ;
        const int le = matrix.lineData[l].eqnend ;

        for(int i=ls;i<le;++i)
          xt[i] = X[i] ;
        
	const int lsz = le-ls ;
	

	for(int i=0;i<lsz;++i) 
	  rhs[i] = B[i] ;


	const int   os = matrix.lineData[l].OuterStart ;
	const int   oe = matrix.lineData[l].OuterEnd ;
	for(int i=os;i<oe;++i) {
	  rhs[matrix.outerCLR[i].first-ls] 
	    -= X[matrix.outerCLR[i].second]*F[i] ;
	}
      
	REAL *restrict x = &X[ls] ;
	// now solve line
	// forward solve 
	x[0] = rhs[0] ;
	for(int i=1;i<lsz;++i) {
	  x[i] = rhs[i] - L[i-1]*D[i-1]*x[i-1] ;
	}
	// backward solve
	x[lsz-1] = D[lsz-1]*x[lsz-1] ;
	for(int i=lsz-2;i>=0;--i) {
	  x[i] = D[i]*(x[i] - x[i+1]*U[i]) ;
	}
      }

      // Backward sweep
      for(int l=matrix.lineData.size()-1;l>=0;--l) { // loop over lines
	const REALFJ *L = &matdata.F[matrix.lineData[l].Lstart] ;
	const REALFJ *U = &matdata.F[matrix.lineData[l].Ustart] ;
	const REALFJ *D = &matdata.D[matrix.lineData[l].eqnstart] ;
	const REAL *B = &rhsdata[matrix.lineData[l].eqnstart] ;
	const REAL *F = &matdata.F[0] ;

	const int lsz = matrix.lineData[l].eqnend-matrix.lineData[l].eqnstart ;

	const int ls = matrix.lineData[l].eqnstart ;

	// First compute rhs ;
	for(int i=0;i<lsz;++i)
	  rhs[i] = B[i] ;

	const int   os = matrix.lineData[l].OuterStart ;
	const int   oe = matrix.lineData[l].OuterEnd ;
	for(int i=os;i<oe;++i) {
	  rhs[matrix.outerCLR[i].first-ls] -= X[matrix.outerCLR[i].second]*F[i] ;
	}

	REAL *restrict x = &X[ls] ;
	// now solve line
	// forward solve 
	x[0] = rhs[0] ;
	for(int i=1;i<lsz;++i) {
	  x[i] = rhs[i] - L[i-1]*D[i-1]*x[i-1] ;
	}
	// backward solve
	x[lsz-1] = D[lsz-1]*x[lsz-1] ;
	for(int i=lsz-2;i>=0;--i) {
	  x[i] = D[i]*(x[i] - x[i+1]*U[i]) ;
	}
      }
      // copy data to send before relaxation
      for(size_t i=0;i<send_data.size();++i)
        send_data[i] = X[matrix.send_entitiesl[i]] ;

      // Relax update and compute current estimate for residual
      double residl = 0 ;
      for(int i=0;i<matrix.nmyeqn;++i) {
        X[i] = w*X[i] + (1.-w)*xt[i] ;
	double r = (X[i]-xt[i]) ;
	residl += r*r ;
      }
      double resid = 0 ;
      MPI_Allreduce(&residl,&resid,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD) ;
      resid = sqrt(resid) ;
      if(0 == iter)
	zero_resid = max(resid,1e-16) ;
      double rel_resid = resid/zero_resid ;
      //      if(Loci::MPI_rank == 0)
      //	cout << iter << ":sresid = " << rel_resid << endl ;
      if(rel_resid > 1.0)
	w *= 1./rel_resid ;

      if(LSGSDiagnosticLevel > 0 && Loci::MPI_rank == 0) {
        cout << "LSGS Matrix " << LSGSMatrixName 
             << ", iter=" << iter+1 
             << ", resid=" << resid << ", rel_resid=" << rel_resid << endl ;
      }
      
      if(resid < LSGSAbsTol)
	break ;
         if(rel_resid < LSGSRelTol)
	break ;

    }
    
    // copy data out
    result.swap(x1) ;  
    //    for(int i=0;i<matrix.nmyeqn;++i) {
    //      $LSGSScalarSolve(X,B)[matrix.eql2g[i]] = x1[i] ; 
    //    }
  }


  template <class REAL, class REALFJ> void
  LSGSScalarMatVecProd(std::vector<REAL> & result,
                       const struct matrixStructure & matrix,
                       const struct scalarMatrixData<REALFJ> & matdata,
                       struct std::vector<REAL> & restrict X) {
    std::vector<REAL> send_data(matrix.send_entitiesl.size()) ;
    std::vector<REAL> recv_data(matrix.recv_entitiesl.size()) ;
    std::vector<MPI_Request> requests(matrix.recv_processors.size()) ;
    std::vector<MPI_Status> status(matrix.recv_processors.size()) ;
 
    std::vector<REAL> Y(matrix.eql2g.size(),0.0) ;
    FATAL(X.size() != matrix.eql2g.size()) ;
      

    int psz = matrix.periodicCopy.size() ;
    for(int i=0;i<psz;++i)
      X[matrix.periodicCopy[i].first] = X[matrix.periodicCopy[i].second] ;
    for(size_t i=0;i<send_data.size();++i)
      send_data[i] = X[matrix.send_entitiesl[i]] ;
      
    for(size_t i=0;i<matrix.recv_processors.size();++i) {
      int rs = matrix.recv_offsets[i] ;
      int re = matrix.recv_offsets[i+1] ;
      int rsz = re-rs ;
      int p = matrix.recv_processors[i] ;
      int tsz = sizeof(REAL) ;
      MPI_Irecv(&recv_data[rs],rsz*tsz,MPI_BYTE,p,99,
                MPI_COMM_WORLD,&requests[i]) ;
    }
    for(size_t i=0;i<matrix.send_processors.size();++i) {
      int ss = matrix.send_offsets[i] ;
      int se = matrix.send_offsets[i+1] ;
      int ssz = se-ss ;
      int p = matrix.send_processors[i] ;
      int tsz = sizeof(REAL) ;
      MPI_Send(&send_data[ss],ssz*tsz,MPI_BYTE,p,99,MPI_COMM_WORLD) ;
    }
    MPI_Waitall(requests.size(),&requests[0],&status[0]) ;
    // copy recieved data 
    for(size_t i=0;i<recv_data.size();++i)
      X[matrix.recv_entitiesl[i]] = recv_data[i] ;
    for(int i=0;i<psz;++i)
      X[matrix.periodicCopy[i].first] = X[matrix.periodicCopy[i].second] ;

    // Diagonal contribution
    for(int i=0;i<matrix.nmyeqn;++i) 
      Y[i] = X[i]/matdata.D[i] ;
    // Off-diagonal contributions (not part of lines)
    for(size_t i=0;i<matrix.outerCLR.size();++i)
      Y[matrix.outerCLR[i].first] += X[matrix.outerCLR[i].second]*matdata.F[i] ;

    // Add in line contributions
    for(size_t l=0;l<matrix.lineData.size();++l) { // loop over lines
      const REALFJ *L = &matdata.F[matrix.lineData[l].Lstart] ;
      const REALFJ *U = &matdata.F[matrix.lineData[l].Ustart] ;

      const int ls = matrix.lineData[l].eqnstart ;
      const int le = matrix.lineData[l].eqnend ;
      
      const int lsz = le-ls ;
      const int offset = matrix.lineData[l].eqnstart ;
      for(int i=0;i<lsz-1;++i) {
        Y[i+offset+1] += X[offset+i]*L[i-1] ;
        Y[i+offset] += X[offset+i+1]*U[i] ;
      }

    }
    result.swap(Y) ;
  }

  template<class A, class B>
  double scalarInnerProduct(const std::vector<A> &v1,
                            const std::vector<B> &v2,
                            int locsize) {
    double sum = 0.0 ;
    for(int i=0;i<locsize;++i)
      sum += v1[i]*v2[i] ;
    double result = 0.0 ;
    MPI_Allreduce(&sum, &result,1, MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD) ;
    return result ;
  }

  class matrixHelpers: public Loci::CPTR_type {
  public:
    virtual void preprocessLine(double *D, 
				const real_fj *L, const real_fj *U, 
				int lsz, int vs) const = 0 ;
    virtual void linerhs(double *rhs, const double *X, const real_fj *B, 
			const pair<int,int> *CLR, const real_fj *F, 
			int ls, int le, int os, int oe, int vs) const = 0;
    virtual void lineSolve(double *X, const double *rhs, 
			   const double *D,
			   const real_fj *L, const real_fj *U, 
			   int lsz, int vs) const = 0 ;
  } ;

  struct blockedMatrixData {
    Loci::CPTR<matrixHelpers> matrixHelp ;
    int vs,vs2  ;
    std::vector<double> D ;
    std::vector<real_fj> F ;
  } ;
  template<class REAL_FJ> struct blockedRHSData {
    std::vector<REAL_FJ> B ;
  } ;

  template <class REAL, class REALV> void
  LSGSBlockedMatVecProd(std::vector<REAL> & result,
                       const struct matrixStructure & matrix,
                       const struct blockedMatrixData & matdata,
                       struct std::vector<REALV> & restrict X) {
    int vs = matdata.vs ;
    int vs2 = vs*vs ;
    int send_sz = matrix.send_entitiesl.size() ;
    int recv_sz = matrix.recv_entitiesl.size() ;
    std::vector<REALV> send_data(send_sz*vs) ;
    std::vector<REALV> recv_data(recv_sz*vs) ;
    std::vector<MPI_Request> requests(matrix.recv_processors.size()) ;
    std::vector<MPI_Status> status(matrix.recv_processors.size()) ;


    vector<double> scratch(vs) ;
    std::vector<REAL> Y(matrix.eql2g.size()*vs,0.0) ;
    FATAL(X.size() != matrix.eql2g.size()*vs) ;
      

    int psz = matrix.periodicCopy.size() ;
    for(int i=0;i<psz;++i)
      for(int j=0;j<vs;++j)
        X[matrix.periodicCopy[i].first*vs+j] =
          X[matrix.periodicCopy[i].second*vs+j] ;
    for(int i=0;i<send_sz;++i)
      for(int j=0;j<vs;++j)
        send_data[i*vs+j] = X[matrix.send_entitiesl[i]*vs+j] ;
      
    for(size_t i=0;i<matrix.recv_processors.size();++i) {
      int rs = matrix.recv_offsets[i]*vs ;
      int re = matrix.recv_offsets[i+1]*vs ;
      int rsz = re-rs ;
      int p = matrix.recv_processors[i] ;
      int tsz = sizeof(REALV) ;
      MPI_Irecv(&recv_data[rs],rsz*tsz,MPI_BYTE,p,99,
                MPI_COMM_WORLD,&requests[i]) ;
    }
    for(size_t i=0;i<matrix.send_processors.size();++i) {
      int ss = matrix.send_offsets[i]*vs ;
      int se = matrix.send_offsets[i+1]*vs ;
      int ssz = se-ss ;
      int p = matrix.send_processors[i] ;
      int tsz = sizeof(REALV) ;
      MPI_Send(&send_data[ss],ssz*tsz,MPI_BYTE,p,99,MPI_COMM_WORLD) ;
    }
    MPI_Waitall(requests.size(),&requests[0],&status[0]) ;
    // copy recieved data 
    for(int i=0;i<recv_sz;++i)
      for(int j=0;j<vs;++j)
        X[matrix.recv_entitiesl[i]*vs+j] = recv_data[i*vs+j] ;
    for(int i=0;i<psz;++i)
      for(int j=0;j<vs;++j)
        X[matrix.periodicCopy[i].first*vs+j] = X[matrix.periodicCopy[i].second*vs+j] ;

    // Diagonal contribution (Note D is LU factorized
    for(int n=0;n<matrix.nmyeqn;++n) {
      const_Mat<double> Dm(&matdata.D[n*vs2],vs) ;
      const REALV *x = &X[n*vs] ;
      // Perform scratch =U*X
      for(int i=0;i<vs;++i) {
        double tmp = 0 ;
        for(int j=i;j<vs;++j)
          tmp += Dm[i][j]*x[j] ;
        scratch[i] = tmp ;
      }
      // Peform Y = L*scratch
      // first the diagonal
      for(int i=0;i<vs;++i)
        Y[n*vs+i] = scratch[i] ;
      for(int i=1;i<vs;++i)
        for(int j=0;j<i;++j)
          Y[n*vs+i] += Dm[i][j]*scratch[j] ;
    }
    // Off-diagonal contributions (not part of lines)
    for(size_t i=0;i<matrix.outerCLR.size();++i) {
      const_Mat<real_fj> Fm(&matdata.F[i*vs2],vs) ;
      size_t first = matrix.outerCLR[i].first*vs ;
      size_t second = matrix.outerCLR[i].second*vs ;
      Fm.dotprod_accum(&X[second],&Y[first]) ;
      //      Y[matrix.outerCLR[i].first] += X[matrix.outerCLR[i].second]*matdata.F[i] ;
    }
    // Add in line contributions
    for(size_t l=0;l<matrix.lineData.size();++l) { // loop over lines
      const real_fj *L = &matdata.F[matrix.lineData[l].Lstart*vs] ;
      const real_fj *U = &matdata.F[matrix.lineData[l].Ustart*vs] ;

      const int ls = matrix.lineData[l].eqnstart ;
      const int le = matrix.lineData[l].eqnend ;
      
      const int lsz = le-ls ;
      const int offset = matrix.lineData[l].eqnstart ;
      for(int i=0;i<lsz-1;++i) {
        const_Mat<real_fj> Lm(L+(i-1)*vs2,vs) ;
        Lm.dotprod_accum(&X[(offset+i)*vs],&Y[(i+offset+1)*vs]) ;
        //Y[i+offset+1] += X[offset+i]*L[i-1] ;
        const_Mat<real_fj> Um(U+i*vs2,vs) ;
        Um.dotprod_accum(&X[(offset+i+1)*vs],&Y[(i+offset)*vs]) ;
        //Y[i+offset] += X[offset+i+1]*U[i] ;
      }

    }
    result.swap(Y) ;
  }
  


  template <class REAL, class REALV> void
  LSGSBlockSolve(std::vector<REAL> &result,
                 const struct matrixStructure & restrict matrix,
                 const struct blockedMatrixData & restrict matdata,
                 const struct std::vector<REALV> & restrict rhsdata,
                 double LSGSRelaxation,
                 double LSGSAbsTol,
                 double LSGSRelTol,
                 int LSGSMaxIter,
                 int LSGSDiagnosticLevel,
                 std::string LSGSMatrixName) {
    using std::cout ;
    using std::cerr ;
    using std::endl ;
    
    int vs = matdata.vs ;
    int vs2 = vs*vs ;
    int vsize = matrix.eql2g.size() ;
    Loci::CPTR<matrixHelpers> matrixHelp = matdata.matrixHelp ;
    
    vector<REAL> x1(vsize*vs) , xt(matrix.nmyeqn*vs) ;
    
    REAL *X = &x1[0] ;
    
    vector<double> rhs(matrix.mxlinelen*vs) ;
    double zero_resid =1;
    // jacobi initializes iteration
    for(size_t l=0;l<matrix.lineData.size();++l) { // loop over lines
      const REALV *B = &rhsdata[0] ;
      int ls = matrix.lineData[l].eqnstart ;
      int le = matrix.lineData[l].eqnend ;

      double * rhsg = &rhs[0] - ls*vs ; 
      for(int i=ls*vs;i<le*vs;++i)
	rhsg[i] = -B[i] ;

      const double *D = &(matdata.D[ls*vs2]) ;
      const real_fj *L = &matdata.F[matrix.lineData[l].Lstart*vs2] ;
      const real_fj *U = &matdata.F[matrix.lineData[l].Ustart*vs2] ;
      int lsz = le-ls ;
      matrixHelp->lineSolve(X+ls*vs,&rhs[0],D,L,U,lsz,vs) ;
    }

    vector<double> send_data(matrix.send_entitiesl.size()*vs) ;
    vector<double> recv_data(matrix.recv_entitiesl.size()*vs) ;
    vector<MPI_Request> requests(matrix.recv_processors.size()+matrix.send_processors.size()) ;
    vector<MPI_Status> status(matrix.recv_processors.size()+matrix.send_processors.size()) ;
    int num_iter = LSGSMaxIter ;
    double w = LSGSRelaxation ;

    // Entry communication 
    for(int iter = 0;iter<num_iter;++iter) {
      // copy send data 
      int psz = matrix.periodicCopy.size() ;
      if(iter == 0) {
        for(int i=0;i<psz;++i) {
          for(int j=0;j<vs;++j) {
            X[matrix.periodicCopy[i].first*vs+j] = X[matrix.periodicCopy[i].second*vs+j] ;
          }
        }
        for(size_t i=0;i<matrix.send_entitiesl.size();++i) {
          int loc = matrix.send_entitiesl[i]*vs ;
          int ii = i*vs ;
          for(int j=0;j<vs;++j)
            send_data[ii+j] = X[loc+j] ;
        }
        for(size_t i=0;i<matrix.recv_processors.size();++i) {
          int rs = matrix.recv_offsets[i]*vs ;
          int re = matrix.recv_offsets[i+1]*vs ;
          int rsz = re-rs ;
          int p = matrix.recv_processors[i] ;
          MPI_Irecv(&recv_data[rs],rsz,MPI_DOUBLE,p,99,
                    MPI_COMM_WORLD,&requests[i]) ;
        }

        int rpsz = matrix.recv_processors.size() ;
        
        for(size_t i=0;i<matrix.send_processors.size();++i) {
          int ss = matrix.send_offsets[i]*vs ;
          int se = matrix.send_offsets[i+1]*vs ;
          int ssz = se-ss ;
          int p = matrix.send_processors[i] ;
          MPI_Isend(&send_data[ss],ssz,MPI_DOUBLE,p,99,MPI_COMM_WORLD,
                    &requests[i+rpsz]) ;
        }
        MPI_Waitall(requests.size(),&requests[0],&status[0]) ;
        // copy recieved data 
        for(size_t i=0;i<matrix.recv_entitiesl.size();++i) {
          int loc = matrix.recv_entitiesl[i]*vs ;
          int ii = i*vs ;
          for(int j=0;j<vs;++j)
            X[loc+j] = recv_data[ii+j] ;
        }
      }
      for(int i=0;i<psz;++i) {
        for(int j=0;j<vs;++j) {
          X[matrix.periodicCopy[i].first*vs+j] = X[matrix.periodicCopy[i].second*vs+j] ;
        }
      }

      // Forward sweep
      for(size_t l=0;l<matrix.lineData.size();++l) { // loop over lines
	const REALV *B = &(rhsdata[0]) ;
	const real_fj *F = &(matdata.F[0]) ;
	const pair<int,int> *CLR = &(matrix.outerCLR[0]) ;
	int ls = matrix.lineData[l].eqnstart ;
	int le = matrix.lineData[l].eqnend ;
	int os = matrix.lineData[l].OuterStart ;
	int oe = matrix.lineData[l].OuterEnd ;
	double *restrict rhsg = &rhs[0] - ls*vs ; 
	matrixHelp->linerhs(rhsg,&X[0],B,CLR,F,ls,le,os,oe,vs) ;

	const double *D = &(matdata.D[ls*vs2]) ;
	const real_fj *L = &matdata.F[matrix.lineData[l].Lstart*vs2] ;
	const real_fj *U = &matdata.F[matrix.lineData[l].Ustart*vs2] ;
	int lsz = le-ls ;
	for(int i=ls*vs;i<le*vs;++i)
	  xt[i] = X[i] ;

	matrixHelp->lineSolve(X+ls*vs,&rhs[0],D,L,U,lsz,vs) ;
      }

      // Backward sweep
      for(int l=matrix.lineData.size()-1;l>=0;--l) { // loop over lines
	const REALV *B = &(rhsdata[0]) ;
	const real_fj *F = &(matdata.F[0]) ;
	const pair<int,int> *CLR = &(matrix.outerCLR[0]) ;
	int ls = matrix.lineData[l].eqnstart ;
	int le = matrix.lineData[l].eqnend ;
	int os = matrix.lineData[l].OuterStart ;
	int oe = matrix.lineData[l].OuterEnd ;
	double * rhsg = &rhs[0] - ls*vs ; 
	matrixHelp->linerhs(rhsg,X,B,CLR,F,ls,le,os,oe,vs) ;

	const double *D = &(matdata.D[ls*vs2]) ;
	const real_fj *L = &matdata.F[matrix.lineData[l].Lstart*vs2] ;
	const real_fj *U = &matdata.F[matrix.lineData[l].Ustart*vs2] ;
	int lsz = le-ls ;
	for(int i=ls*vs;i<le*vs;++i)
	  xt[i] = X[i] ;
	matrixHelp->lineSolve(X+ls*vs,&rhs[0],D,L,U,lsz,vs) ;
      }
      
      // Copy data for communication
      for(int i=0;i<psz;++i) {
        for(int j=0;j<vs;++j) {
          X[matrix.periodicCopy[i].first*vs+j] = X[matrix.periodicCopy[i].second*vs+j] ;
        }
      }
      for(size_t i=0;i<matrix.send_entitiesl.size();++i) {
        int loc = matrix.send_entitiesl[i]*vs ;
        int ii = i*vs ;
        for(int j=0;j<vs;++j)
          send_data[ii+j] = X[loc+j] ;
      }
      
      // Initiate interprocessor communication
      for(size_t i=0;i<matrix.recv_processors.size();++i) {
        int rs = matrix.recv_offsets[i]*vs ;
        int re = matrix.recv_offsets[i+1]*vs ;
        int rsz = re-rs ;
        int p = matrix.recv_processors[i] ;
        MPI_Irecv(&recv_data[rs],rsz,MPI_DOUBLE,p,99,
                  MPI_COMM_WORLD,&requests[i]) ;
      }

      int rpsz = matrix.recv_processors.size() ;
        
      for(size_t i=0;i<matrix.send_processors.size();++i) {
        int ss = matrix.send_offsets[i]*vs ;
        int se = matrix.send_offsets[i+1]*vs ;
        int ssz = se-ss ;
        int p = matrix.send_processors[i] ;
        MPI_Isend(&send_data[ss],ssz,MPI_DOUBLE,p,99,MPI_COMM_WORLD,
                  &requests[i+rpsz]) ;
      }

      // Compute residual
      double residl = 0 ;
      for(int i=0;i<matrix.nmyeqn*vs;++i) {
	X[i] = w*X[i] + (1.-w)*xt[i] ;
	double r = X[i]-xt[i] ;
	residl += r*r ;
      }

      // Finish communication after local line solve
      MPI_Waitall(requests.size(),&requests[0],&status[0]) ;
      // copy recieved data 
      for(size_t i=0;i<matrix.recv_entitiesl.size();++i) {
        int loc = matrix.recv_entitiesl[i]*vs ;
        int ii = i*vs ;
        for(int j=0;j<vs;++j)
          X[loc+j] = recv_data[ii+j] ;
      }

      for(int i=0;i<psz;++i) {
        for(int j=0;j<vs;++j) {
          X[matrix.periodicCopy[i].first*vs+j] = X[matrix.periodicCopy[i].second*vs+j] ;
        }
      }

      // Check residual
      double resid = 0 ;
      MPI_Allreduce(&residl,&resid,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD) ;
      resid = sqrt(resid) ;
      if(0 == iter)
	zero_resid = max(resid,1e-16) ;
      double rel_resid = resid/zero_resid ;
      if(rel_resid > 1.0)
	w *= 1./rel_resid ;

      if(LSGSDiagnosticLevel > 0 && Loci::MPI_rank == 0 ) {
        cout << "LSGS Matrix " << LSGSMatrixName
             << ", iter=" << iter+1 
             << ", resid=" << resid << ", rel_resid=" << rel_resid << endl ;

      }
      if(LSGSDiagnosticLevel > 1) {
        vector<double> resid ;
        LSGSBlockedMatVecProd(resid,matrix,matdata,x1) ;
        double sum = 0 ;
        for(int i=0;i<matrix.nmyeqn*vs;++i) {
          double diff = resid[i]-rhsdata[i] ;
          sum += diff*diff ;
        }
        double gsum = 0 ;
        MPI_Allreduce(&gsum,&sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD) ;
        if(Loci::MPI_rank == 0) 
          cout << "true residual = " << sqrt(gsum) << endl ;
      }
      
      if(resid < LSGSAbsTol) 
	break ;
      if(rel_resid < LSGSRelTol)
	break ;
    }
    result.swap(x1) ;
  }

}
#endif
