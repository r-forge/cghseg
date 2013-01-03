#ifndef CGHSEG_COMPACTEMINITREG_H
#define CGHSEG_COMPACTEMINITREG_H

#include "numlib.h"

#include <iostream>
#include <limits>
#include <fstream>
#include <cmath>

#define  Dist(i,j) _D[i-2][j-1]
#define  nk(i)     _nk[i-1]
#define  mxk(i)    _mxk[i-1]
#define  myk(i)    _myk[i-1]
#define  mx2k(i)   _mx2k[i-1]
#define  my2k(i)   _my2k[i-1]
#define  mxyk(i)   _mxyk[i-1]
#define  cxyk(i)   _cxyk[i-1]
#define  vxk(i)    _vxk[i-1]
#define  vyk(i)    _vyk[i-1]
#define  Dtmp(i)   _Dtmp[i-1]

namespace cghseg
{
  void compute_compactEM_init_Reg(numlib_vector *xk,numlib_vector *x2k, numlib_vector *yk,numlib_vector *y2k,numlib_vector *xyk,numlib_vector *nk, int P , numlib_vector **phiend, int OMP_NUM_THREADS = 1);
       
  class compactEM_init_Reg{
  private:
    inline void CAHminReg(int k, int& imin, int& jmin, double& dmin);
    inline void CAHcoreReg(int i, int ntmp, double mxtmp, double mx2tmp, double vxtmp, double mytmp, double my2tmp, double vytmp,double mxytmp, double cxytmp);
    inline void CAHcopyReg(int k, int imin, int jmin);
    
  public:
    int       _lengthx;        
    int       _K;              
    int       _P;              
    double   *_phi;            
    double   *_xk;             
    double   *_x2k;            
    double   *_yk;             
    double   *_y2k;            
    double   *_xyk;            
    double   *_mxk;            
    double   *_mx2k;           
    double   *_vxk;            
    double   *_myk;            
    double   *_my2k;           
    double   *_vyk;            
    double   *_mxyk;           
    double   *_cxyk;           
    double   *_nk;             
    double   *_mxk0;           
    double   *_mx2k0;          
    double   *_vxk0;           
    double   *_myk0;           
    double   *_my2k0;          
    double   *_vyk0;           
    double   *_mxyk0;          
    double   *_cxyk0;          
    int      *_nk0;            
    double   **_D;             
    double   *_Dtmp;           
    double   **_tau;           

    compactEM_init_Reg(int nbseg, int nbclust, int OMP_NUM_THREADS=1);
    void CAH_Reg();
    void compute_phi_Reg();
    void Init(double *Dataxk, double *Datax2k,double *Datayk,double *Datay2k,double *Dataxyk,double *Datank);    
    friend std::ostream & operator << (std::ostream &s, const compactEM_init_Reg & compactEMiReg);
    ~compactEM_init_Reg();
  };
  
  void compactEM_init_Reg::CAHcoreReg(int i, int ntmp, double mxtmp, double mx2tmp, double vxtmp, double mytmp, double my2tmp, double vytmp,double mxytmp,double cxytmp){
    double mxpool    = (nk(i)*mxk(i)+ntmp*mxtmp)/(nk(i)+ntmp);
    double mypool    = (nk(i)*myk(i)+ntmp*mytmp)/(nk(i)+ntmp);
    double mx2pool   = (nk(i)*mx2k(i)+ntmp*mx2tmp)/(nk(i)+ntmp);
    double my2pool   = (nk(i)*my2k(i)+ntmp*my2tmp)/(nk(i)+ntmp);
    double mxypool   = (nk(i)*mxyk(i)+ntmp*mxytmp)/(nk(i)+ntmp);
    double cxypool   = mxypool-mxpool*mypool;
    double vxpool    = mx2pool-mxpool*mxpool;
    double vypool    = my2pool-mypool*mypool;
    double varpool   = vypool- cxypool*cxypool/vxpool;
    double vark      = vyk(i)- cxyk(i)*cxyk(i)/vxk(i);
    double varr      = vytmp- cxytmp*cxytmp/vxtmp;
    Dtmp(i)          = (nk(i)+ntmp)*varpool -nk(i)*vark-ntmp*varr;
  }

  void compactEM_init_Reg::CAHminReg(int k, int& imin, int& jmin, double& dmin){
    int nbD = ((k-1)*k)/2; // 1 -> _K-1 elements 
    int iimin;    
    int iimin_shared;
    double dmin_p = std::numeric_limits<double>::max();
    double dmin_shared = std::numeric_limits<double>::max();

#pragma omp parallel if (nbD>1000) shared(dmin_shared, iimin_shared) private(iimin) firstprivate(dmin_p)
    {
#pragma omp for
      for(int ii=0; ii<nbD; ++ii)
	{
	  double distij = *(_D[0]+ii);
	  if(distij<=dmin_p){
	    dmin_p = distij;
	    iimin = ii;
	  }
	}
#pragma omp critical 
      {
	if(dmin_p<=dmin_shared){
	  dmin_shared = dmin_p;
	  iimin_shared = iimin;	
	}
      }
    } 
    dmin = dmin_shared;
    //std::cerr<<dmin<<" "<<iimin_shared<<std::endl;
    imin = int(floor(-0.5+sqrt(0.25+2*iimin_shared)));
    jmin = iimin_shared-(imin+1)*imin/2;
    imin += 2; jmin += 1;
    //std::cerr<<dmin<<" "<<imin<<" "<<jmin<<std::endl<<std::endl;    
  }

  void compactEM_init_Reg::CAHcopyReg(int k, int imin, int jmin){
    for (int i=1; i<(jmin-1+1);i++)
      Dist(jmin,i) = Dtmp(i);
    for (int i=jmin+1;i<(k-1+1);i++)
      Dist(i,jmin) = Dtmp(i);
    for (int j=1;j<(jmin-1+1);j++)
      Dist(imin,j) = Dist(k,j);    
    for (int j=jmin+1;j<(imin-1+1);j++)
      Dist(imin,j) = Dist(k,j);
    for (int i=imin+1;i<(k-1+1);i++)
      Dist(i,imin) = Dist(k,i);
  }

}
#endif



