#include "compactEM_init_Reg.h"
#include <fstream>
#include <list>
#include <cmath>
#include <math.h>
#include <cstdlib>
#include <memory>
#include <vector>
#include <cstring>

#include <Rconfig.h>
#ifdef SUPPORT_OPENMP
  #include <omp.h>
#endif

//#define  Dist(i,j) _D[i-2][j-1]
//#define  nk(i)     _nk[i-1]
//#define  mk(i)     _mk[i-1]
//#define  vk(i)     _vk[i-1]
//#define  Dtmp(i)   _Dtmp[i-1]


using namespace std;

namespace cghseg
{

  void compactEM_init_Reg::compute_phi_Reg(){
    for (int p=0;p<_P;p++){
      _phi[p]      = _myk[p]-_mxk[p]*_cxyk[p]/_vxk[p];
      _phi[p+_P]   = _cxyk[p]/_vxk[p];
      _phi[p+2*_P] =  _vyk[p]- _cxyk[p]*_cxyk[p]/_vxk[p];      
      _phi[p+3*_P] = double(_nk[p])/_lengthx;
    }
  }

  void compactEM_init_Reg::CAH_Reg(){
    for (int k=_K;k>=_P+1;k--){
      int imin    = 2; 
      int jmin    = 1; 
      double dmin = Dist(2,1);
      CAHminReg(k, imin, jmin, dmin);

      int    ntmp     =  nk(imin)+nk(jmin);    
      double mxtmp    = (nk(imin)*mxk(imin)+nk(jmin)*mxk(jmin))/ntmp;
      double mytmp    = (nk(imin)*myk(imin)+nk(jmin)*myk(jmin))/ntmp;
      double mx2tmp   = (nk(imin)*mx2k(imin)+nk(jmin)*mx2k(jmin))/ntmp;
      double my2tmp   = (nk(imin)*my2k(imin)+nk(jmin)*my2k(jmin))/ntmp;
      double mxytmp   = (nk(imin)*mxyk(imin)+nk(jmin)*mxyk(jmin))/ntmp;
      double cxytmp   = mxytmp-mxtmp*mytmp;
      double vxtmp    = mx2tmp-mxtmp*mxtmp;
      double vytmp    = my2tmp-mytmp*mytmp;
      double vartmp   = vytmp- cxytmp*cxytmp/vxtmp;
      double vark     = vyk(imin)- cxyk(imin)*cxyk(imin)/vxk(imin);
      double varr     = vyk(jmin)- cxyk(jmin)*cxyk(jmin)/vxk(jmin);
      mxk(jmin)  = mxtmp;
      mx2k(jmin) = mx2tmp;
      vxk(jmin)  = vxtmp;
      myk(jmin)  = mytmp;
      my2k(jmin) = my2tmp;
      vyk(jmin)  = vytmp;
      mxyk(jmin) = mxytmp;
      cxyk(jmin) = cxytmp;
      nk(jmin)   = ntmp;

      int i1 = 0;
      int i2 = 0;
      int i3 = 0;

      vector<int> index;

      if (jmin>1){ // 1:(jmin-1)
	i1     = jmin-1;     
	for (int h = 1; h<i1+1 ;h++){
	  CAHcoreReg(h,ntmp,mxtmp,mx2tmp,vxtmp,mytmp,my2tmp,vytmp,mxytmp,cxytmp);
	}
      } else {
	i1 = 1;
	CAHcoreReg(1,ntmp,mxtmp,mx2tmp,vxtmp,mytmp,my2tmp,vytmp,mxytmp,cxytmp);
      }
      if (jmin+1<=imin-1){ // (jmin+1):(imin-1)
	i2     = imin-1-jmin-1+1;
	for (int h=1;h<(i2+1);h++)
	  CAHcoreReg(h+jmin,ntmp,mxtmp,mx2tmp,vxtmp,mytmp,my2tmp,vytmp,mxytmp,cxytmp);
      } 
      if (imin+1<=(k)){// (imin+1):k
	i3     = k-imin-1+1;
	for (int h=1;h<(i3+1);h++)
	  CAHcoreReg(h+imin,ntmp,mxtmp,mx2tmp,vxtmp,mytmp,my2tmp,vytmp,mxytmp,cxytmp);
      }  
      
      double auxmx  = mxk(imin);  mxk(imin) = mxk(k);   mxk(k) = auxmx;
      double auxvx  = vxk(imin);  vxk(imin) = vxk(k);   vxk(k) = auxvx;
      double auxmy  = myk(imin);  myk(imin) = myk(k);   myk(k) = auxmy;
      double auxvy  = vyk(imin);  vyk(imin) = vyk(k);   vyk(k) = auxvy;
      double auxmxy = mxyk(imin); mxyk(imin) = mxyk(k); mxyk(k) = auxmxy;
      double auxcxy = cxyk(imin); cxyk(imin) = cxyk(k); cxyk(k) = auxcxy;
      int    auxn   = nk(imin);   nk(imin) = nk(k);     nk(k) = auxn;
      double auxk   = Dtmp(k);

      Dtmp(k)     = Dtmp(imin);
      Dtmp(imin)  = auxk; 
      CAHcopyReg(k, imin, jmin);
    
    } // end k
  } //end CAH



  void compactEM_init_Reg::Init(double *Dataxk, double *Datax2k,double *Datayk, double *Datay2k,double *Dataxyk,double *Datank){

    memcpy(_mxk, Dataxk, sizeof(double)*_K);
    memcpy(_mx2k,Datax2k,sizeof(double)*_K);
    memcpy(_myk, Datayk, sizeof(double)*_K);
    memcpy(_my2k,Datay2k,sizeof(double)*_K);
    memcpy(_mxyk,Dataxyk,sizeof(double)*_K);
    memcpy(_nk,  Datank, sizeof(double)*_K);
  
    for (int k = 0 ;  k < _K; k++){
	_vxk[k]   = _mx2k[k]-_mxk[k]*_mxk[k];
	_vyk[k]   = _my2k[k]-_myk[k]*_myk[k];
	_cxyk[k]  = _mxyk[k]-_myk[k]*_mxk[k];
	_lengthx += _nk[k];
    }

    for (int k = 0 ;  k < _K; k++){
      _mxk0[k]  = _mxk[k];
      _vxk0[k]  = _vxk[k];
      _myk0[k]  = _myk[k];
      _vyk0[k]  = _vyk[k];
      _mxyk0[k] = _mxyk[k];
      _cxyk0[k] = _cxyk[k];
      _nk0[k]   = _nk[k];
    }
  
    for (int k = 0 ;  k < _K-1; k++){
      for (int r = 0 ;  r < k+1; r++){
	double mxpool    = (_nk[k+1]*_mxk[k+1]+_nk[r]*_mxk[r])/(_nk[k+1]+_nk[r]);
	double mypool    = (_nk[k+1]*_myk[k+1]+_nk[r]*_myk[r])/(_nk[k+1]+_nk[r]);
	double mx2pool   = (_nk[k+1]*_mx2k[k+1]+_nk[r]*_mx2k[r])/(_nk[k+1]+_nk[r]);
	double my2pool   = (_nk[k+1]*_my2k[k+1]+_nk[r]*_my2k[r])/(_nk[k+1]+_nk[r]);
	double mxypool   = (_nk[k+1]*_mxyk[k+1]+_nk[r]*_mxyk[r])/(_nk[k+1]+_nk[r]);
	double cxypool   = mxypool-mxpool*mypool;
	double vxpool    = mx2pool-mxpool*mxpool;
	double vypool    = my2pool-mypool*mypool;
	double varpool   = vypool- cxypool*cxypool/vxpool;
	double vark      = _vyk[k+1]- _cxyk[k+1]*_cxyk[k+1]/_vxk[k+1];
	double varr      = _vyk[r]  - _cxyk[r]*_cxyk[r]/_vxk[r];
	_D[k][r]         = (_nk[k+1]+_nk[r])*varpool -_nk[k+1]*vark-_nk[r]*varr;
      }
    }
  } //end Init


  compactEM_init_Reg::compactEM_init_Reg(int nbsegments, int nbclusters, int OMP_NUM_THREADS){
#ifdef SUPPORT_OPENMP
    omp_set_num_threads(OMP_NUM_THREADS);
#endif

    _K       = nbsegments; 
    _P       = nbclusters;
    _D       = new double *[_K-1];
    _Dtmp    = new double [_K];
    _phi     = new double [6*_P];
    _tau     = new double *[_K];
    _mxk     = new double[_K];
    _mx2k    = new double[_K];
    _myk     = new double[_K];
    _my2k    = new double[_K];
    _mxyk    = new double[_K];
    _cxyk    = new double[_K];
    _vxk     = new double[_K];
    _vyk     = new double[_K];
    _nk      = new double[_K];
    _mxk0    = new double[_K];
    _mx2k0   = new double[_K];
    _vxk0    = new double[_K];
    _myk0    = new double[_K];
    _my2k0   = new double[_K];
    _vyk0    = new double[_K];
    _mxyk0   = new double[_K];
    _cxyk0   = new double[_K];
    _nk0     = new int[_K];
  
    for (int p = 0; p < 6*_P; p++)
      _phi[p] = 0;
  
    _tau[0] = new double[_K*_P];
    for (int k =1; k<_K; k++)
      _tau[k] =  _tau[k-1] + _P;    
    for (int k =0; k<_K; k++){
      for (int p =0; p<_P; p++){
	_tau[k][p] = 0;
      }
    }

    int nbD = ((_K-1)*_K)/2; // 1 -> _K-1 elements
    _D[0] = new double[nbD];
    for (int k =1; k<_K-1; k++){
      _D[k] = _D[k-1]+ k;
    }
    // for (int k =0; k<_K-1; k++){
    //   for (int r=0; r<k+1; r++){
    // 	_D[k][r] = 0;
    //   }
    // }

#pragma omp parallel if (nbD>1000)
    { // page placement by first touch
#pragma omp for
      for(int ii=0; ii<nbD; ++ii)
	{
	  *(_D[0]+ii) = 0.;
	}
    }
  
    for (int k =0; k<_K; k++){
      _Dtmp[k]  = 0;
      _mxk[k]   = 0;
      _mx2k[k]  = 0;
      _myk[k]   = 0;
      _my2k[k]  = 0;
      _mxyk[k]  = 0;
      _cxyk[k]  = 0;
      _vxk[k]   = 0;
      _vyk[k]   = 0;
      _nk[k]    = 0;
      _mxk0[k]  = 0;
      _myk0[k]  = 0;
      _mx2k0[k] = 0;
      _my2k0[k] = 0;
      _mxyk0[k] = 0;
      _cxyk0[k] = 0;
      _vxk0[k]  = 0;
      _vyk0[k]  = 0;
      _nk0[k]   = 0;
    }
  
  } // end constructor


  // destructor

  compactEM_init_Reg::~compactEM_init_Reg(){
  
    delete[] _phi;
    delete[] _mxk;
    delete[] _mx2k;
    delete[] _myk;
    delete[] _my2k;
    delete[] _mxyk;
    delete[] _cxyk;
    delete[] _vxk;
    delete[] _vyk;
    delete[] _nk;
    delete[] _mxk0;
    delete[] _mx2k0;
    delete[] _myk0;
    delete[] _my2k0;
    delete[] _mxyk0;
    delete[] _cxyk0;
    delete[] _vxk0;
    delete[] _vyk0;
    delete[] _nk0;
    delete[] _Dtmp;
  
    delete[] _tau[0];
    delete[] _tau;
    delete[] _D[0];
    delete[] _D; 
  } // end destructor

  void   
  compute_compactEM_init_Reg(numlib_vector *xk,numlib_vector *x2k, numlib_vector *yk,numlib_vector *y2k, numlib_vector *xyk, numlib_vector *nk, int P , numlib_vector **phiend, int OMP_NUM_THREADS)
  {
        
    double *dataxk=new double[xk->size];
    for (int i=0;i<xk->size;i++)
      dataxk[i]=numlib_vector_get(xk,i);

    double *datax2k=new double[x2k->size];
    for (int i=0;i<x2k->size;i++)
      datax2k[i]=numlib_vector_get(x2k,i);

    double *datayk=new double[yk->size];
    for (int i=0;i<yk->size;i++)
      datayk[i]=numlib_vector_get(yk,i);

    double *datay2k=new double[y2k->size];
    for (int i=0;i<y2k->size;i++)
      datay2k[i]=numlib_vector_get(y2k,i);
    
    double *dataxyk=new double[xyk->size];
    for (int i=0;i<xyk->size;i++)
      dataxyk[i]=numlib_vector_get(xyk,i);

    double *datank=new double[nk->size];
    for (int i=0;i<nk->size;i++)
      datank[i]=numlib_vector_get(nk,i);

    int K   = nk->size;

    compactEM_init_Reg compactEMiReg(K,P,OMP_NUM_THREADS);
    compactEMiReg.Init(dataxk,datax2k,datayk,datay2k,dataxyk,datank);
    compactEMiReg.CAH_Reg();
    compactEMiReg.compute_phi_Reg();
           
    (*phiend)=numlib_vector_calloc(6*P);
    for (int p=0;p<6*P;p++)
      numlib_vector_set((*phiend),p,compactEMiReg._phi[p]);
    
    delete[] dataxk;
    delete[] datax2k;
    delete[] datayk;
    delete[] datay2k;
    delete[] dataxyk;
    delete[] datank;

  }  // end compute_EM_init

} // end namespace cghseg






