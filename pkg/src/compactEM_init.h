#ifndef CGHSEG_COMPACTEMINIT_H
#define CGHSEG_COMPACTEMINIT_H

#include "numlib.h"

#include <iostream>
#include <fstream>

namespace cghseg
{
  void compute_compactEM_init(numlib_vector *xk,numlib_vector *x2k, numlib_vector *nk, int P , numlib_vector **phiend);
       
  class compactEM_init{
    
  public:
    int       _lengthx;        // size of the data
    int       _K;              // Number of segments 
    int       _P;              // Number of clusters
    double   *_phi;            // parameters
    double   *_xk;              // data
    double   *_mk;             // empirical means
    double   *_m2k;             // empirical means
    double   *_vk;             // empirical var
    double   *_nk;             // length of segments
    double   *_mk0;             // empirical means
    double   *_vk0;             // empirical var
    int      *_nk0;             // length of segments
    double   **_D;             // distance matrix
    double   *_Dtmp;           // distance vector of the merged groups
    double   **_tau;           // posterior
    
    compactEM_init(int nbseg, int nbclust);
    void CAH();
    void compute_phi();
    void Init(double *Datak, double *Data2k,double *Datank);    
    friend std::ostream & operator << (std::ostream &s, const compactEM_init & compactEMi);
    ~compactEM_init();
  };
}
#endif



