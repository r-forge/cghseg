//g++ numlib.cc -c && g++ compactEM_init.cc -c && g++ main_CAH.cpp compactEM_init.o numlib.o -o main_CAH
//g++ -O3 numlib.cc -c && g++ -O3 compactEM_init.cc -c && g++ -O3 main_CAH.cpp compactEM_init.o numlib.o -o main_CAH
//g++ -g -pg numlib.cc -c && g++ -g -pg compactEM_init.cc -c && g++ -g -pg main_CAH.cpp compactEM_init.o numlib.o -o main_CAH
//g++ -O3 -fopenmp numlib.cc -c && g++ -O3 -fopenmp compactEM_init.cc -c && g++ -O3 -fopenmp main_CAH.cpp compactEM_init.o numlib.o -o main_CAH

#include <iostream>
#include <fstream>
#include <cstdio>
#include <time.h>
#include "compactEM_init.h"

using namespace cghseg;
using namespace std;

int main(int argc, char *argv[])
{
  ifstream ifile("sc_EMinit_input.txt");
  long K;
  ifile>>K;
  long P;
  ifile>>P;
  bool vh     = true;
  int vhInt;
  ifile>>vhInt;
  if (!vhInt)
    vh=false;
  
  double *datak = new double[K];
  int lengthxkR;
  ifile>>lengthxkR;
  for (int i=0;i<lengthxkR;i++){
    ifile>>datak[i];
  }
  
  double *data2k = new double[K];
  int lengthx2kR;
  ifile>>lengthx2kR;
  for (int i=0;i<lengthx2kR;i++){
    ifile>>data2k[i];
  }

  double *datank = new double[K];
  int lengthnkR;
  ifile>>lengthnkR;
  for (int i=0;i<lengthnkR;i++){
    ifile>>datank[i];
  }
  ifile.close();
  
  compactEM_init compactEMi(K,P);  
  compactEMi.Init(datak,data2k,datank);

  time_t start,end;
  time (&start);
  compactEMi.CAH();
  time (&end);
  double dif = difftime (end,start);
  printf ("It took %.2lf seconds to do CAH.\n", dif);
}
