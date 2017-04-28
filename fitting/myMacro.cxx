#include <iostream>
#include "TROOT.h"

using namespace std;

int myMacro(int k=0){

  gROOT->ProcessLine(".!ls"); // on Linux

  cout << "The input parameter was " << k << endl;
  return k;

}
