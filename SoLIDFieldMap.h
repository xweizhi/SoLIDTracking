#ifndef ROOT_SoLID_Field_Map
#define ROOT_SoLID_Field_Map
//c++
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
//ROOT
#include "TVector3.h"
//SoLIDTracking
#include "SoLIDUtility.h"

#define ZSIZE 1201
#define RSIZE 501
#define ZSHIFT 600
#define ZSTEP 1
#define RSTEP 1

using namespace std;

class SoLIDFieldMap
{
  public: 
  ~SoLIDFieldMap() {;}
  static SoLIDFieldMap * GetInstance() {
    if (fInstance == NULL) fInstance = new SoLIDFieldMap();
    return fInstance;
  }
  
  TVector3 & GetBField(double x, double y, double z);
  
  protected:
  SoLIDFieldMap();
  static SoLIDFieldMap *fInstance;
  void LoadFieldMap();
  
  Double_t  Bz[ZSIZE][RSIZE]; //hard coded for now, array to store SoLID B field
  Double_t  Br[ZSIZE][RSIZE]; //hard coded for now, array to store SoLID B field
  TVector3  fField;
  
};  

#endif
