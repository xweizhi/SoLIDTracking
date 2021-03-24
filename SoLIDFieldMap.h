#ifndef ROOT_SoLID_Field_Map
#define ROOT_SoLID_Field_Map
//c++
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
//ROOT
#include "TVector3.h"
#include "TRandom3.h"
//SoLIDTracking
#include "SoLIDUtility.h"

#define ZSIZE 1201
#define RSIZE 501
#define ZSHIFT 600
#define ZSTEP 1
#define RSTEP 1
#define ZSIZE_TARGET 401
#define RSIZE_TARGET 200
#define ZSTEP_TARGET 1
#define RSTEP_TARGET 1
#define ZSHIFT_TARGET 200
#define FIELD_PRECISION 0.0

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
  void LoadTargetFieldMap();  

  protected:
  SoLIDFieldMap();
  static SoLIDFieldMap *fInstance;
  void LoadFieldMap();
  
  Double_t  Bz[ZSIZE][RSIZE]; //hard coded for now, array to store SoLID B field
  Double_t  Br[ZSIZE][RSIZE]; //hard coded for now, array to store SoLID B field
  Double_t  Bz_target[ZSIZE_TARGET][RSIZE_TARGET];
  Double_t  Br_target[ZSIZE_TARGET][RSIZE_TARGET];
  TVector3  fField;
  bool      targetFieldFlag;
  TRandom3  fFieldRNG;   
};  

#endif
