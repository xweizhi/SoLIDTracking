#ifndef ROOT_SoLID_Progressive_Tracking
#define ROOT_SoLID_Progressive_Tracking
//c++
#include <map>
#include <cassert>
#include <vector>
#include <iostream>
//Hall A Analyzer
//#include "VarDef.h"
//ROOT
#include "TClonesArray.h"
#include "TVector3.h"
//SoLIDTracking
#include "SoLIDTrack.h"
#include "SoLIDGEMHit.h"
#include "SoLIDUtility.h"

using namespace std;

class ProgressiveTracking
{
  public:
  ProgressiveTracking(Int_t ntracker, Bool_t isMC);
  ~ProgressiveTracking();
  
  void               ProcessHits(map<Int_t, vector<TSeqCollection*> > *theHitMap, 
                           TClonesArray* theTracks);
  Int_t              ReadDataBase();
  void               Clear( Option_t* opt="" );
  
  Bool_t             IsIterBackward() const { return fIsIterBackward; }
  Bool_t             HasCaloHit()     const { return fHasCaloHit; }
  
  void               SetCaloHit(Double_t xpos, Double_t ypos, Int_t plane, Double_t edp){
                      fCaloHits.push_back(SoLIDCaloHit(xpos, ypos, plane, edp)); 
                     }
  void               SetIterBackward(Bool_t is) { fIsIterBackward = is; }
  void               SetNElectron(Int_t n) { fNElectron = n; }
  void               SetNHadron(Int_t n) { fNHadron = n; }
  
  private:
  void               FindTrack(Int_t angleflag, Int_t type, map<Int_t, vector<TSeqCollection*> > *theHitMap);
  void               CheckTracks();
  void               CombineTrackRoad(TClonesArray* theTracks);
  Int_t              FindMomRange(Int_t layer1, Double_t phi0, Int_t layer2, Double_t phi1, 
                                  Double_t* mom_min, Double_t* mom_max,Int_t angleflag);
  Int_t              FindThetaRange(Double_t r1, Int_t layer1, Double_t r2, Int_t layer2,
                                    Double_t* thetaMin, Double_t* thetaMax, Int_t angleflag);
  void               CalCircle(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Double_t x3, 
                               Double_t y3, Double_t* R, Double_t* Xc, Double_t* Yc);
  Double_t           CalDeltaPhi(Double_t phi1, Double_t phi2);                             
  Double_t           CalChargeAsy(Double_t qu, Double_t qv);
  Bool_t             CalVertexZ(SoLIDGEMHit* hit1, SoLIDGEMHit* hit2, SoLIDGEMHit *hit3, Double_t *Rmom, Int_t angleFlag,
                                Double_t *reconVertexZ, Double_t *reconTheta);
  
  Bool_t             fDoMC;
  Int_t              fNTracker;
  Int_t              fNElectron;
  Int_t              fNHadron;
  Int_t              fNTrack;
  Bool_t             fIsIterBackward;
  Bool_t             fHasCaloHit;
  TClonesArray*      fCoarseTracks; //internal track array of progressive tracking 
  Int_t              fNGoodHits[6];
  TClonesArray*      fGoodHits[6];
  
  struct SoLIDCaloHit{
  Double_t fXPos;
  Double_t fYPos;
  Int_t    fECID;
  Double_t fEdp;
  SoLIDCaloHit() {}
  SoLIDCaloHit(Double_t xpos, Double_t ypos, Int_t plane, Double_t edp) :
  fXPos(xpos), fYPos(ypos), fECID(plane), fEdp(edp) {}
  };
  vector<SoLIDCaloHit> fCaloHits;
  
};

#endif
