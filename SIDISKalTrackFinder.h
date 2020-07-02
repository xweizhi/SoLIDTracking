#ifndef ROOT_SIDIS_KAL_TRACK_FINDER
#define ROOT_SIDIS_KAL_TRACK_FINDER
//c++
#include <iostream>
#include <map>
#include <cassert>
#include <vector>
//ROOT
#include "TClonesArray.h"
#include "TSeqCollection.h"
#include "TMath.h"
#include "TVector3.h"
#include "TVector2.h"

//Hall A analyzer
#include "THaAnalysisObject.h"

//SoLIDTracking
#include "SoLIDTrack.h"
#include "SoLIDGEMHit.h"
#include "SoLIDUtility.h"
#include "SoLIDFieldMap.h"
#include "SoLIDGEMTracker.h"
#include "SoLIDGEMChamber.h"
#include "SoLIDECal.h"
#include "SoLKalMatrix.h"
#include "SoLKalFieldStepper.h"
#include "SoLKalTrackFinder.h"


using namespace std;

class SIDISKalTrackFinder : public SoLKalTrackFinder, public THaAnalysisObject
{
  public:
  SIDISKalTrackFinder() {;}
  SIDISKalTrackFinder(bool isMC, const char* name = "SIDISTrackFinder", int detConf = 0);
  ~SIDISKalTrackFinder();
  
  void ProcessHits(TClonesArray* theTracks);
                           
  void Clear( Option_t* opt="" );
  Int_t ReadDatabase (const TDatime& date);
#ifdef MCDATA
  bool GetSeedEfficiency(int i) const { return fSeedEfficiency[i]; }
  bool GetMCTrackEfficiency(int i) const { return fMcTrackEfficiency[i]; }
#endif
  
  
  protected:
  
  //Main analysis functions
  void FindDoubletSeed(Int_t planej, Int_t planek, ECType type = kFAEC);
  void MergeSeed();
  void TrackFollow();
  void CoarseCheckVertex();
  void FindandAddVertex();
  void FinalSelection(TClonesArray* theTracks);
  void CopyTrack(SoLIDTrack* soltrack, SoLKalTrackSystem* kaltrack);
  void ECalFinalMatch();
  
  
  //assistent functions
  double CalDeltaPhi(const double & phi1, const double & phi2); 
  double CalDeltaR(const double & r1, const double & r2);
  SoLKalTrackSite & SiteInitWithSeed(DoubletSeed* thisSeed);
  Double_t TriggerCheck(SoLIDGEMHit* theHit, ECType type);
  SoLIDGEMHit* FindCloestHitInWindow(double &x, double &y);
  double PredictR(Int_t &plane, SoLIDGEMHit* hit1, SoLIDGEMHit* hit2);
  int GetHitsInWindow(int plane, double x, double wx, double y, double wy, bool flag = false);
  Double_t FindVertexZ(SoLKalTrackState* thisState);
  Bool_t CheckChargeAsy(SoLKalTrackSystem* theSystem);
  void GetHitChamberList(vector<Int_t> &theList, Int_t thisChamber, Int_t size);
  Int_t GetChamIDFromPos(Double_t &x, Double_t &y, Int_t TrackerID);
  Bool_t CalInitParForPair(SoLIDGEMHit* hita, SoLIDGEMHit* hitb, Double_t &charge, 
                           Double_t& mom, Double_t& theta, Double_t& phi, ECType& type);
  Int_t BinarySearchForR(TSeqCollection* array, Double_t &lowr);
    
  bool fIsMC;
  bool fSeedEfficiency[2];
  bool fMcTrackEfficiency[2];
  vector<SoLIDGEMHit*> fWindowHits;
  map< Int_t, vector<SoLIDGEMHit*> > fGoodHits;
  Int_t fNGoodTrack;
  Int_t fDetConf;
  
  Double_t fMomMinCut[2];
  Double_t fMomMaxCut[2];
  Double_t fThetaMinCut[2];
  Double_t fThetaMaxCut[2];
  Double_t fCellEdgeCut[2];
  Double_t fCoarseCellEdgeCut[2];
  Double_t fRlimitMin[3][2];
  Double_t fRlimitMax[3][2];
  Double_t fDeltaRMin[5][2];
  Double_t fDeltaRMax[5][2];
  Double_t fDeltaPhiMin[5][2];
  Double_t fDeltaPhiMax[5][2];
  Double_t fCoarseECPosCut[2];
  Double_t fECPosCut;
  Int_t    fECEnergyMatch[2];
  
  
  ClassDef(SIDISKalTrackFinder,0)
};

#endif

