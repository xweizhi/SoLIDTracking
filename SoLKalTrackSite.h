#ifndef ROOT_SOL_KAL_TRACK_SITE
#define ROOT_SOL_KAL_TRACK_SITE
//ROOT
#include "TObjArray.h"
#include "TMath.h"
//SoLIDTracking
#include "SoLKalMatrix.h"
#include "SoLKalTrackState.h"
#include "SoLIDUtility.h"
#include "SoLIDGEMHit.h"
class SoLKalTrackState;

class SoLKalTrackSite : public TObjArray {
  friend class SoLKalTrackSystem;
  
  public:
   enum EStType { kPredicted = 0,
                  kFiltered,
                  kSmoothed,
                  kInvFiltered };
                  
  SoLKalTrackSite(Int_t m = kMdim, Int_t p = kSdim, Double_t chi2 = 60.);
  SoLKalTrackSite(SoLIDGEMHit* ht, Int_t m = kMdim, Int_t p = kSdim , Double_t chi2 = 60.);
  ~SoLKalTrackSite();
  
  Int_t   CalcExpectedMeasVec  (const SoLKalTrackState &a,
                                      SoLKalMatrix &m);
  Int_t   CalcMeasVecDerivative(const SoLKalTrackState &a,
                                      SoLKalMatrix &H);
  Bool_t  Filter();
  void    Smooth(SoLKalTrackSite &pre);
  void    InvFilter();
  void    Add(TObject *obj);
  Bool_t  IsAccepted();
  void    SetHitResolution(Double_t ex, Double_t ey);
  void    SetMeasurement(Double_t x, Double_t y);       

  inline Int_t        GetDimension() const { return fM.GetNrows(); }
  inline SoLKalTrackState & GetCurState ()       { return *fCurStatePtr; }
  inline SoLKalTrackState & GetCurState () const { return *fCurStatePtr; }
  SoLKalTrackState & GetState (EStType t);
  inline SoLKalMatrix & GetMeasVec      ()   { return fM;            }
  inline SoLKalMatrix & GetMeasNoiseMat ()   { return fV;            }
  inline SoLKalMatrix & GetResVec       ()   { return fResVec;       }
  inline SoLKalMatrix & GetCovMat       ()   { return fR;            }
  inline Double_t     GetDeltaChi2() const { return fDeltaChi2;    }
  inline Double_t     GetZ() const { return fZ0; }
  inline void SetZ(Double_t z) { fZ0 = z; }
  SoLKalMatrix   GetResVec (EStType t);
         
  inline const SoLIDGEMHit * GetHit    () const { return fGEMHit; }    
  SoLIDGEMHit * GetPredInfoHit();   
  private:
   // Private utility methods

   SoLKalTrackState & CreateState(const SoLKalMatrix &sv, Int_t type = 0);
   SoLKalTrackState & CreateState(const SoLKalMatrix &sv, const SoLKalMatrix &c,
                                    Int_t type = 0);
   inline  Int_t CalcXexp(const SoLKalTrackState &a, TVector3   &xx) const;
  private:

   // private data member ------------------------------------------- 
   SoLKalTrackState*  fCurStatePtr; // pointer to current best state
   SoLKalMatrix       fM;           // measurement vector: M(m,1)
   SoLKalMatrix       fV;           // noise matrix: M(m,m)
   SoLKalMatrix       fH;           // H = (@h/@a): M(m,p)
   SoLKalMatrix       fHt;          // H^t = (@h/@a)^t: M(p,m)
   SoLKalMatrix       fResVec;      // m - h(a): M(m,1)
   SoLKalMatrix       fR;           // covariance matrix: M(m,m)
   Double_t           fDeltaChi2;   // chi2 increment
   Double_t           fMaxDeltaChi2;
   SoLIDGEMHit* fGEMHit;
   Double_t           fZ0;
   ClassDef(SoLKalTrackSite,1)      // Base class for measurement vector objects

};


#endif
