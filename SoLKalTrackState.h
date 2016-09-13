#ifndef ROOT_SOL_KAL_TRACK_STATE
#define ROOT_SOL_KAL_TRACK_STATE
//SoLIDTracking
#include "SoLKalMatrix.h"
#include "SoLIDUtility.h"
#include "SoLKalTrackSite.h"
class SoLKalTrackSite;
class SoLKalFieldStepper;

class SoLKalTrackState : public SoLKalMatrix {
  public:
  
  SoLKalTrackState(Int_t type = 0, Int_t p = kSdim);
  SoLKalTrackState(const SoLKalMatrix &sv, Int_t type = 0, Int_t p = kSdim);
  SoLKalTrackState(const SoLKalMatrix &sv, const SoLKalMatrix &c,
              Int_t type = 0, Int_t p = kSdim);
  SoLKalTrackState(const SoLKalMatrix &sv, const SoLKalTrackSite &site,
              Int_t type = 0, Int_t p = kSdim);
  SoLKalTrackState(const SoLKalMatrix &sv, const SoLKalMatrix &c,
              const SoLKalTrackSite &site, Int_t type = 0, Int_t p = kSdim);
  ~SoLKalTrackState();

  virtual SoLKalTrackState * MoveTo(SoLKalTrackSite  &to,
                               SoLKalMatrix &F,
                               SoLKalMatrix *QPtr = 0) const;
  virtual SoLKalTrackState & MoveTo(SoLKalTrackSite  &to,
                              SoLKalMatrix &F,
                              SoLKalMatrix &Q) const;
  virtual SoLKalTrackState * MoveToZ(Double_t z, SoLKalMatrix &F,
                                     SoLKalMatrix *QPtr = 0) const;
  virtual SoLKalTrackState & MoveToZ(Double_t z, SoLKalMatrix &F,
                                     SoLKalMatrix &Q) const;
  virtual void         Propagate(SoLKalTrackSite &to);
  virtual SoLKalTrackState * PredictSVatZ(Double_t &z);
  virtual SoLKalTrackState * PredictSVatNextZ(Double_t &z);
  virtual void InitPredictSV();

  inline void  ClearAttemptSV() { fAttemptState = nullptr; }
  inline Int_t GetDimension                () const { return GetNrows(); }
  inline const SoLKalTrackSite  & GetSite        () const { return *fSitePtr; }
  inline const SoLKalMatrix & GetCovMat      () const { return fC; }
  inline const SoLKalMatrix & GetProcNoiseMat() const { return fQ; }
  inline const SoLKalMatrix & GetPropMat     (const Char_t *t = "") const { return (t[0] == 'T' ? fFt : fF); }
  inline const Int_t & GetType() const { return fType; }
  inline const Double_t & GetZ0() const { return fZ0; }

  inline void SetStateVec    (const SoLKalMatrix &c) { TMatrixD::operator=(c); }
  inline void SetCovMat      (const SoLKalMatrix &c) { fC       = c; }
  inline void SetProcNoiseMat(const SoLKalMatrix &q) { fQ       = q; }
  inline void SetSitePtr     (SoLKalTrackSite  *s)   { fSitePtr = s; }
  inline void SetZ0          (Double_t z)            { fZ0 = z; }
  
  virtual void   CalcDir (TVector3 &dir) const;
  static  void   CalcDir (TVector3 &dir, const SoLKalMatrix &sv);
  virtual void   CalcMomVec (TVector3 &dir) const;
  static  void   CalcMomVec (TVector3 &dir, const SoLKalMatrix &sv);
  
  private:

   // private data members -------------------------------------------

   Int_t            fType;    // (0,1,2,3) = (uninited,predicted,filtered,smoothed)
   SoLKalTrackSite  *fSitePtr; // pointer to corresponding KalSite
   SoLKalMatrix     fF;       // propagator matrix to next site (F = @f/@a)
   SoLKalMatrix     fFt;      // transposed propagator matrix (F^T = (@f/@a)^T)
   SoLKalMatrix     fQ;       // process noise from this to the next sites
   SoLKalMatrix     fC;       // covariance matrix
   SoLKalTrackState* fAttemptState; //a pointer to the next predicted state vector
                                    //for finding hits on the next detector
                                    
   Double_t         fZ0;
   
   //SoLKalFieldStepper* fFieldStepper;
   ClassDef(SoLKalTrackState,1)      // Base class for state vector objects

};


#endif
