#ifndef ROOT_SOL_KAL_TRACK_SYSTEM
#define ROOT_SOL_KAL_TRACK_SYSTEM
//*************************************************************************
//* ===================
//*  TVKalSystem Class
//* ===================
//*
//* (Description)
//*   Base class for Kalman filtering class
//* (Requires)
//* 	TObjArray
//* (Provides)
//* 	class TVKalSystem
//* (Update Recored)
//*   2003/09/30  K.Fujii	Original version.
//*   2005/08/25  A.Yamaguchi	Added fgCurInstancePtr and its getter & setter.
//*   2009/06/18  K.Fujii       Implement inverse Kalman filter
//*   2016/07/20  W.Xiong   modified for SoLID
//*************************************************************************
//ROOT
#include "TObjArray.h"
#include "TMath.h"
//SoLIDTracking
#include "SoLKalTrackSite.h"
#include "SoLKalTrackState.h"
#include "SoLIDUtility.h"
class SoLKalMatrix;
class SoLKalTrackState;

class SoLKalTrackSystem : public TObjArray {
  friend class SoLKalTrackSite;
  
  public:
  SoLKalTrackSystem(Int_t n = 1);
  ~SoLKalTrackSystem();
  
  Bool_t AddAndFilter(SoLKalTrackSite &next);
  void   SmoothBackTo(Int_t k);
  void   SmoothAll();
  void   InvFilter(Int_t k);

  void   Add(TObject *obj);
  inline SoLKalTrackSite  & GetCurSite() { return *fCurSitePtr; }
  inline SoLKalTrackState & GetState(SoLKalTrackSite::EStType t)
                            { return fCurSitePtr->GetState(t); }
  void   UpdateChi2();
  Double_t GetChi2() const {return fChi2; };
  void IncreaseChi2(Double_t dchi2) { fChi2 += dchi2; } 
  Int_t    GetNDF (Bool_t self = kTRUE) const;
  Double_t GetChi2perNDF() const;

  static  SoLKalTrackSystem *GetCurInstancePtr() { return fgCurInstancePtr; }
  void SetCurInstancePtr(SoLKalTrackSystem *ksp) { fgCurInstancePtr = ksp; }
  inline void SetSitePtrToLastSite() { fCurSitePtr = static_cast<SoLKalTrackSite*>(this->Last()); }
  
  inline void      SetMomentum(Double_t m)       { fMomentum = m;    }
  inline Double_t  GetMomentum()         const   { return fMomentum; }
  inline void      SetTheta(Double_t t)          { fTheta = t;    }
  inline Double_t  GetTheta()            const   { return fTheta; }
  inline void      SetPhi(Double_t f)            { fPhi = f;    }
  inline Double_t  GetPhi()              const   { return fPhi; }
  inline void      SetVertexZ(Double_t z)        { fVertexZ = z;    }
  inline Double_t  GetVertexZ()          const   { return fVertexZ; }
  
  inline void      SetMass(Double_t m)           { fMass = m;    }
  inline Double_t  GetMass()             const   { return fMass; }
  inline void      SetCharge(Double_t q)         { fCharge = q;  }
  inline Double_t  GetCharge()           const   { return fCharge; }
  inline void      SetElectron(Bool_t tof)         { fIsElectron = tof;}
  inline Bool_t    IsElectron()         const   { return fIsElectron; }
  inline Bool_t    GetTrackStatus()      const   { return fIsGood; }
  inline void      SetTrackStatus(Bool_t a)      { fIsGood = a; }
  inline Int_t     GetAngleFlag()      const   { return fAngleFlag; }
  inline void      SetAngleFlag(Int_t a)      { fAngleFlag = a; }
  inline SeedType  GetSeedType()      const   { return fSeedType; }
  inline void      SetSeedType(SeedType a)      { fSeedType = a; }

  void             CheckTrackStatus();
  inline const Int_t & GetNMissingHits() const { return fNMissingHits; }
  inline void      AddMissingHits() { fNMissingHits++; }
  inline void      SetMissingHits(Int_t a) {fNMissingHits = a;}
  Int_t  GetNHits() const { return fNHits; }
  
  virtual Int_t Compare( const TObject* obj ) const;
  virtual Bool_t IsSortable () const { return kTRUE; }

  Double_t     fDeltaECX;
  Double_t     fDeltaECY;
  Double_t     fDeltaECE;

  private:
  SoLKalTrackSite   *fCurSitePtr;  // pointer to current site
  
  
  static SoLKalTrackSystem *fgCurInstancePtr;  //! currently active instance

  Double_t     fMass;        // mass [GeV]
  Double_t     fCharge;      // in unit of proton charge
  Bool_t       fIsElectron;  // kTRUE if the particle is a lepton, otherwise kFALSE
  Bool_t       fIsGood;      // for track fitting monitering 
  Int_t        fAngleFlag;   // 0 for the large angle 1 for the forward angle
  Int_t        fNMissingHits;
  Int_t        fNHits;
  Int_t        fNDF;
  Double_t     fMomentum;
  Double_t     fTheta;
  Double_t     fPhi;
  Double_t     fVertexZ;
  SeedType     fSeedType;
  
  Double_t     fChi2;
       
  
  ClassDef(SoLKalTrackSystem,1)  // Base class for Kalman Filter

};
#endif
