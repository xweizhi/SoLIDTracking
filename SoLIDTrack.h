#ifndef ROOT_SoLID_Track
#define ROOT_SoLID_Track
//c++
#include <vector>
#include <iostream>
//SoLIDTracking
#include "SoLIDGEMHit.h"
#define kINFINITY 1e12
class SoLIDTrack : public TObject
{
  public:
  SoLIDTrack():fCoarseChi2(kINFINITY), fFineChi2(kINFINITY), fIsCoarseFitted(kFALSE), 
  fIsFineFitted(kFALSE), fStatus(1), fCharge(-1.), fMass(0.51e-3), 
  fPID(11), fAngleFlag(0),fVertexZ(0.),fTheta(0.) {}
  ~SoLIDTrack(){}
  
  virtual Int_t Compare( const TObject* obj ) const;
  virtual Bool_t IsSortable () const { return kTRUE; }
  
  void         SetCoarseFitStatus(Bool_t is) { fIsCoarseFitted = is; }
  void         SetCoarseChi2(Double_t chi2)  { fCoarseChi2 = chi2; }
  void         SetFineFirStatus(Bool_t is)   { fIsFineFitted = is; }
  void         SetFineChi2(Double_t chi2)    { fFineChi2 = chi2; }
  void         SetStatus(Bool_t is)          { fStatus = is; }
  void         SetCharge(Double_t charge)    { fCharge = charge; }
  void         SetMaxx(Double_t mass)        { fMass = mass; }
  void         SetPID(Int_t pid)             { fPID = pid; }
  void         SetAngleFlag(Int_t angleflag) { fAngleFlag = angleflag; }
  void         SetMomMax(Double_t mommax)    { fMomMax = mommax; }
  void         SetMomMin(Double_t mommin)    { fMomMin = mommin; }
  void         SetThetaMax(Double_t thetamax){ fThetaMax = thetamax; }
  void         SetThetaMin(Double_t thetamin){ fThetaMin = thetamin; }
  void         SetMomentum(Double_t momentum){ fMomentum = momentum; }
  void         SetVertexZ(Double_t vertexz)  { fVertexZ = vertexz; }
  void         SetTheta(Double_t theta)      { fTheta = theta; }
  void         SetPhi(Double_t phi)          { fPhi = phi; } 
   
  Bool_t       IsCoarseFitted()        const { return fIsCoarseFitted; }
  Bool_t       IsFineFitted()          const { return fIsFineFitted; }
  Double_t     GetChi2()               const { return fIsFineFitted ? fFineChi2 : fCoarseChi2; }
  UInt_t       GetNHits()              const { return fHits.size(); }
  SoLIDGEMHit* GetHit(UInt_t i)        const { return i<fHits.size() ?
                                               fHits.at(i) : 0; } 
  Bool_t       GetStatus()             const { return fStatus; }
  Double_t     GetCharge()             const { return fCharge; }
  Double_t     GetMass()               const { return fMass; }
  Double_t     GetPID()                const { return fPID; }
  Int_t        GetAngleFlag()          const { return fAngleFlag; }
  Double_t     GetMomMax()             const { return fMomMax; }
  Double_t     GetMomMin()             const { return fMomMin; }
  Double_t     GetThetaMax()           const { return fThetaMax; }
  Double_t     GetThetaMin()           const { return fThetaMin; }
  Double_t     GetMomentum()           const { return fMomentum; }
  Double_t     GetVertexZ()            const { return fVertexZ; }
  Double_t     GetTheta()              const { return fTheta; }
  Double_t     GetPhi()                const { return fPhi; }
  
  void         AddHit(SoLIDGEMHit* theHit) { fHits.push_back(theHit); }
  
  Double_t     GetHitInfo(UInt_t i, UInt_t type) const;
  
  Double_t     GetHitX0()              const { return GetHitInfo(0, 0); }
  Double_t     GetHitX1()              const { return GetHitInfo(1, 0); }
  Double_t     GetHitX2()              const { return GetHitInfo(2, 0); }
  Double_t     GetHitX3()              const { return GetHitInfo(3, 0); }
  Double_t     GetHitX4()              const { return GetHitInfo(4, 0); }
  Double_t     GetHitX5()              const { return GetHitInfo(5, 0); }
  
  Double_t     GetHitY0()              const { return GetHitInfo(0, 1); }
  Double_t     GetHitY1()              const { return GetHitInfo(1, 1); }
  Double_t     GetHitY2()              const { return GetHitInfo(2, 1); }
  Double_t     GetHitY3()              const { return GetHitInfo(3, 1); }
  Double_t     GetHitY4()              const { return GetHitInfo(4, 1); }
  Double_t     GetHitY5()              const { return GetHitInfo(5, 1); }
  
  Double_t     GetHitTracker0()        const { return GetHitInfo(0, 5); }
  Double_t     GetHitTracker1()        const { return GetHitInfo(1, 5); }
  Double_t     GetHitTracker2()        const { return GetHitInfo(2, 5); }
  Double_t     GetHitTracker3()        const { return GetHitInfo(3, 5); }
  Double_t     GetHitTracker4()        const { return GetHitInfo(4, 5); }
  Double_t     GetHitTracker5()        const { return GetHitInfo(5, 5); }
          
  protected:
  Double_t fCoarseChi2;
  Double_t fFineChi2;
  Bool_t   fIsCoarseFitted;
  Bool_t   fIsFineFitted;
  Bool_t   fStatus;                   //status of the track, whether it is still a good track
  Double_t fCharge;
  Double_t fMass;
  Int_t    fPID;
  Int_t    fAngleFlag;               //0 for SIDIS large angle, 1 for SIDIS forwardangle, 0 for PVDIS
  Double_t fMomMax;
  Double_t fMomMin;
  Double_t fThetaMax;
  Double_t fThetaMin;
  Double_t fMomentum;
  Double_t fVertexZ;
  Double_t fTheta;
  Double_t fPhi;
  std::vector<SoLIDGEMHit*> fHits;
  
  ClassDef(SoLIDTrack, 1)
};

//----------------------------------------------------------------------------------------//
#ifdef MCDATA
class SoLIDMCTrack : public SoLIDTrack
{
  public:
  SoLIDMCTrack():SoLIDTrack() {}
  ~SoLIDMCTrack() {}
  
  Double_t    GetMCHitInfo(UInt_t i, UInt_t type) const;
  
  Double_t    IsSignalHit0()         const { return GetMCHitInfo(0, 0); }
  Double_t    IsSignalHit1()         const { return GetMCHitInfo(1, 0); }
  Double_t    IsSignalHit2()         const { return GetMCHitInfo(2, 0); }
  Double_t    IsSignalHit3()         const { return GetMCHitInfo(3, 0); }
  Double_t    IsSignalHit4()         const { return GetMCHitInfo(4, 0); }
  Double_t    IsSignalHit5()         const { return GetMCHitInfo(5, 0); }
  
  Int_t       GetNMCHits();
  ClassDef(SoLIDMCTrack, 1)
};

#endif //MCDATA
#endif //ROOT_SoLID_Track























