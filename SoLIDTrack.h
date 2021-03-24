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
  fPID(11), fAngleFlag(0),fVertexZ(0.),fTheta(0.), fPhi(0.), fNDF(1), fSelectFlag(0), fIsElectron(true) {}
  ~SoLIDTrack(){fHits.clear();}
  
  virtual Int_t Compare( const TObject* obj ) const;
  virtual Bool_t IsSortable () const { return kTRUE; }
  
  void         SetCoarseFitStatus(Bool_t is) { fIsCoarseFitted = is; }
  void         SetCoarseChi2(Double_t chi2)  { fCoarseChi2 = chi2; }
  void         SetFineFirStatus(Bool_t is)   { fIsFineFitted = is; }
  void         SetFineChi2(Double_t chi2)    { fFineChi2 = chi2; }
  void         SetStatus(Bool_t is)          { fStatus = is; }
  void         SetCharge(Double_t charge)    { fCharge = charge; }
  void         SetMass(Double_t mass)        { fMass = mass; }
  void         SetPID(Int_t pid)             { fPID = pid; }
  void         SetAngleFlag(Int_t angleflag) { fAngleFlag = angleflag; }
  void         SetECX(Double_t ecx)          { fECX = ecx; }
  void         SetECY(Double_t ecy)          { fECY = ecy; }
  void         SetECE(Double_t ece)          { fECE = ece; }
  void         SetECEx(Double_t ecex)        { fECEx = ecex; }
  void         SetECEy(Double_t ecey)        { fECEy = ecey; }
  void         SetMomentum(Double_t momentum){ fMomentum = momentum; }
  void         SetVertexZ(Double_t vertexz)  { fVertexZ = vertexz; }
  void         SetTheta(Double_t theta)      { fTheta = theta; }
  void         SetPhi(Double_t phi)          { fPhi = phi; } 
  void         SetNDF(Int_t ndf)             { fNDF = ndf; }
  void         SetSelectFlag(Int_t f)        { fSelectFlag = f; }
  
  void         SetBackTheta(Double_t theta)  { fBackTheta = theta; }
  void         SetBackPhi(Double_t phi)      { fBackPhi = phi; }
  void         SetBackX(Double_t x)          { fBackX = x; }
  void         SetBackY(Double_t y)          { fBackY = y; }
  void         SetBackMom(Double_t mom)      { fBackMom = mom; }

  void         SetElectron(Bool_t b)         { fIsElectron = b; }
  
 
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
  Double_t     GetECX()                const { return fECX; }
  Double_t     GetECY()                const { return fECY; }
  Double_t     GetECE()                const { return fECE; }
  Double_t     GetECEx()               const { return fECEx; }
  Double_t     GetECEy()               const { return fECEy; }
  Double_t     GetMomentum()           const { return fMomentum; }
  Double_t     GetVertexZ()            const { return fVertexZ; }
  Double_t     GetTheta()              const { return fTheta; }
  Double_t     GetPhi()                const { return fPhi; }
  Int_t        GetNDF()                const { return fNDF; }
  Double_t     GetChi2perNDF()         const { return fIsFineFitted ? 
               fFineChi2/(Double_t)fNDF : fCoarseChi2/(Double_t)fNDF; }
  Int_t        GetSelectFlag()         const { return fSelectFlag; }
  
  void         AddHit(SoLIDGEMHit* theHit) { fHits.push_back(theHit); }
  
  Double_t     GetHitInfo(UInt_t i, UInt_t type) const;
  
  Double_t     GetHitX0()              const { return GetHitInfo(0, 0); }
  Double_t     GetHitX1()              const { return GetHitInfo(1, 0); }
  Double_t     GetHitX2()              const { return GetHitInfo(2, 0); }
  Double_t     GetHitX3()              const { return GetHitInfo(3, 0); }
  Double_t     GetHitX4()              const { return GetHitInfo(4, 0); }
  Double_t     GetHitX5()              const { return GetHitInfo(5, 0); }
  
  Double_t     GetPredHitX0()          const { return GetHitInfo(0, 7); }
  Double_t     GetPredHitX1()          const { return GetHitInfo(1, 7); }
  Double_t     GetPredHitX2()          const { return GetHitInfo(2, 7); }
  Double_t     GetPredHitX3()          const { return GetHitInfo(3, 7); }
  Double_t     GetPredHitX4()          const { return GetHitInfo(4, 7); }
  Double_t     GetPredHitX5()          const { return GetHitInfo(5, 7); }
  
  Double_t     GetHitY0()              const { return GetHitInfo(0, 1); }
  Double_t     GetHitY1()              const { return GetHitInfo(1, 1); }
  Double_t     GetHitY2()              const { return GetHitInfo(2, 1); }
  Double_t     GetHitY3()              const { return GetHitInfo(3, 1); }
  Double_t     GetHitY4()              const { return GetHitInfo(4, 1); }
  Double_t     GetHitY5()              const { return GetHitInfo(5, 1); }
  
  Double_t     GetPredHitY0()          const { return GetHitInfo(0, 8); }
  Double_t     GetPredHitY1()          const { return GetHitInfo(1, 8); }
  Double_t     GetPredHitY2()          const { return GetHitInfo(2, 8); }
  Double_t     GetPredHitY3()          const { return GetHitInfo(3, 8); }
  Double_t     GetPredHitY4()          const { return GetHitInfo(4, 8); }
  Double_t     GetPredHitY5()          const { return GetHitInfo(5, 8); }
  
  Double_t     GetHitTracker0()        const { return GetHitInfo(0, 5); }
  Double_t     GetHitTracker1()        const { return GetHitInfo(1, 5); }
  Double_t     GetHitTracker2()        const { return GetHitInfo(2, 5); }
  Double_t     GetHitTracker3()        const { return GetHitInfo(3, 5); }
  Double_t     GetHitTracker4()        const { return GetHitInfo(4, 5); }
  Double_t     GetHitTracker5()        const { return GetHitInfo(5, 5); }
  
  Double_t     GetPredHiteX0()          const { return GetHitInfo(0, 9); }
  Double_t     GetPredHiteX1()          const { return GetHitInfo(1, 9); }
  Double_t     GetPredHiteX2()          const { return GetHitInfo(2, 9); }
  Double_t     GetPredHiteX3()          const { return GetHitInfo(3, 9); }
  Double_t     GetPredHiteX4()          const { return GetHitInfo(4, 9); }
  Double_t     GetPredHiteX5()          const { return GetHitInfo(5, 9); }
  
  Double_t     GetPredHiteY0()          const { return GetHitInfo(0, 10); }
  Double_t     GetPredHiteY1()          const { return GetHitInfo(1, 10); }
  Double_t     GetPredHiteY2()          const { return GetHitInfo(2, 10); }
  Double_t     GetPredHiteY3()          const { return GetHitInfo(3, 10); }
  Double_t     GetPredHiteY4()          const { return GetHitInfo(4, 10); }
  Double_t     GetPredHiteY5()          const { return GetHitInfo(5, 10); }

  Double_t     GetHitPX0()              const { return GetHitInfo(0, 11); }
  Double_t     GetHitPX1()              const { return GetHitInfo(1, 11); }
  Double_t     GetHitPX2()              const { return GetHitInfo(2, 11); }
  Double_t     GetHitPX3()              const { return GetHitInfo(3, 11); }
  Double_t     GetHitPX4()              const { return GetHitInfo(4, 11); }
  Double_t     GetHitPX5()              const { return GetHitInfo(5, 11); }

  Double_t     GetHitPY0()              const { return GetHitInfo(0, 12); }
  Double_t     GetHitPY1()              const { return GetHitInfo(1, 12); }
  Double_t     GetHitPY2()              const { return GetHitInfo(2, 12); }
  Double_t     GetHitPY3()              const { return GetHitInfo(3, 12); }
  Double_t     GetHitPY4()              const { return GetHitInfo(4, 12); }
  Double_t     GetHitPY5()              const { return GetHitInfo(5, 12); }

  Double_t     GetHitPZ0()              const { return GetHitInfo(0, 13); }
  Double_t     GetHitPZ1()              const { return GetHitInfo(1, 13); }
  Double_t     GetHitPZ2()              const { return GetHitInfo(2, 13); }
  Double_t     GetHitPZ3()              const { return GetHitInfo(3, 13); }
  Double_t     GetHitPZ4()              const { return GetHitInfo(4, 13); }
  Double_t     GetHitPZ5()              const { return GetHitInfo(5, 13); }
  
  Double_t     GetHit0Chi2()             const { return GetHitInfo(0, 14); }
  Double_t     GetHit1Chi2()             const { return GetHitInfo(1, 14); }
  Double_t     GetHit2Chi2()             const { return GetHitInfo(2, 14); }
  Double_t     GetHit3Chi2()             const { return GetHitInfo(3, 14); }
  Double_t     GetHit4Chi2()             const { return GetHitInfo(4, 14); }
  Double_t     GetHit5Chi2()             const { return GetHitInfo(5, 14); }

  Double_t     GetHit0QU()               const { return GetHitInfo(0, 15); }
  Double_t     GetHit1QU()               const { return GetHitInfo(1, 15); }
  Double_t     GetHit2QU()               const { return GetHitInfo(2, 15); }
  Double_t     GetHit3QU()               const { return GetHitInfo(3, 15); }
  Double_t     GetHit4QU()               const { return GetHitInfo(4, 15); }
  Double_t     GetHit5QU()               const { return GetHitInfo(5, 15); }

  Double_t     GetHit0QV()               const { return GetHitInfo(0, 16); }
  Double_t     GetHit1QV()               const { return GetHitInfo(1, 16); }
  Double_t     GetHit2QV()               const { return GetHitInfo(2, 16); }
  Double_t     GetHit3QV()               const { return GetHitInfo(3, 16); }
  Double_t     GetHit4QV()               const { return GetHitInfo(4, 16); }
  Double_t     GetHit5QV()               const { return GetHitInfo(5, 16); }

  Double_t     GetBackTheta()           const { return fBackTheta; }
  Double_t     GetBackPhi()             const { return fBackPhi; }
  Double_t     GetBackX()               const { return fBackX; }
  Double_t     GetBackY()               const { return fBackY; }
  Double_t     GetBackMom()             const { return fBackMom; }

  Bool_t       IsElectron()               const { return fIsElectron; }
  
  
  static bool SortHitZ(const SoLIDGEMHit* a, const SoLIDGEMHit* b);
  void SortHits();
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
  Double_t fECX;
  Double_t fECY;
  Double_t fECE;
  Double_t fECEx;
  Double_t fECEy;
  Double_t fMomentum;
  Double_t fVertexZ;
  Double_t fTheta;
  Double_t fPhi;
  Int_t fNDF;
  Int_t fSelectFlag;                  //whether the track is the best within tracks that share
                                      // common hits

  Bool_t fIsElectron;  
  //back track info
  Double_t fBackTheta;
  Double_t fBackPhi;
  Double_t fBackX;
  Double_t fBackY;
  Double_t fBackMom;
  
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























