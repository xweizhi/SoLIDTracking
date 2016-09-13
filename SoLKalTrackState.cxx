//c++
#include <cassert>
//ROOT
#include "TMath.h"
//SoLIDTracking
#include "SoLKalTrackState.h"
#include "SoLKalFieldStepper.h"

ClassImp(SoLKalTrackState)

//_________________________________________________________________
SoLKalTrackState::SoLKalTrackState(Int_t type, Int_t p)
: SoLKalMatrix(p,1), fType(type), fSitePtr(nullptr),
  fF(p,p), fFt(p,p), fQ(p,p), fC(p,p), fAttemptState(nullptr)
{
   fF.UnitMatrix();
   fFt.UnitMatrix();
   fZ0 = 0.;
   //fFieldStepper = SoLKalFieldStepper::GetInstance();
}
//_________________________________________________________________
SoLKalTrackState::SoLKalTrackState(const SoLKalMatrix &sv,
                       Int_t type, Int_t p)
: SoLKalMatrix(sv), fType(type), fSitePtr(nullptr),
  fF(p,p), fFt(p,p), fQ(p,p), fC(p,p), fAttemptState(nullptr)
{
   fF.UnitMatrix();
   fFt.UnitMatrix();
   fZ0 = 0.;
   //fFieldStepper = SoLKalFieldStepper::GetInstance();
}
//_________________________________________________________________
SoLKalTrackState::SoLKalTrackState(const SoLKalMatrix &sv, const SoLKalMatrix &c,
                       Int_t type, Int_t p)
: SoLKalMatrix(sv), fType(type), fSitePtr(nullptr), fF(p,p), fFt(p,p),
  fQ(p,p), fC(c), fAttemptState(nullptr)
{
   fF.UnitMatrix();
   fFt.UnitMatrix();
   fZ0 = 0.;
   //fFieldStepper = SoLKalFieldStepper::GetInstance();
}
//__________________________________________________________________
SoLKalTrackState::SoLKalTrackState(const SoLKalMatrix &sv, const SoLKalTrackSite &site,
                       Int_t type, Int_t p)
: SoLKalMatrix(sv), fType(type), fSitePtr((SoLKalTrackSite *)&site),
  fF(p,p), fFt(p,p), fQ(p,p), fC(p,p), fAttemptState(nullptr)
{
   fF.UnitMatrix();
   fFt.UnitMatrix();
   fZ0 = ((SoLKalTrackSite *)&site)->GetZ();
   //fFieldStepper = SoLKalFieldStepper::GetInstance();
}

//__________________________________________________________________
SoLKalTrackState::SoLKalTrackState(const SoLKalMatrix &sv, const SoLKalMatrix &c,
                       const SoLKalTrackSite &site, Int_t type, Int_t p)
: SoLKalMatrix(sv), fType(type), fSitePtr((SoLKalTrackSite *)&site),
  fF(p,p), fFt(p,p), fQ(p,p), fC(c), fAttemptState(nullptr)
{
   fF.UnitMatrix();
   fFt.UnitMatrix();
   fZ0 = ((SoLKalTrackSite *)&site)->GetZ();
   //fFieldStepper = SoLKalFieldStepper::GetInstance();
}
//___________________________________________________________________
SoLKalTrackState::~SoLKalTrackState()
{
  if (fAttemptState != nullptr && (&fAttemptState->GetSite()) == nullptr){
    delete fAttemptState;
  }
}
//___________________________________________________________________
void SoLKalTrackState::Propagate(SoLKalTrackSite &to)
{
   // Calculate 
   //    prea:  predicted state vector      : a^k-1_k = f_k-1(a_k-1)
   //    fF:    propagator derivative       : F_k-1   = (@f_k-1/@a_k-1)
   //    fQ:    process noise from k-1 to k : Q_k-1)

   SoLKalTrackState &prea    = MoveTo(to,fF,fQ);
   SoLKalTrackState *preaPtr = &prea;

   fFt = SoLKalMatrix(SoLKalMatrix::kTransposed, fF);

   // Calculate covariance matrix

   SoLKalMatrix preC = fF * fC * fFt + fQ;
   if (preC.Determinant() == 0) {
     cout<<"SoLKalTrackState::Propagate"<<endl;
     fF.Print();
     fC.Print();
     fQ.Print();
   }
   // Set predicted state vector and covariance matrix to next site

   prea.SetCovMat(preC);
   to.Add(preaPtr);
   to.SetOwner();
}
//____________________________________________________________________
SoLKalTrackState * SoLKalTrackState::MoveTo(SoLKalTrackSite  &to,
                                        SoLKalMatrix &F,
                                        SoLKalMatrix *QPtr) const
{
   if (QPtr) {
      const SoLKalTrackSite &from   = static_cast<const SoLKalTrackSite &>(GetSite());
            SoLKalTrackSite &siteto = static_cast<SoLKalTrackSite &>(to);

      SoLKalMatrix sv(kSdim,1);
      SoLKalFieldStepper::GetInstance()->Transport(from, to, sv, F, *QPtr);
      return new SoLKalTrackState(sv, siteto, SoLKalTrackSite::kPredicted, kSdim);
   } else {
     return nullptr;
   }
}
//_____________________________________________________________________
SoLKalTrackState & SoLKalTrackState::MoveTo(SoLKalTrackSite  &to,
                                        SoLKalMatrix &F,
                                        SoLKalMatrix &Q) const
{
   return *MoveTo(to, F, &Q);
}
//______________________________________________________________________
SoLKalTrackState* SoLKalTrackState::PredictSVatZ(Double_t &z)
{
   //simply pass the current filtered state vector to the next detector
   //it is up to the track finder to decide whether we have a hit in the 
   //next measurement layer
   
   SoLKalTrackState &prea    = MoveToZ(z,fF,fQ);
   SoLKalTrackState *preaPtr = &prea;

   fFt = SoLKalMatrix(SoLKalMatrix::kTransposed, fF);

   // Calculate covariance matrix

   SoLKalMatrix preC = fF * fC * fFt + fQ;

   prea.SetCovMat(preC);
   
   //only a filtered state can be used to predict hit position on the next layer
   if (fType == SoLKalTrackSite::kFiltered) fAttemptState = preaPtr;
   
   return preaPtr;
}
//______________________________________________________________________
void SoLKalTrackState::InitPredictSV()
{
  assert(fType != SoLKalTrackSite::kPredicted);
  if (fAttemptState == nullptr){
    SoLKalMatrix sv(kSdim,1);
    for (Int_t i=0; i<kSdim; i++) { sv(i, 0) = (*this)(i, 0); }
    fAttemptState = new SoLKalTrackState(sv, SoLKalTrackSite::kPredicted, kSdim);
    fAttemptState->SetZ0(this->GetZ0());
  }else{
    return;
  }
}
//______________________________________________________________________
SoLKalTrackState* SoLKalTrackState::PredictSVatNextZ(Double_t &z)
{
  //this function can be called only if the PredictSVatZ has been called
  //which is the first attempt to find hits on the next measurement site
  assert(fAttemptState != nullptr && fType != SoLKalTrackSite::kPredicted);
  SoLKalMatrix thisSV (kSdim, 1);
  SoLKalMatrix thisF (kSdim, kSdim);
  SoLKalMatrix thisQ (kSdim, kSdim);
  
  SoLKalFieldStepper::GetInstance()->Transport(*fAttemptState, z, thisSV, thisF, thisQ);
  
  for (Int_t i=0; i<kSdim; i++) { (*fAttemptState)(i, 0) = thisSV(i, 0); }
  
  fF = thisF*fF;
  fFt = SoLKalMatrix(SoLKalMatrix::kTransposed, fF);
  fQ = fQ + thisQ;
  
  SoLKalMatrix preC = fF * fC * fFt + fQ;
  fAttemptState->SetCovMat(preC);
  fAttemptState->SetZ0(SoLKalFieldStepper::GetInstance()->GetTrackPosAtZ());
  return fAttemptState;
}
//______________________________________________________________________
SoLKalTrackState * SoLKalTrackState::MoveToZ(Double_t z,
                                        SoLKalMatrix &F,
                                        SoLKalMatrix *QPtr) const
{
   if (QPtr) {
     const SoLKalTrackSite &from   = static_cast<const SoLKalTrackSite &>(GetSite());
     SoLKalMatrix sv(kSdim,1); 
     SoLKalFieldStepper::GetInstance()->Transport(from, z, sv, F, *QPtr);
     SoLKalTrackState* thisState =  new SoLKalTrackState(sv, SoLKalTrackSite::kPredicted, kSdim);
     thisState->SetZ0(SoLKalFieldStepper::GetInstance()->GetTrackPosAtZ());
     return thisState;
   } else {
     return nullptr;
   }
}
//_____________________________________________________________________
SoLKalTrackState & SoLKalTrackState::MoveToZ(Double_t z,
                                        SoLKalMatrix &F,
                                        SoLKalMatrix &Q) const
{
   return *MoveToZ(z, F, &Q);
}
//______________________________________________________________________
void SoLKalTrackState::CalcDir(TVector3 &dir) const {
    // Calculates a direction unit vector from this state.
    // dir: track direction (return value).

  CalcDir(dir, (*this));
}
//______________________________________________________________________
void SoLKalTrackState::CalcDir(TVector3 &dir, const SoLKalMatrix &sv) {
    // Calculates a direction unit vector from a state vector given by function
    // parameter sv.
    // dir: track direction (return value)
    // sv:  state vector to calculate direction from.

  Double_t tanx = sv(2,0);
  Double_t tany = sv(3,0);
  Double_t qp   = sv(4,0);

    dir.SetZ( 1./(TMath::Abs(qp) * TMath::Sqrt(tanx*tanx + tany*tany + 1. )) );
    dir.SetX(tanx * dir.Z());
    dir.SetY(tany * dir.Z());
    dir = dir.Unit();
}
//________________________________________________________________________
void SoLKalTrackState::CalcMomVec(TVector3 &dir) const {
    // Calculates the momentum vector from this state.
    // dir: track momentum (return value).

  CalcMomVec(dir, (*this));
}
//_______________________________________________________________________
void SoLKalTrackState::CalcMomVec(TVector3 &dir, const SoLKalMatrix &sv) {
    // Calculates momentum vector from a state vector given by function
    // parameter sv.
    // dir: track momentum (return value)
    // sv:  state vector to calculate direction from.

  Double_t tanx = sv(2,0);
  Double_t tany = sv(3,0);
  Double_t qp   = sv(4,0);

    dir.SetZ( 1./(TMath::Abs(qp) * TMath::Sqrt(tanx*tanx + tany*tany + 1)) );
    dir.SetX(tanx * dir.Z());
    dir.SetY(tany * dir.Z());
}
//________________________________________________________________________                         






















