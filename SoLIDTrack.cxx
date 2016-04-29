//SoLIDTracking
#include "SoLIDTrack.h"

ClassImp(SoLIDTrack)
#ifdef MCDATA
ClassImp(SoLIDMCTrack)
#endif


//----------------------SoLIDTrack------------------------------//
//_______________________________________________________________
Int_t SoLIDTrack::Compare( const TObject* obj ) const
{
   // Used for sorting tracks in a TSeqCollection (e.g. TList, TClonesArray).
   //a track is considered better than the other if it has more hits 
   //and has a smaller chi2 per degree of freedom
   assert( dynamic_cast<const SoLIDTrack*>(obj) );
   const SoLIDTrack* rhs = static_cast<const SoLIDTrack*>(obj);

    if (fAngleFlag < rhs->GetAngleFlag()){
    return -1;
  }else if (fAngleFlag > rhs->GetAngleFlag()){
    return 1;
  }else{
    if (fCharge < rhs->GetCharge()){
      return -1;
    }else if (fCharge > rhs->GetCharge()){
      return 1;
    }else{
      if (GetNHits() < rhs->GetNHits()){
        return 1;
      }else if (GetNHits() > rhs->GetNHits()){
        return -1;
      }else{
        if (GetChi2() > rhs->GetChi2()){
          return 1;
        }else{
          return -1;
        }
      }
    }
  }

}
//________________________________________________________________
Double_t SoLIDTrack::GetHitInfo(UInt_t i, UInt_t type) const
{
  while (i>= fHits.size()) { i -= fHits.size(); } 
  
  if (type == 0){ //x coordinate info for the ith hit
    return fHits.at(i)->GetX();
  }
  else if (type == 1){ //y coordinate info for ith hit
    return fHits.at(i)->GetY();
  }
  else if (type == 2){ //z coordinate info for the ith hit
    return fHits.at(i)->GetZ();
  }
  else if (type == 3){ //r coordinate info for the ith hit
    return fHits.at(i)->GetR();
  }
  else if (type == 4){ //phi coordinate info for the ith hit
    return fHits.at(i)->GetPhi();
  }
  else if (type == 5){ //the tracker ID of the ith hit
    return fHits.at(i)->GetTrackerID();
  }
  else if (type == 6){ //the chamber ID of the ith hit
    return fHits.at(i)->GetChamberID();
  }
  else { return 0; }
  
}

//----------------------SoLIDMCTrack--------------------------------//
//__________________________________________________________________
#ifdef MCDATA
Double_t SoLIDMCTrack::GetMCHitInfo(UInt_t i, UInt_t type) const 
{
  while (i>= fHits.size()) { i -= fHits.size(); } 
  
  if (type == 0){ //x coordinate info for the ith hit
    return dynamic_cast<SoLIDMCGEMHit*>(fHits.at(i))->IsSignalHit();
  }
  else { return 0; }
}
//__________________________________________________________________
Int_t SoLIDMCTrack::GetNMCHits()
{
  UInt_t count = 0;
  for (UInt_t i=0; i<fHits.size(); i++){
    if (dynamic_cast<SoLIDMCGEMHit*>(fHits.at(i))->IsSignalHit()) {
      count++;
    }
  }
  return count;
}
#endif
