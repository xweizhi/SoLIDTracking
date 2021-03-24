//c++
#include <cassert>
#include <sstream>
#include <iostream>

//ROOT
#include "TMath.h"
#include "TVector2.h"
#include "TVector3.h"
//SoLIDTracking
#include "SoLIDECal.h"
#include "SoLIDTrackerSystem.h"
ClassImp(SoLIDECal)

//________________________________________________________________________________________
SoLIDECal::SoLIDECal( const char* name, const char* description,
                                 THaDetectorBase* parent)
:THaSubDetector(name,description,parent),
fIsLAECTriggered(kFALSE), fIsFAECTriggered(kFALSE), fLAECEdpCut(0.), fFAECEdpCut(0.)
{
  static const char* const here = "SoLIDECal";
  assert( name && parent );
  fCaloHits = new TClonesArray("SoLIDCaloHit", 10);

}
//________________________________________________________________________________________
SoLIDECal::~SoLIDECal()
{
  if( fIsSetup )
    RemoveVariables();
  delete fCaloHits;
}
//________________________________________________________________________________________
void SoLIDECal::Clear(Option_t* opt)
{
  fNLAECHits = 0;
  fNFAECHits = 0;
  fIsLAECTriggered = kFALSE;
  fIsFAECTriggered = kFALSE;
  fCaloHits->Clear(opt);
}
//________________________________________________________________________________________
Int_t SoLIDECal::Decode(const THaEvData& evdata)
{
  vector<Float_t> laXPos;
  vector<Float_t> laYPos;
  vector<Float_t> laEdp;
  vector<Int_t>   laPID;
  
  vector<Float_t> faXPos;
  vector<Float_t> faYPos;
  vector<Float_t> faEdp;
  vector<Int_t>   faPID;
  
  if( GetNLAECHits() == 0 && GetNFAECHits() == 0) {
  

    for( Int_t imod = 0; imod < fDetMap->GetSize(); ++imod ) {
      // Decode data
      THaDetMap::Module* d = fDetMap->GetModule(imod);
      // Read active channels of this module
      Int_t nchan = evdata.GetNumChan( d->crate, d->slot );
      for( Int_t ichan = 0; ichan < nchan; ++ichan ) {
        Int_t chan = evdata.GetNextChan( d->crate, d->slot, ichan );
        if( chan < d->lo or chan > d->hi ) continue; // not part of this detector


        Int_t nhit = evdata.GetNumHits( d->crate, d->slot, chan );
        for( Int_t ihit = 0; ihit < nhit; ++ihit ) {
          // The hit's data and raw data words hold the x and y coordinates,
          // respectively
          union FloatIntUnion {
            Float_t f;
            Int_t   i;
          } datx, daty;
          datx.i = evdata.GetData( d->crate, d->slot, chan, ihit );
          daty.i = evdata.GetRawData( d->crate, d->slot, chan, ihit );
          
          switch(imod){
            case kLAECPos:
              {
                laXPos.push_back(datx.f);
                laYPos.push_back(daty.f);
                break;
              }
            case kLAECEdp:
              if (datx.f > fLAECEdpCut) fIsLAECTriggered = kTRUE;
              laEdp.push_back(datx.f);
              laPID.push_back(TMath::Nint(daty.f));
              break;
            case kFAECPos:
              {
                faXPos.push_back(datx.f);
                faYPos.push_back(daty.f);
                break;
              }
            case kFAECEdp:
              if (datx.f > fFAECEdpCut) fIsFAECTriggered = kTRUE;
              faEdp.push_back(datx.f);
              faPID.push_back(TMath::Nint(daty.f));
              break;
            default:
              break;
          }
          
        }
      }
    } 
  }
  assert(laXPos.size() == laEdp.size() && faXPos.size() == faEdp.size());
  fNLAECHits = 0;
  fNFAECHits = 0;
  
  for (unsigned int i=0; i<laXPos.size(); i++) {
    Int_t thisID = 0;
    Int_t thisPID = TMath::Nint(laPID[i]);
    Float_t thisX  = laXPos[i];
    Float_t thisY  = laYPos[i];
    Float_t thisEdep = laEdp[i];
    
    if (abs(thisPID) == 11 || thisPID == 22) {
        SmearPosition(&thisX, &thisY, &thisID, 0);
        SmearEnergy(&thisEdep);
    }else{
        thisID = -1;
    }
    
    if (thisID < 0) continue; 
    new ( (*fCaloHits)[i]) SoLIDCaloHit( thisX, thisY, kLAEC, thisEdep, thisID);
    fNLAECHits++;
  }
  
  for (unsigned int i=0; i<faXPos.size(); i++){
    Int_t thisID = 0;
    Int_t thisPID = TMath::Nint(faPID[i]);
    Float_t thisX  = faXPos[i];
    Float_t thisY  = faYPos[i];
    Float_t thisEdep = faEdp[i];
    
    if (abs(thisPID) == 11 || thisPID == 22) {
        SmearPosition(&thisX, &thisY, &thisID, 0);
        SmearEnergy(&thisEdep);
    }else{
        SmearPosition(&thisX, &thisY, &thisID, 1);
        thisEdep = 0.; //no edep info if not ec
    }
    if (thisID < 0) continue; 
    new ( (*fCaloHits)[i+fNLAECHits]) SoLIDCaloHit(thisX, thisY, kFAEC, thisEdep, thisID);
    fNFAECHits++;
  }
  
  return kOK;
}
//________________________________________________________________________________________
THaAnalysisObject::EStatus SoLIDECal::Init( const TDatime& date )
{

  EStatus status = THaAnalysisObject::Init(date);

  if( status ){
  return fStatus = status;
  }
  return fStatus = kOK;
}
//_________________________________________________________________________________________
void SoLIDECal::Print( Option_t* /*opt*/ ) const
{

}
//__________________________________________________________________________________________
void SoLIDECal::PrintDataBase() const
{

}
//__________________________________________________________________________________________
Int_t SoLIDECal::Begin( THaRunBase* /*r*/ )
{
  return 0;
}
//_________________________________________________________________________________________
Int_t SoLIDECal::End( THaRunBase* /*r*/ )
{
  return 0;
}
//________________________________________________________________________________________
Int_t SoLIDECal::ReadDatabase( const TDatime& date )
{
  static const char* const here = "SoLIDECal::ReadDatabase";
  fIsInit = kFALSE;
  FILE* file = OpenFile( date );
  if( !file ) return kFileError;
  
  Int_t status = -1;
  fLAECZ = 0;
  fFAECZ = 0;
  fPosReso = -1.;//m
  fEReso = -1.;
#ifdef SIDIS
  fUseMRPC = 0;
  fUseSPD = 0;
#endif
  vector<Int_t>* laec_detmap_pos = 0;
  vector<Int_t>* laec_detmap_edp = 0;
  vector<Int_t>* faec_detmap_pos = 0;
  vector<Int_t>* faec_detmap_edp = 0;
  vector<Double_t>* spd_r_segment = 0;
  try{
   laec_detmap_pos = new vector<Int_t>;
   laec_detmap_edp = new vector<Int_t>;
   faec_detmap_pos = new vector<Int_t>;
   faec_detmap_edp = new vector<Int_t>;
   spd_r_segment   = new vector<Double_t>;
   const DBRequest request[] = {
         { "laec_detmap_pos",      laec_detmap_pos,   kIntV },
         { "laec_detmap_edp",      laec_detmap_edp,   kIntV },
         { "faec_detmap_pos",      faec_detmap_pos,   kIntV },
         { "faec_detmap_edp",      faec_detmap_edp,   kIntV },
         { "laec_z",               &fLAECZ,           kDouble,  0, 1 },
         { "faec_z",               &fFAECZ,           kDouble,  0, 1 },
         { "ec_pos_reso",          &fPosReso,         kDouble,  0, 1 },
         { "ec_energy_reso",       &fEReso,           kDouble,  0, 1 },
#ifdef SIDIS
         { "use_mrpc",             &fUseMRPC,         kInt,     0, 1 },
	     { "mrpc_pitch_width",     &fMRPCPitchWidth,  kDouble,  0, 1 },
         { "mrpc_n_sectors",       &fMRPCNSectors,    kInt,     1, 1 },
         { "mrpc_phi_reso",        &fMRPCPhiReso,     kDouble,  0, 1 },
	     { "mrpc_rmin",            &fMRPCRmin,        kDouble,  0, 1 },
	     { "use_spd",              &fUseSPD,         kInt,     0, 1 },
	     { "n_r_segment",          &fSPDNRSegment,    kInt,     0, 1 },
	     { "n_phi_segment",        &fSPDNPhiSegment,  kInt,     0, 1 },
	     { "spd_r_segment",        spd_r_segment,     kDoubleV},
#endif
         { 0 }
       };
   status = LoadDB( file, date, request, fPrefix );
   assert(fPosReso >= 0 && fEReso >= 0);
#ifdef SIDIS
    assert((fUseMRPC > 0 || fUseSPD > 0) && fUseMRPC*fUseSPD == 0);
    assert(fSPDRSegment.size() == 0);
    for (unsigned int i=0; i<spd_r_segment->size(); i++) fSPDRSegment.push_back(spd_r_segment->at(i));
    //fSPDRSegemtn contains the end point of these segments, so he size is +1 of fSPDNRSegment
    assert(fSPDRSegment.size() == fSPDNRSegment+1);
#endif
    if (status == kOK){
      if( FillDetMap( *laec_detmap_pos, THaDetMap::kDoNotClear, here ) <= 0 ) status = kInitError;
      if( FillDetMap( *laec_detmap_edp, THaDetMap::kDoNotClear, here ) <= 0 ) status = kInitError;
      if( FillDetMap( *faec_detmap_pos, THaDetMap::kDoNotClear, here ) <= 0 ) status = kInitError;
      if( FillDetMap( *faec_detmap_edp, THaDetMap::kDoNotClear, here ) <= 0 ) status = kInitError;
    }
    delete laec_detmap_pos;
    delete laec_detmap_edp;
    delete faec_detmap_pos;
    delete faec_detmap_edp;
  }catch(...) {
      delete laec_detmap_pos;
      delete laec_detmap_edp;
      delete faec_detmap_pos;
      delete faec_detmap_edp;
      fclose(file);
      throw;
  }
  
  fclose(file);
  if( status != kOK ) return status;
  
  for( Int_t imod = 0; imod < fDetMap->GetSize(); ++imod ) {
    THaDetMap::Module* d = fDetMap->GetModule(imod);
    assert( dynamic_cast<SoLIDTrackerSystem*>(GetMainDetector()) );
    SoLIDTrackerSystem *thisSystem = dynamic_cast<SoLIDTrackerSystem*>( GetMainDetector() );
    thisSystem->LoadDAQmodel(d);
    thisSystem->LoadDAQresolution(d);
    //only ADC for now
    d->MakeADC();
    UInt_t nchan = thisSystem->GetDAQnchan(d);
    if( d->hi >= nchan ) {
      Error( Here(here), "Detector map channel out of range for module "
          "cr/sl/lo/hi = %u/%u/%u/%u. Must be < %u. Fix database.",
          d->crate, d->slot, d->lo, d->hi, nchan );
      return kInitError;
    }
  }
#ifdef SIDIS
  fMRPCPhiCover = 2.*TMath::Pi() / fMRPCNSectors;
  fSPDPhiCover  = 2.*TMath::Pi() / fSPDNPhiSegment;
#endif
  fIsInit = kTRUE;

  return kOK;
}
//___________________________________________________________________________________________________
Int_t SoLIDECal::DefineVariables( EMode mode )
{
  if( mode == kDefine && fIsSetup ) return kOK;
  fIsSetup = ( mode == kDefine );
  
  Int_t ret;
  RVarDef nonmcvars[] = {
      { "hitCalo.x",       "hit x coordinate on EC",          "fCaloHits.SoLIDCaloHit.fXPos"         },
      { "hitCalo.y",       "hit y coordinate on EC",          "fCaloHits.SoLIDCaloHit.fYPos"         },
      { "hitCalo.id",      "EC id of the hit",                "fCaloHits.SoLIDCaloHit.fECID"         },
      { "hitCalo.E",       "hit energy measurement on EC",    "fCaloHits.SoLIDCaloHit.fEdp"          },
      { 0 }
    };
  ret = DefineVarsFromList( nonmcvars, mode );
  
  return ret;
}
//___________________________________________________________________________________________________
void SoLIDECal::SmearPosition(Float_t *x, Float_t *y, Int_t* id, Int_t mode)
{
#ifdef SIDIS
  //only in SIDIS forward angle, whre ec hit position will be replaced by hit on MRPC
  if (mode != 0){
    if (fUseMRPC){
        //if FAEC, use hit on MRPC to replace hit on EC

        //which MRPC sector
        //assume the first sector is centered at 0 deg
        double phi = atan2(*y, *x) + fMRPCPhiCover/2.;
        phi = TVector2::Phi_0_2pi(phi);
        int sector = (int)(phi / fMRPCPhiCover);
        assert(sector < fMRPCNSectors);

        double phiCenter = sector*fMRPCPhiCover;

        //rotate to the frame that with x-axis parallel to the symmetric axis of the sector

        double tmpx = cos(-1.*phiCenter)*(*x) + -1.*sin(-1.*phiCenter)*(*y);
        double tmpy = sin(-1.*phiCenter)*(*x) +     cos(-1.*phiCenter)*(*y);

        assert(tmpx > fMRPCRmin - 0.01); //should not happen

        int ipitch = (tmpx - fMRPCRmin)/fMRPCPitchWidth;
        tmpx  = fMRPCRmin + (ipitch + 0.5)*fMRPCPitchWidth; //using the pitch center
        tmpy += gRandom->Gaus(0., fMRPCPhiReso);

        //rotate back to the lab frame
        *x = cos(phiCenter)*tmpx + -1.*sin(phiCenter)*tmpy;
        *y = sin(phiCenter)*tmpx +     cos(phiCenter)*tmpy;
        *id = 1;
    }
    else if (fUseSPD){
        double thisx = *x;
        double thisy = *y;
        double thisr = sqrt(thisx*thisx + thisy*thisy);
        *id = -1;
        for (unsigned int i=0; i<fSPDRSegment.size()-1; i++){
            if (thisr < fSPDRSegment[i+1] && thisr >= fSPDRSegment[i]){
                *id = i+1;
                thisr = (fSPDRSegment[i+1] + fSPDRSegment[i])/2.;
                break;
            }
        }
        int phiid = 0;
        double phi0 = 0;
        double phiCover = fSPDPhiCover*180./TMath::Pi();
        double thisPhi = atan2(thisy, thisx);
        thisPhi = TVector2::Phi_mpi_pi(thisPhi)*180./TMath::Pi();
        
        if (thisPhi>=90) {
            phiid=int((thisPhi-90)/phiCover);
            phi0 = (phiid+0.5)*phiCover + 90;
        }
        else {
            phiid=int((thisPhi+360-90)/phiCover);
            phi0 = (phiid+0.5)*phiCover - 360 + 90;
        }
        
        phi0 *= TMath::Pi()/180.;
        
        *x = thisr*cos(phi0);
        *y = thisr*sin(phi0);
    }
    return; //no need to continue
  }
#endif

  *x += gRandom->Gaus(0, fPosReso);
  *y += gRandom->Gaus(0, fPosReso);
  *id = 0;
}
//___________________________________________________________________________________________________
inline void SoLIDECal::SmearEnergy(Float_t *energy)
{
  *energy *= gRandom->Gaus(1., fEReso/TMath::Sqrt(*energy)) ;

}
























