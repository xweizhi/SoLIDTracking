//c++
#include <cassert>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>

//ROOT
#include "TVector2.h"

//SoLID Tracking
#include "SoLIDGEMTracker.h"
#include "SoLIDTrackerSystem.h" 
#include "SoLIDUtility.h"

#define OUTMAX 4096
ClassImp(SoLIDGEMTracker)

using namespace std;
SoLIDGEMTracker::SoLIDGEMTracker(Int_t iGEM, const char* name, const char* description,
                                 THaDetectorBase* parent)
  : THaSubDetector(name,description,parent), fTrackerID(iGEM), fNHits(0)
{
  static const char* const here = "SoLIDGEMTracker";
  assert( name && parent );
  
  
  try {
#ifdef MCDATA
     if( dynamic_cast<SoLIDTrackerSystem*>(GetMainDetector())->TestBit(SoLIDTrackerSystem::kMCData) )
       fHits = new TClonesArray("SoLIDMCGEMHit", OUTMAX);
     else
#endif
       fHits = new TClonesArray("SoLIDGEMHit", OUTMAX);
     }
     catch( std::bad_alloc ) {
     Error( Here(here), "Out of memory allocating hit array for tracker "
            " %d. Call expert.", fTrackerID );
     MakeZombie();
   }
}
//_________________________________________________________________________________________
SoLIDGEMTracker::~SoLIDGEMTracker()
{

  if( fIsSetup )
    RemoveVariables();
    
  DeleteContainer(fGEMChamber);
  delete fHits;
}
//_________________________________________________________________________________________
void SoLIDGEMTracker::Clear( Option_t* opt)
{
 if( !opt or *opt != 'I' ) {
    for (Int_t i=0; i<fNChamber; ++i){
      fGEMChamber[i]->Clear(opt);
    }
  }
  fHits->Clear(opt);
  fNHits = 0;
}
//_________________________________________________________________________________________
Int_t SoLIDGEMTracker::Decode( const THaEvData& evdata)
{
  //cout<<"GEM tracker "<<fTrackerID<<endl;
  for (Int_t i=0; i<fNChamber; i++){
    fGEMChamber[i]->Decode(evdata);
  }
  //if (!fDoCombineHits) return 0; //fNHits = CombineChamberHits();
  
  fNHits = 0;
  for (Int_t i=0; i<fNChamber; i++){
    TSeqCollection* hits = fGEMChamber[i]->GetHits();
    fNHits +=  hits->GetEntries();
  }
  

#ifdef MCDATA
  Int_t nHit = 0;
  Int_t countMC = 0;
  
  fMCDecoder = dynamic_cast<const Podd::SimDecoder*>(&evdata);
  
  if (dynamic_cast<SoLIDTrackerSystem*>(GetMainDetector())->TestBit(SoLIDTrackerSystem::kMCData)){
    for (Int_t i = 0; i < fMCDecoder->GetNMCPoints(); i++){
      Podd::MCTrackPoint* pt = fMCDecoder->GetMCPoint(i);
      if (pt->fPlane != fTrackerID) continue; //not from this tracker
      if (pt->fType != 0) continue;
      countMC++;
      Double_t mcX = (pt->fMCPoint).X();
      Double_t mcY = (pt->fMCPoint).Y();
      
      //cout<<fTrackerID<<" "<< (pt->fMCPoint).X()<<" "<<(pt->fMCPoint).Y()<<" "<<(pt->fMCPoint).Z()<<" "<<pt->fType<<endl; 
      //find which chamber the MC hit belongs to
      Double_t phi = atan2( mcY, mcX );
      Int_t chamberList[3] = {-1, -1, -1};
      //for PVDIS there is no need to search for neighboring chamber
      Int_t searchChamber = 3 < fNChamber ? 3 : 1 ;
      
      for (Int_t j=0; j<fNChamber; j++){
        Double_t dphi = phi - fGEMChamber[j]->GetPhiInLab();
        dphi = TVector2::Phi_mpi_pi(dphi);
        if ( dphi < fGEMChamber[j]->GetPhiCover()/2. && dphi > -1*fGEMChamber[j]->GetPhiCover()/2.){
          for (Int_t k=0; k<3; k++){
            chamberList[k] = j-1+k;
            if (chamberList[k] < 0) chamberList[k] += fNChamber;
            if (chamberList[k] > fNChamber-1) chamberList[k] -= fNChamber;
          }
          break;
        }
      }
        
      SoLIDMCGEMHit *closeHit = NULL;
      Double_t dist = 1.e9;
        
      for (Int_t j=0; j<searchChamber; j++){
          TSeqCollection* hits = fGEMChamber[chamberList[j]]->GetHits();
          Double_t mcU, mcV;
          if (!fGEMChamber[chamberList[j]]->CartesianToUV(mcX, mcY, mcU, mcV) ) continue;
          
          for (Int_t k=0; k<hits->GetLast()+1; k++){
            SoLIDMCGEMHit *thisHit = (SoLIDMCGEMHit*)hits->At(k); 
            
            if ( fabs(mcU - thisHit->GetUPos()) > 3.* fGEMChamber[chamberList[j]]->GetReadOut(0)->GetPitch() || 
                 fabs(mcV - thisHit->GetVPos()) > 3.* fGEMChamber[chamberList[j]]->GetReadOut(1)->GetPitch() )
                 continue;
            
            Double_t thisDist = sqrt( pow(mcX - thisHit->GetX(), 2) 
                                   +  pow(mcY - thisHit->GetY(), 2) );
            if ( thisDist < dist){
                dist = thisDist;
                closeHit = thisHit;
            }
          }
      }
      if (closeHit != NULL && closeHit->IsSignalHit()) closeHit->SetGoodMCHit(1);
      if (!fDoCombineHits && closeHit != NULL){
        new ( (*fHits)[nHit++]) SoLIDMCGEMHit(*closeHit);
      }
        
    }
  }
  if (!fDoCombineHits) assert(nHit <= countMC); //for each MC hit, there is at most one cloest hit
#endif
  if (fDoCombineHits) CombineChamberHits();
  return 1;
}
//_________________________________________________________________________________________
THaAnalysisObject::EStatus SoLIDGEMTracker::Init( const TDatime& date )
{
  cout<<"Initializing tracker "<<fTrackerID<<endl;
  EStatus status = THaAnalysisObject::Init(date);
  if (status == kOK){
    //fGEMChamber = new SoLIDGEMChamber *[fNChamber];
    for (Int_t i=0; i<fNChamber; i++){
      stringstream sn, sd;
      sn <<i;
      sd << "Chamber " << i << " on GEM Tracker "<<fTrackerID;
      SoLIDGEMChamber *theChamber = new SoLIDGEMChamber(i, sn.str().c_str(), 
                                         sd.str().c_str(), this);
      fGEMChamber.push_back(theChamber);
      status = fGEMChamber[i]->Init(date);
      if (status) break;
    }
  }
  
  if( status ){
  return fStatus = status;
  }
  return fStatus = kOK;
}
//_________________________________________________________________________________________
void SoLIDGEMTracker::Print( Option_t* opt ) const
{

}
//_________________________________________________________________________________________
void SoLIDGEMTracker::PrintDataBase(Int_t level) const
{
  if (level>=0 && level<3){
    if (level == 0){
      
      int parent_systemID = dynamic_cast<SoLIDTrackerSystem*>(GetParent())->GetSystemID();
      string out_prefix = Form("solid.trackersystem.%d.", parent_systemID);
      cout<<"******parameter from database for tracker "<<fTrackerID<<" in tracker system "<<
      parent_systemID<<"******"<<endl;
      cout<<out_prefix<<fTrackerID<<".nchamber = "<<fNChamber<<endl;
      cout<<out_prefix<<fTrackerID<<".tracker_z = "<<fTrackerZ<<endl;
      cout<<"************************************************************************"<<endl;
    }else if(level > 0){
      level--;
      for (Int_t i=0; i<fNChamber; i++){
        fGEMChamber[i]->PrintDataBase(level);
      }
    }
  }
} 
//_________________________________________________________________________________________
Int_t SoLIDGEMTracker::Begin( THaRunBase* /*r*/ )
{
  return 0;
}
//_________________________________________________________________________________________
Int_t SoLIDGEMTracker::End( THaRunBase* /*r*/ )
{
  return 0;
}
//_________________________________________________________________________________________
Int_t SoLIDGEMTracker::ReadDatabase( const TDatime& date )
{
  static const char* const here = "SoLIDGEMTracker::ReadDatabase";
  fIsInit = kFALSE;
  FILE* file = OpenFile( date );
  if( !file ) return kFileError;
  
  fNChamber = -1;
  fTrackerZ = -100;
  fDoCombineHits = -1;
  const DBRequest request[] = {
        { "nchamber",          &fNChamber,        kInt,     0, 0 },
        { "tracker_z",         &fTrackerZ,        kDouble,  0, 0 },
        { "combine_hits",      &fDoCombineHits,   kInt,     0, 0 },
        { "kill_cross_talk",   &fKillCrossTalk,   kInt,     0, 0 },
        { "cross_talk_thres",  &fCrossTalkThres,  kDouble,  0, 0 },
        { "cross_strip_apart", &fCrossStripApart, kInt,     0, 0 },
        { 0 }
      };
  Int_t status = LoadDB( file, date, request, fPrefix );
  assert(fNChamber > 0 && fTrackerZ != -100 && fDoCombineHits>=0);
  
  fIsInit = kTRUE;
  return kOK;
}
//__________________________________________________________________________________________
Int_t SoLIDGEMTracker::DefineVariables( EMode mode )
{
  // initialize global variables
  
  if( mode == kDefine && fIsSetup ) return kOK;
  fIsSetup = ( mode == kDefine );
  // Register variables in global list
  Int_t ret;

  if( !dynamic_cast<SoLIDTrackerSystem*>(GetMainDetector())->TestBit(SoLIDTrackerSystem::kMCData) ) {
    RVarDef vars[] = {
      { "nhits",       "total number of hits on this tracker", "GetNHits()"        },
      { 0 },
    };
  ret = DefineVarsFromList( vars, mode );
  }else{
#ifdef MCDATA
    //Monte-Carlo hit data
    RVarDef mcvars[] = {
      { "nhits",       "total number of hits on this tracker", "GetNHits()"        },
      { "hit2D.x",       "2D hit x coordinate",            "fHits.SoLIDMCGEMHit.GetX()"        },
      { "hit2D.y",       "2D hit y coordinate",            "fHits.SoLIDMCGEMHit.GetY()"        },
      { "hit2D.z",       "2D hit z coordinate",            "fHits.SoLIDMCGEMHit.GetZ()"        },
      { "hit2D.r",       "2D hit r coordinate",            "fHits.SoLIDMCGEMHit.GetR()"        },
      { "hit2D.qu",      "2D hit charge deposition on u",  "fHits.SoLIDMCGEMHit.GetQU()"       },
      { "hit2D.qv",      "2D hit charge deposition on v",  "fHits.SoLIDMCGEMHit.GetQV()"       },
      { "hit2D.upos",    "2D hit position on u",           "fHits.SoLIDMCGEMHit.GetUPos()"     },
      { "hit2D.vpos",    "2D hit position on v",           "fHits.SoLIDMCGEMHit.GetVPos()"     },
      { "hit2D.usize",   "2D hit cluster size on u",      "fHits.SoLIDMCGEMHit.GetUSize()"  },
      { "hit2D.vsize",   "2D hit cluster size on v",      "fHits.SoLIDMCGEMHit.GetVSize()"  },
      { "hit2D.umcpos",  "2D hit mc position on u",        "fHits.SoLIDMCGEMHit.GetUPosMC()"   },
      { "hit2D.vmcpos",  "2D hit mc position on v",        "fHits.SoLIDMCGEMHit.GetVPosMC()"   },
      { "hit2D.phi",     "2D hit phi coordinate",          "fHits.SoLIDMCGEMHit.GetPhi()"      },
      { "hit2D.signal",  "if this hit is signal",          "fHits.SoLIDMCGEMHit.IsSignalHit()" },
      { "hit2D.chamber", "2D hit chamber ID",              "fHits.SoLIDMCGEMHit.GetChamberID()"  },
      { "hit2D.isgoodmc","if it is the closest to MC hit", "fHits.SoLIDMCGEMHit.IsGoodMCHit()" },
      { "u.occu",        "occupancy pass noise cut",       "GetUMeanOccu()"                   },
      { "v.occu",        "occupancy pass noise cut",       "GetVMeanOccu()"                   },
      { "u.rawoccu",    "deconvoluted occupancy for u",   "GetUMeanHitOccu()"                },
      { "v.rawoccu",    "deconvoluted occupancy for v",   "GetVMeanHitOccu()"                },
      { 0 }
    };
    ret = DefineVarsFromList( mcvars, mode );
#endif
  }
  
  return ret;
}
//___________________________________________________________________________________________
Int_t SoLIDGEMTracker::CombineChamberHits()
{
  //collect all the MC hits info here on this GEM tracker

  assert(fDoCombineHits == 1);
  Int_t nTotalHits = 0;
  Int_t nHit = 0;
  for (Int_t i=0; i<fNChamber; i++){
    TSeqCollection* hits = fGEMChamber[i]->GetHits();
    nTotalHits += hits->GetEntries();
    for (Int_t j=0; j<fGEMChamber[i]->GetNHits(); j++){
#ifdef MCDATA
       SoLIDMCGEMHit *thisHit = (SoLIDMCGEMHit*)hits->At(j);
       //if (!thisHit->IsSignalHit()) continue;
       new ( (*fHits)[nHit++]) SoLIDMCGEMHit(*thisHit);
#endif
    }
    
  }
  return nTotalHits;
}
//_____________________________________________________________________________________________
Double_t SoLIDGEMTracker::GetUMeanOccu()
{
  Double_t meanUOccu = 0.;
  for (Int_t i=0; i<fNChamber; i++){
    meanUOccu += fGEMChamber[i]->GetUOccupancy(); 
  }
  return meanUOccu/(Double_t)fNChamber;
}
//_____________________________________________________________________________________________
Double_t SoLIDGEMTracker::GetVMeanOccu()
{
  Double_t meanVOccu = 0.;
  for (Int_t i=0; i<fNChamber; i++){
    meanVOccu += fGEMChamber[i]->GetVOccupancy(); 
  }
  return meanVOccu/(Double_t)fNChamber;
}
//______________________________________________________________________________________________
Double_t SoLIDGEMTracker::GetUMeanHitOccu()
{
  Double_t meanUOccu = 0.;
  for (Int_t i=0; i<fNChamber; i++){
    meanUOccu += fGEMChamber[i]->GetUHitOccupancy(); 
  }
  return meanUOccu/(Double_t)fNChamber;
}
//______________________________________________________________________________________________
Double_t SoLIDGEMTracker::GetVMeanHitOccu()
{
  Double_t meanVOccu = 0.;
  for (Int_t i=0; i<fNChamber; i++){
    meanVOccu += fGEMChamber[i]->GetVHitOccupancy(); 
  }
  return meanVOccu/(Double_t)fNChamber;
}













