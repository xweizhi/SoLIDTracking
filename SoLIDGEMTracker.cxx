//c++
#include <cassert>
#include <iostream>
#include <sstream>
#include <string>
//SoLID Tracking
#include "SoLIDGEMTracker.h"
#include "SoLIDTrackerSystem.h" 
#include "SoLIDUtility.h"

#define OUTMAX 4096
ClassImp(SoLIDGEMTracker)

using namespace std;
SoLIDGEMTracker::SoLIDGEMTracker(Int_t iGEM, const char* name, const char* description,
                                 THaDetectorBase* parent)
  : THaSubDetector(name,description,parent), fTrackerID(iGEM)
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
}
//_________________________________________________________________________________________
Int_t SoLIDGEMTracker::Decode( const THaEvData& evdata)
{
  for (Int_t i=0; i<fNChamber; i++){
    fGEMChamber[i]->Decode(evdata);
  }
  if (fDoCombineHits) CombineChamberHits();
  return 1;
}
//_________________________________________________________________________________________
THaAnalysisObject::EStatus SoLIDGEMTracker::Init( const TDatime& date )
{
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
Int_t SoLIDGEMTracker::Begin( THaRunBase* r )
{
  return 0;
}
//_________________________________________________________________________________________
Int_t SoLIDGEMTracker::End( THaRunBase* r )
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
  cout<<"Initializing GEM tracker "<<fTrackerID<<endl;
  const DBRequest request[] = {
        { "nchamber",       &fNChamber,        kInt,     0, 1 },
        { "tracker_z",      &fTrackerZ,        kDouble,  0, 1 },
        { "combine_hits",   &fDoCombineHits,   kInt,     0, 1 },
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
#ifdef MCDATA//FIXME:probably a bug
  if( !dynamic_cast<SoLIDTrackerSystem*>(GetMainDetector())->TestBit(SoLIDTrackerSystem::kMCData) ) {
#endif
     // Non-Monte Carlo hit data
    RVarDef nonmcvars[] = {
      { "hit2D.x",       "2D hit x coordinate",          "fHits.SoLIDGEMHit.GetX()"         },
      { "hit2D.y",       "2D hit y coordinate",          "fHits.SoLIDGEMHit.GetY()"         },
      { "hit2D.z",       "2D hit z coordinate",          "fHits.SoLIDGEMHit.GetZ()"         },
      { "hit2D.r",       "2D hit r coordinate",          "fHits.SoLIDGEMHit.GetR()"         },
      { "hit2D.phi",     "2D hit phi coordinate",        "fHits.SoLIDGEMHit.GetPhi()"       },
      { "hit2D.chamber", "2D hit chamber ID",            "fHits.SoLIDGEMHit.GetChamberID()" },
      { 0 }   
    };
    ret = DefineVarsFromList( nonmcvars, mode );
  }else{
    //Monte-Carlo hit data
    RVarDef mcvars[] = {
      { "hit2D.x",       "2D hit x coordinate",            "fHits.SoLIDMCGEMHit.GetX()"        },
      { "hit2D.y",       "2D hit y coordinate",            "fHits.SoLIDMCGEMHit.GetY()"        },
      { "hit2D.z",       "2D hit z coordinate",            "fHits.SoLIDMCGEMHit.GetZ()"        },
      { "hit2D.r",       "2D hit r coordinate",            "fHits.SoLIDMCGEMHit.GetR()"        },
      { "hit2D.qu",      "2D hit charge deposition on u",  "fHits.SoLIDMCGEMHit.GetQU()"       },
      { "hit2D.qv",      "2D hit charge deposition on v",  "fHits.SoLIDMCGEMHit.GetQV()"       },
      { "hit2D.phi",     "2D hit phi coordinate",          "fHits.SoLIDMCGEMHit.GetPhi()"      },
      { "hit2D.signal",  "if this hit is signal",          "fHits.SoLIDMCGEMHit.IsSignalHit()" },
      { "hit2D.chamber", "2D hit chamber ID",              "fHits.SoLIDGEMHit.GetChamberID()"  },
      { 0 }
    };
    ret = DefineVarsFromList( mcvars, mode );
  }
  
  return ret;
}
//___________________________________________________________________________________________
Int_t SoLIDGEMTracker::CombineChamberHits()
{
  //first check if the total hits on all chamber is too large
  assert(fDoCombineHits == 1);
  Int_t totalHit = 0;
  for (Int_t i=0; i<fNChamber; i++){
    totalHit += fGEMChamber[i]->GetNHits();
  }
  if (totalHit>=OUTMAX) cout<<"too many hits on tracker: "<<fTrackerID<<" ("<<totalHit<<")"<<endl;
  if (totalHit>=OUTMAX) return 0;
  Bool_t mc_data = dynamic_cast<SoLIDTrackerSystem*>(GetMainDetector())->TestBit(SoLIDTrackerSystem::kMCData);
  Int_t nHit = 0;
  for (Int_t i=0; i<fNChamber; i++){
    TSeqCollection* hits = fGEMChamber[i]->GetHits();
    for (Int_t j=0; j<fGEMChamber[i]->GetNHits(); j++){
       if (!mc_data){
         SoLIDGEMHit *thisHit = (SoLIDGEMHit*)hits->At(j);
         new ( (*fHits)[nHit++]) SoLIDGEMHit(*thisHit);
       }else{
#ifdef MCDATA
         SoLIDMCGEMHit *thisHit = (SoLIDMCGEMHit*)hits->At(j);
         new ( (*fHits)[nHit++]) SoLIDMCGEMHit(*thisHit);
#endif
       }
    }
    
  }
  
  assert(nHit == totalHit);
  return 1;
}














