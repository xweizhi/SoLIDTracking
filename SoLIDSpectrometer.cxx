//*-- Author :    Ole Hansen  03-December-2012

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// SoLID::SoLIDSpectrometer                                             //
//                                                                      //
// SoLID spectrometer class containing GEM trackers                     //
// Used for testing.                                                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
//c++
#include <string>
#include <sstream>
#include <exception>
#include <cassert>
//root
#include "TList.h"
//Hall A analyzer
#include "THaGlobals.h"
#include "THaTextvars.h"
#include "THaTrack.h"
//SoLID tracking
#include "SoLIDSpectrometer.h"
#include "EProjType.h"

using namespace std;

ClassImp(SoLIDSpectrometer)

//_____________________________________________________________________________
SoLIDSpectrometer::SoLIDSpectrometer( const char* name, const char* description,
		  UInt_t nsystem, UInt_t ntracker, UInt_t nsector, UInt_t nreadout)
  : THaSpectrometer( name, description ), fNSystem(nsystem), fNTracker(ntracker),
    fNSector(nsector), fNReadOut(nreadout)
{
  // Constructor. Define a GEM tracker detector for each of the 'nsectors'
  // sectors

 static const char* const here = "SoLIDSpectrometer";
 cout<<"create SoLIDSpectrometer"<<endl;

  if( nsector > 64 ) {
    Error( Here(here), "Number of sectors = %u too large. Must be <= 64. "
	   "Creating SoLID spectrometer failed.", nsector );
    MakeZombie();
    throw range_error("SoLIDSpectrometer: number of sectors too large");
  }
  if( nsector == 0 ) {
    Warning( Here(here), "Creating SoLID spectrometer with zero sectors" );
  }

  //For SIDIS, we have 6 GEM detector and currently assume 30 sectors each
  //we need to consider hits in these 30 sectors all at once
  //For PVDIS, WE HAVE 5 GEM detector nad currently assume 30 sectors each
  //However, each of them works as individual system, so effectively only 
  //1 sector
  
  for (UInt_t isystem=0; isystem<nsystem; ++isystem){
    stringstream sn, sd;
    sn << "trackersystem." << isystem;
    sd << "SoLID track system " << isystem;
    SoLIDTrackerSystem *theSystem = new SoLIDTrackerSystem(sn.str().c_str(), sd.str().c_str());
	  theSystem->SetSystemID(isystem);
	  fTrackerSystem.push_back(theSystem);
    Int_t ret = AddDetector( fTrackerSystem.at(isystem) );
    if( ret != 0) {
      stringstream s;
      s << "Error adding GEM track system \"" << sn.str() << "\" "
	      << "to SoLID spectrometer \"" << GetName() << "\"";
      MakeZombie();
      throw logic_error(s.str());
    }
  }

  fSolTrackInfo = new TClonesArray( "SolTrackInfo", kInitTrackMultiplicity );

  // For now, don't require run database
  // We might need it later to read things like the field setting
  fProperties &= ~kNeedsRunDB;
}

//_____________________________________________________________________________
SoLIDSpectrometer::~SoLIDSpectrometer()
{
  // Destructor

  delete fSolTrackInfo; fSolTrackInfo = 0;
  DefineVariables( kDelete );
  DeleteContainer(fTrackerSystem);
}

//_____________________________________________________________________________
void SoLIDSpectrometer::Clear( Option_t* opt )
{
  // Clear event-by-event data

  THaSpectrometer::Clear(opt);
  fSolTrackInfo->Clear();
}

//_____________________________________________________________________________
THaAnalysisObject::EStatus SoLIDSpectrometer::Init( const TDatime& run_time )
{
  // Extra initialization for the SoLID spectrometer.
  //
  // To simplify database access, automatically define the following
  // gHaTextvars macros:
  // "allsectors": all sector numbers
  // "plane1": all supported projection types suffixed by "1"
  // "plane2"..."plane5": dto. with suffixes "2"..."5"

  //Int_t nsect = fDetectors->GetSize();
  Int_t nsect = fNSector;
  if( gHaTextvars != 0 ) {
    if( nsect > 0 ) {
      stringstream s;
      for( Int_t i = 0; i < nsect; ++i ) {
	s << i;
	if( i+1 < nsect )
	  s << ",";
      }
      assert( s && !s.str().empty() );
      gHaTextvars->Set( "allchambers", s.str() );
    }
   
    if (fNSystem>0){
     stringstream s;
      for( Int_t i = 0; i < fNSystem; ++i ) {
	      s << i;
	      if( i+1 < fNSystem )
	      s << ",";
      }
      assert( s && !s.str().empty() );
      gHaTextvars->Set( "allsystems", s.str() );
    }
    
    if (fNTracker>0){
     stringstream s;
      for( Int_t i = 0; i < fNTracker; ++i ) {
	      s << i;
	      if( i+1 < fNTracker )
	      s << ",";
      }
      assert( s && !s.str().empty() );
      gHaTextvars->Set( "alltrackers", s.str() );
    }

    if (fNReadOut>0){
     stringstream s;
      for( Int_t i = 0; i < fNReadOut; ++i ) {
	      s << i;
	      if( i+1 < fNReadOut )
	      s << ",";
      }
      assert( s && !s.str().empty() );
      gHaTextvars->Set( "allreadouts", s.str() );
    }
  }

  // Proceed with normal spectrometer initialization
  // EStatus ret = THaSpectrometer::Init( run_time );
  // if( ret != kOK )
  //   return ret;
  return THaSpectrometer::Init( run_time );

  // TODO: set up text variables like "plane1" etc. for output
  // definitions & cuts, using actual plane names defined
  // (quite a job)

  // return kOK;
}

//_____________________________________________________________________________
Int_t SoLIDSpectrometer::DefineVariables( EMode mode )
{
  // Define/delete standard spectrometer variables (tracks etc.) plus
  // SoLID-specific additional track variables (cylindrical coords etc.)
  if( mode == kDefine && fIsSetup ) return kOK;
  // Define standard spectrometer variables (track array)
  if( mode == kDefine )
    THaSpectrometer::DefineVariables(mode);
    
   RVarDef vars[] = {
    { 0 }
   };

  return DefineVarsFromList( vars, mode );

  //return kOK;
}

//_____________________________________________________________________________
Int_t SoLIDSpectrometer::ReadRunDatabase( const TDatime& date )
{
  // Dummy run database reader to override the THaSpectrometer function.

  // This is a bit of a kludge. SoLIDSpectrometer actually shouldn't inherit from
  // THaSpectrometer, but from a different spectrometer base class that
  // doesn't assume small-angle focusing optics.

  return THaApparatus::ReadRunDatabase( date );
}

//_____________________________________________________________________________
Int_t SoLIDSpectrometer::FindVertices( TClonesArray& /*tracks*/ )
{

  return 0;
}

//_____________________________________________________________________________
Int_t SoLIDSpectrometer::TrackCalc()
{
  // Additional track calculations


  return 0;
}

//_____________________________________________________________________________
void SoLIDSpectrometer::PrintDataBase(Int_t level) const
{
  for (Int_t i=0; i<fNSystem; i++){
    fTrackerSystem.at(i)->PrintDataBase(level);
  }
}

