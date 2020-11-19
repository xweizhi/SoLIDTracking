//c++
#include <cassert>
#include <iostream>
#include <sstream>
//ROOT
#include "TMath.h"
//Hall A Analyzer
#ifndef MCDATA
#include "THaEvData.h"
#else
#include "SimDecoder.h"
#endif
//SoLID Tracking
#include "SoLIDGEMReadOut.h" 
#include "SoLIDGEMChamber.h"
#include "SoLIDGEMTracker.h"
#include "SoLIDTrackerSystem.h"
#include "SoLIDGEMHit.h"
#define kBig 1e38
ClassImp(SoLIDGEMReadOut)

using namespace std;
using namespace Podd;
static const Int_t kMaxNChan = 20000;

SoLIDGEMReadOut::SoLIDGEMReadOut(Int_t ireadout, const char* name, const char* description,
                                 THaDetectorBase* parent)
  : THaSubDetector(name,description,parent), fReadOutID(ireadout), fADCraw(0), fADC(0), 
  fHitTime(0), fADCcor(0), fGoodHit(0), fDNoise(0), fStripsSeen(0), 
  fNRawStrips(0), fNHitStrips(0), fHitOcc(0), fOccupancy(0), fHits(0), fMapType(kOneToOne), 
  fADCMap(0)
#ifdef MCDATA
  , fMCHitInfo(0), fHitMap(0)
#endif

{
  static const char* const here = "SoLIDGEMReadOut";
  assert( name && parent ); 
  assert( dynamic_cast<SoLIDTrackerSystem*>(GetMainDetector()) );  

  try {
#ifdef MCDATA
    if( dynamic_cast<SoLIDTrackerSystem*>(GetMainDetector())->TestBit(SoLIDTrackerSystem::kMCData) )
      fHits = new TClonesArray("SoLIDMCRawHit", 200);
    else
#endif
      fHits = new TClonesArray("SoLIDRawHit", 200);
  }
  catch( std::bad_alloc ) {
    Error( Here(here), "Out of memory allocating hit array for readout "
           "plane %s. Call expert.", name );
    MakeZombie();
    return;
  }
  fSigStrips[0].clear();
  fSigStrips[1].clear();

}
//_________________________________________________________________________________________
SoLIDGEMReadOut::~SoLIDGEMReadOut()
{
  // Destructor.
  // Histograms in fHist should be automatically deleted by ROOT when
  // the output file is closed

  if( fIsSetup ){
    RemoveVariables();
  }
  // fHits deleted in base class
  delete fGoodHit;
  delete fADCcor;
  delete fHitTime;
  delete fADC;
  delete fADCraw;

  delete fHits;
}
//_________________________________________________________________________________________
void SoLIDGEMReadOut::Clear( Option_t* opt)
{
  assert( fIsInit );
  assert( fADCraw and fADC and fHitTime and fADCcor and fGoodHit );
  memset( fADCraw, 0, fNChan*sizeof(Float_t) );
  memset( fADC, 0, fNChan*sizeof(Float_t) );
  memset( fHitTime, 0, fNChan*sizeof(Float_t) );
  memset( fADCcor, 0, fNChan*sizeof(Float_t) );
  memset( fGoodHit, 0, fNChan*sizeof(Byte_t) );
  fSigStrips[0].clear();
  fSigStrips[1].clear();
  fStripsSeen.assign( fNChan, false );
  
  fNHitStrips = fNRawStrips = 0;
  fHitOcc = fOccupancy = fDNoise = 0.0;
  
   // Clear event-by-event data (hits)

  fHits->Clear(opt);
#ifdef MCDATA
  if( fMCHitInfo ) {
    assert( dynamic_cast<SoLIDTrackerSystem*>(GetMainDetector())->TestBit(SoLIDTrackerSystem::kMCData) );
    for( Vint_t::size_type i = 0; i < fMCHitList.size(); ++i ) {
      fMCHitInfo[i].MCClear();
    }
    fMCHitList.clear();
  }
#endif


}
//_________________________________________________________________________________________
Int_t SoLIDGEMReadOut::Decode( const THaEvData& evdata)
{
   // Extract this plane's hit data from the raw evdata.
  //
  // This routine decodes the front-end readout data.
  // Finds clusters of active strips (=above software threshold) and
  // computes weighted average of position. Each such cluster makes one "Hit".
  
  const char* const here = "SoLIDGEMReadOut::Decode";

#ifdef MCDATA
  
  bool mc_data = dynamic_cast<SoLIDTrackerSystem*>(GetMainDetector())
                 ->TestBit(SoLIDTrackerSystem::kMCData);
  assert( !mc_data || dynamic_cast<const SimDecoder*>(&evdata) != 0 );
#endif
  assert( fADCraw and fADC and fADCcor and fHitTime );

#ifdef TESTCODE
  if( fDoHisto )
    assert( fHitMap != 0 and fADCMap != 0 );
#endif
  assert( fPed.empty() or
          fPed.size() == static_cast<Vflt_t::size_type>(fNChan) );
  assert( fSigStrips[0].empty() );
  assert( fSigStrips[1].empty() );
  assert( fStripsSeen.size() == static_cast<Vbool_t::size_type>(fNChan) );

  UInt_t nHits = 0;

  // Set up pedestal and noise corrections
  Double_t noisesum = 0.0;
  UInt_t   n_noise = 0;

  bool do_pedestal_subtraction = !fPed.empty();
  bool do_noise_subtraction    = fDoNoise;

#ifdef MCDATA
  const SimDecoder* simdata = 0;
  if( mc_data ) {
    assert( dynamic_cast<const SimDecoder*>(&evdata) );
    simdata = static_cast<const SimDecoder*>(&evdata);
  }
#endif
  Vflt_t samples;
  if( fMaxTimeSample > 1 )
    samples.reserve(fMaxTimeSample);

   // Decode data
  for( Int_t imod = 0; imod < fDetMap->GetSize(); ++imod ) {
    THaDetMap::Module * d = fDetMap->GetModule(imod);

    // Read the active channels
    Int_t nchan = evdata.GetNumChan( d->crate, d->slot );
    for( Int_t ichan = 0; ichan < nchan; ++ichan ) {
      Int_t chan = evdata.GetNextChan( d->crate, d->slot, ichan );
      if( chan < d->lo or chan > d->hi ) continue; // not part of this detector

      // Map channel number to strip number
      Int_t istrip =
	    MapChannel( d->first + ((d->reverse) ? d->hi - chan : chan - d->lo) );
      // Test for duplicate istrip, if found, warn and skip
      assert( istrip >= 0 and istrip < fNChan );
      if( fStripsSeen[istrip] ) {
      const char* inp_source = "DAQ";
#ifdef MCDATA
	    if( mc_data )
	    inp_source = "digitization";
#endif
	Warning( Here(here), "Duplicate strip number %d in plane %s, event %d. "
		 "Ignorning it. Fix your %s.",
		 istrip, GetName(), evdata.GetEvNum(), inp_source );
	continue;
      }
      fStripsSeen[istrip] = true;

      // For the APV25 analog pipeline, multiple "hits" on a decoder channel
      // correspond to time samples 25 ns apart
      Int_t nsamp = evdata.GetNumHits( d->crate, d->slot, chan );
      assert( nsamp > 0 );
      ++fNRawStrips;
      nsamp = TMath::Min( nsamp, static_cast<Int_t>(fMaxTimeSample) );
      if (fDeconMode == 1) nsamp = 5;
      // Integrate the signal over time and analyze pulse shape
      StripData_t stripdata;
      if( nsamp > 1 ) {
	samples.clear();
	for( Int_t isamp = 0; isamp < nsamp; ++isamp ) {
	  Float_t fsamp = static_cast<Float_t>
	    ( evdata.GetData(d->crate, d->slot, chan, isamp) );
          //if (istrip >= fNStrip) fsamp = 0.; 
	  samples.push_back( fsamp );
	}
        //cout<<istrip<<" "<<samples[2]<<" "<<fReadOutID<<endl;
	// Analyze the pulse shape
	if (fDeconMode == 1){
            stripdata = ChipChargeDep(samples);
        }else if (fDeconMode == 2){
            stripdata = SAMPAChargeDep(samples);
        }
        else if (fDeconMode == 3){
            stripdata = VMMChargeDep(samples);
        }
        else{
	    stripdata = ChargeDep(samples);
        }
      }
      else {
        stripdata.adcraw = stripdata.adc =
          static_cast<Float_t>( evdata.GetData(d->crate, d->slot, chan, 2) );
        stripdata.time = 0;
        stripdata.pass = true;
      }
      // Skip null data
      if( stripdata.adcraw == 0 )
	continue;

      // Save results for cluster finding later
      fADCraw[istrip]  = stripdata.adcraw;
      fADC[istrip]     = stripdata.adc;
      fHitTime[istrip] = stripdata.time;

      // Strip-by-strip pedestal subtraction
      Float_t adc = stripdata.adc;
      if( do_pedestal_subtraction )
	adc -= fPed[istrip];

      fADCcor[istrip] = adc;
      fGoodHit[istrip] = not fCheckPulseShape or stripdata.pass;

      if( do_noise_subtraction ) {
	// Sum up ADCs that are likely not a hit
	if( adc < fADCMin ) {
	  noisesum += adc;
	  n_noise++;
	}
      }
      else {
	// If no noise subtraction is done, then we can finish up with this
	// strip number right here. Otherwise we need a second iteration below
	AddStrip( istrip );
      }

#ifdef MCDATA
      // If doing MC data, save the truth information for each strip
      if( mc_data ) {
	fMCHitList.push_back(istrip);
	fMCHitInfo[istrip] = simdata->GetMCHitInfo(d->crate,d->slot,chan);
      }
#endif
    }  // chans
  }    // modules
  
   // Calculate average common-mode noise and subtract it from corrected
  // ADC values, if requested
  if( do_noise_subtraction ) {
    if ( n_noise > 0 ) {
      fDNoise = noisesum/n_noise;
      assert( fDNoise < fADCMin );
    }
    // Save strip numbers of corrected ADC data above threshold. Fill histograms.
    assert( fSigStrips[0].empty() );
    assert( fSigStrips[1].empty() );
    for( Int_t i = 0; i < fNChan; i++ ) {
      fADCcor[i] -= fDNoise;
      AddStrip( i );
    }
  }

  fHitOcc    = static_cast<Double_t>(fNHitStrips) / fNChan;
  fOccupancy = static_cast<Double_t>(GetNSigStrips()) / fNChan;
  
   // Find and analyze clusters. Clusters of active strips are considered
  // a "Hit".
  //
  // The cluster analysis is a critical part of the GEM analysis. Various
  // things can and probably need to be done right here already: splitting
  // oversized clusters, detecting noise hits/bogus clusters, detecting and
  // fitting overlapping clusters etc.
  //
  // This analysis may even need to be re-done after preliminary tracking to
  // see if the clustering can be improved using candidate tracks.
  // Additionally, correlated amplitude information from a second readout
  // direction in the same readout plane could be used here. These advanced
  // procedures would require significant redesign of the code:
  // all raw strip info will have to be saved and prcessed at a later point,
  // similar to the finding of hit pairs in like-oriented planes of the MWDC.
  //
  // For the moment, we implement a very simple algorithm: any cluster of
  // strips larger than what a single cluster should be is assumed to be two or
  // more overlapping hits, and the cluster will be split as follows: anything
  // that looks like a local peak followed by a valley will be considered an
  // actual cluster. The parameter frac = fSplitFrac (0.0 ... 1.0) can
  // be used for some crude tuning. frac > 0.0 means that a peak is
  // only a peak if the amplitude drops below (1-frac), so
  // frac = 0.1 means: trigger on a drop below 90% etc. Likewise for the
  // following valley: the bottom is found if the amplitude rises again
  // by (1+frac), so frac = 0.1 means: trigger on a rise above 110% etc.

  // The active strip numbers must be sorted for the clustering algorithm
  if (fSigStrips[0].size() == 0 && fSigStrips[1].size() == 0) return 1;
  sort( ALL(fSigStrips[0]) );
  sort( ALL(fSigStrips[1]) );
  typedef Vint_t::iterator viter_t;
  Double_t frac_down = 1.0 - fSplitFrac, frac_up = 1.0 + fSplitFrac;

  /*cout<<"print strip info"<<endl;
  for (int sc = 0; sc < 2; sc++){
    for (unsigned int i=0; i<fSigStrips[sc].size(); i++)
    cout<<sc<<" "<<fSigStrips[sc][i]<<" "<<fADCcor[fSigStrips[sc][i]]<<endl;
  }
  cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl; 
  */

  for (Int_t sc = 0; sc<2; sc++){
      if (fSigStrips[sc].size() == 0) continue;

      #ifndef NDEBUG
      SoLIDRawHit* prevHit = 0;
    #endif
      Vint_t splits;  // Strips with ampl split between 2 clusters
      viter_t next = fSigStrips[sc].begin();
      while( next != fSigStrips[sc].end() ) {
        viter_t start = next, cur = next;
        ++next;

        assert( next == fSigStrips[sc].end() or *next > *cur );
        while( next != fSigStrips[sc].end() and (*next - *cur) == 1  ) {
          ++cur;
          ++next;
        }
        // Now the cluster candidate is between start and cur
        assert( *cur >= *start );
        // The "type" parameter indicates the result of the cluster analysis:
        // 0: clean (i.e. smaller than fMaxClusterSize, no further analysis)
        // 1: large, maximum at right edge, not split
        // 2: large, no clear minimum on the right side found, not split
        // 3: split, well-defined peak found (may still be larger than maxsize)
        Int_t  type = 0;
        UInt_t size = *cur - *start + 1;
        if( size > (UInt_t)fMaxClusterSize ) {
          Double_t maxadc = 0.0, minadc = kBig;
          viter_t it = start, maxpos = start, minpos = start;
          enum EStep { kFindMax = 1, kFindMin, kDone };
          EStep step = kFindMax;
          while( step != kDone and it != next ) {
            Double_t adc = fADCcor[*it];
            switch( step ) {
              case kFindMax:
                // Looking for maximum
                if( adc > maxadc ) {
                  maxpos = it;
                  maxadc = adc;
                } else if( adc < maxadc * frac_down) {
                  assert( maxadc > 0.0 );
                  step = kFindMin;
                  continue;
                }
                break;
              case kFindMin:
                // Looking for minimum
                if( adc < minadc ) {
                  minpos = it;
                  minadc = adc;
                } else if( adc > minadc * frac_up) {
                  assert( minadc < kBig );
                  step = kDone;
                }
                break;
              case kDone:
                assert( false );  // should never get here
                break;
            }
            ++it;
          }
          if( step == kDone ) {
            // Found maximum followed by minimum
            assert( minpos != start );
            assert( minpos != cur );
            assert( *minpos > *maxpos );
            // Split the cluster at the position of the minimum, assuming that
            // the strip with the minimum amplitude is shared between both clusters
            cur  = minpos;
            next = minpos;
            // In order not to double-count amplitude, we split the signal height
            // of that strip evenly between the two clusters. This is a very
            // crude way of doing what we really should be doing: "fitting" a peak
            // shape and using the area and centroid of the curve
	    fADCcor[*minpos] /= 2.0;
	    splits.push_back(*minpos);
          }
          type = step;
          size = *cur - *start + 1;
          assert( *cur >= *start );
        }
        assert( size > 0 );
        // Compute weighted position average. Again, a crude (but fast) substitute
        // for fitting the centroid of the peak.
        Double_t xsum = 0.0, adcsum = 0.0;
        Short_t  maxID = -1;
        Double_t maxSave = 1e-9;
    #ifdef MCDATA
        Double_t mcpos = 0.0, mctime = kBig;
        Int_t mctrack = 0, num_bg = 0;
    #endif
        for( ; start != next; ++start ) {
          Int_t istrip = *start;
          Double_t pos = GetChanPos(istrip); //GetStart() + istrip * GetPitch();
          Double_t adc = fADCcor[istrip];
          //cout<<istrip<<" "<<adc<<" "<<sc<<endl;
          xsum   += pos * adc;
          adcsum += adc;
          if (adc > maxSave){
              maxSave = adc;
              maxID   = istrip;
          }
    #ifdef MCDATA
          // If doing MC data, analyze the strip truth information
          if( mc_data ) {
	    MCHitInfo& mc = fMCHitInfo[istrip];
	    // This may be smaller than the actual total number of background hits
	    // contributing to the entire cluster, but counting them would involve
	    // lists of secondary particle numbers ... overkill for now
	    num_bg = TMath::Max( num_bg, mc.fContam );
	    // All primary particle hits in the cluster are from the same track
	    //assert( mctrack == 0 || mc.fMCTrack == 0 || mctrack == mc.fMCTrack );

        //fMCTrack > 0 means signal hit, the value is actually the track ID
        //fMCTrack == 0 means background hit
        //fMCTrack == -1 means it is induced hit on the readout strip

	    if( mctrack == 0 ) {
	      if( mc.fMCTrack > 0 ) {
	        // If the cluster contains a signal hit, save its info and be done
	        mctrack = mc.fMCTrack;
	        mcpos   = mc.fMCPos;
	        mctime  = mc.fMCTime;
	      }
	      else if (mc.fMCTrack == 0){
	        // If background hits only, compute position average
	        mcpos  += mc.fMCPos;
	        mctime  = TMath::Min( mctime, mc.fMCTime );
	      }else{
            assert(mc.fMCTrack == -1); //induced strip, cannot be anything else
            mcpos  += pos; //no other avaiable information for induced strip
            //no timing information right now for induced strip as well
          }
	    }
          }
    #endif // MCDATA
        }
        assert( adcsum > 0.0 );
        Double_t pos = xsum/adcsum;
        //cout<<pos<<" "<<sc<<endl;
    #ifdef MCDATA
        if( mc_data && mctrack == 0 ) {
          mcpos /= static_cast<Double_t>(size);
        }
    #endif
        // The resolution (sigma) of the position measurement depends on the
        // cluster size. In particular, if the cluster consists of only a single
        // hit, the resolution is much reduced
        Double_t resolution = fResolution;
        if (type == 3) resolution *= 2.;
        if( size == 1 ) {
          resolution = TMath::Max( 0.25*GetPitch(), fResolution );
          // The factor of 1/2*pitch is just a guess. Since with real GEMs
          // there _should_ always be more than one strip per cluster, we must
          // assume that the other strip(s) did not fire due to inefficiency.
          // As a result, the error is bigger than it would be if only ever one
          // strip fired per hit.
    //       resolution = TMath::Max( 0.5*GetPitch(), 2.0*fResolution );
    //     } else if( size == 2 ) {
    //       // Again, this is a guess, to be quantified with Monte Carlo
    //       resolution = 1.2*fResolution;
        }

    Bool_t splitType = sc == 0 ? true : false;
        // Make a new hit
    #ifndef NDEBUG
        SoLIDRawHit* theHit = 0;
    #endif
    #ifdef MCDATA
        if( !mc_data ) {
    #endif
    #ifndef NDEBUG
          theHit =
    #endif
	    new( (*fHits)[nHits++] ) SoLIDRawHit( pos,
					     adcsum,
					     size,
					     type,
					     resolution,
					     splitType,
                                             maxID,
					     this
					     );
    #ifdef MCDATA
        } else {
          // Monte Carlo data
    #ifndef NDEBUG
          theHit =
    #endif
	    new( (*fHits)[nHits++] ) SoLIDMCRawHit( pos,
					       adcsum,
					       size,
					       type,
					       resolution,
					       splitType,
					       maxID,
					       this,
					       mctrack,
					       mcpos,
					       mctime,
					       num_bg
					       );
					       //cout<<mctrack<<" "<<mcpos<<" "<<pos<<" "<<sc<<endl;
        }
    #endif // MCDATA
    #ifndef NDEBUG
        // Ensure hits are ordered by position (should be guaranteed by std::map)
        assert( prevHit == 0 or theHit->Compare(prevHit) > 0 );
        prevHit = theHit;
    #endif
      }
  }
  if (fKillCrossTalk) KillCrossTalk();

  return 1;
}
//_________________________________________________________________________________________
void SoLIDGEMReadOut::KillCrossTalk()
{
    //to identify potential cluster that is induced by a larger neighboring signal
    //This induction happens when the signal goes into the APV chip
    //The order of channels inside APV chips does not correspond to order of channel on 
    //the readout board. So that the induced signal appear on a channel that is about 32 strips
    //away from the strip that the main signal is on. And the amplitude is about 10% of 
    //the main signal
    
    //array should be ordered by position (in increaseing order)
    
    Int_t nhit = GetNhits();
    for (Int_t i=0; i<nhit; i++){
        Int_t rIndex = i >= nhit -1 ? i : i+1;
        Int_t lIndex = i == 0 ? 0 : i-1;
        
        SoLIDRawHit* cHit = (SoLIDRawHit*)fHits->At(i);
        Int_t cid = (cHit->GetPos() - GetStart()) / GetPitch();
        //look cluster on the right side of the current cluster
        while (rIndex < nhit){
            SoLIDRawHit* rHit = (SoLIDRawHit*)fHits->At(rIndex);
            Int_t did = (rHit->GetPos() - GetStart()) / GetPitch() - cid;
            assert(did >= 0);
            
            if (did > fCrossStripApart) { break; }
            else if ( did == fCrossStripApart ) {
                    if (rHit->GetADCsum() < fCrossTalkThres*cHit->GetADCsum())
                    { rHit->Deactivate(); }
            }
            rIndex++;
        }
        //do the same for left side
        while (lIndex >= 0){
            SoLIDRawHit* lHit = (SoLIDRawHit*)fHits->At(lIndex);
            Int_t did = (lHit->GetPos() - GetStart()) / GetPitch() - cid;
            assert(did <= 0);
            
            if (did < -1*fCrossStripApart) { break; }
            else if ( did == -1*fCrossStripApart ) {
                    if (lHit->GetADCsum() < fCrossTalkThres*cHit->GetADCsum())
                    { lHit->Deactivate(); }
            }
            lIndex--;
        }
    }
}
//_________________________________________________________________________________________
THaAnalysisObject::EStatus SoLIDGEMReadOut::Init( const TDatime& date )
{
  EStatus status = THaAnalysisObject::Init(date);
  if( status ){
  return fStatus = status;
  }
  return fStatus = kOK;
}
//_________________________________________________________________________________________
void SoLIDGEMReadOut::Print( Option_t* /*opt*/ ) const
{

}
//_________________________________________________________________________________________
void SoLIDGEMReadOut::PrintDataBase(Int_t level) const
{
  if (level == 0){ //this is the bottom level
    //oh god I know
    SoLIDGEMChamber *parent = dynamic_cast<SoLIDGEMChamber*>(GetParent());
    SoLIDGEMTracker *grand_parent = dynamic_cast<SoLIDGEMTracker*>(parent->GetParent());
    SoLIDTrackerSystem *great_grand_parent = dynamic_cast<SoLIDTrackerSystem*>(grand_parent->GetParent());
    Int_t great_grand_parent_systemID = great_grand_parent->GetSystemID();
    Int_t grand_parent_trackerID = grand_parent->GetTrackerID();
    Int_t parent_chamberID = parent->GetChamberID();
    
    string out_prefix = Form("solid.trackersystem.%d.%d.%d.%d", great_grand_parent_systemID, 
                            grand_parent_trackerID, parent_chamberID, fReadOutID);
    cout<<"******parameter from database for readout "<<fReadOutID<<" of chamber "
      <<parent_chamberID<<" in tracker "<<grand_parent_trackerID<<" in tracker system"<<
      great_grand_parent_systemID<<"******"<<endl;
    cout<<out_prefix<<".nstrips = "<<fNStrip<<endl;
    cout<<out_prefix<<".strip_angle = "<<fStripAngle<<endl;
    cout<<out_prefix<<".xp_res = "<<fResolution<<endl;
    cout<<out_prefix<<".maxclustsiz = "<<fMaxClusterSize<<endl;
    cout<<out_prefix<<".adc_min = "<<fADCMin<<endl;
    cout<<out_prefix<<".split_frac = "<<fSplitFrac<<endl;
    cout<<out_prefix<<".maxhits = "<<fMaxHits<<endl;
    cout<<out_prefix<<".maxsamp = "<<fMaxTimeSample<<endl;
    cout<<out_prefix<<".adc_sigma = "<<fADCSigma<<endl;
    cout<<out_prefix<<".do_noise = "<<fDoNoise<<endl;
    cout<<out_prefix<<".check_pulse_shape = "<<fCheckPulseShape<<endl;
    cout<<out_prefix<<".do_histos = "<<fDoHisto<<endl;
    cout<<out_prefix<<".strip_pitch = "<<fStripPitch<<endl;
    cout<<out_prefix<<".dz = "<<fDz<<endl;
    cout<<out_prefix<<".dphi = "<<fDPhi<<endl;
    cout<<out_prefix<<".start = "<<fStartPos<<endl;
    for( Int_t imod = 0; imod < fDetMap->GetSize(); ++imod ) {
      THaDetMap::Module* d = fDetMap->GetModule(imod);
      cout<<out_prefix<<".detmap = "<<d->crate<<" "<<d->slot<<" "<<d->lo<<" "<<d->hi<<endl;
    }
  }
} 
//_________________________________________________________________________________________
Int_t SoLIDGEMReadOut::Begin( THaRunBase* /*r*/ )
{
  return 0;
}
//_________________________________________________________________________________________
Int_t SoLIDGEMReadOut::End( THaRunBase* /*r*/ )
{
  return 0;
}
//_________________________________________________________________________________________
Int_t SoLIDGEMReadOut::ReadDatabase( const TDatime& date )
{
  static const char* const here = "SoLIDGEMReadOut::ReadDataBase";
  fIsInit = kFALSE;
  
  FILE* file = OpenFile( date );
  if( !file ) return kFileError;
  Int_t status = -1;
  fNStrip = 0;
  fNDivStrip = 0;
  fDivStart = -1;
  fStripAngle  = 361.;
  fResolution = -1.;
  fMaxClusterSize = -1;
  fADCMin = -1;
  fSplitFrac = -1.;
  fMaxHits = -1;
  fMaxTimeSample = -1;
  fADCSigma = 0.36; // default, an educated guess
  fDoNoise = -1;
  fCheckPulseShape = -1;
  fDoHisto = -1;
  fStripPitch = -1.;
  fDz = -100.;
  fDPhi = -100.;
  fStartPos = -100.;
  fPed.clear();
  fChanMap.clear();
  TString mapping;
  vector<Int_t>* detmap = 0;
  
  //TODO: some of these parameter can go to the tracker database, lots of
  //repititions here
  
  try{
    detmap = new vector<Int_t>;
    const DBRequest request[] = {
          { "detmap",               detmap,           kIntV },
          { "nstrips",             &fNStrip,          kInt,     0, 0 },
          { "strip_angle",         &fStripAngle,      kDouble,  0, 0 },
          { "xp_res",              &fResolution,      kDouble,  0, 0 },
          { "maxclustsiz",         &fMaxClusterSize,  kInt,     0, 0 },
          { "adc_min",             &fADCMin,          kInt,     0, 0 },
          { "split_frac",          &fSplitFrac,       kDouble,  0, 0 },  
          { "maxhits",             &fMaxHits,         kInt,     0, 0 },
          { "maxsamp",             &fMaxTimeSample,   kInt,     0, 0 },
          { "adc_sigma",           &fADCSigma,        kDouble,  0, 0 },
          { "do_noise",            &fDoNoise,         kInt,     0, 0 },
          { "check_pulse_shape",   &fCheckPulseShape, kInt,     0, 0 },
          { "do_histos",           &fDoHisto,         kInt,     0, 0 },
          { "strip_pitch",         &fStripPitch,      kDouble,  0, 0 },
          { "dz",                  &fDz,              kDouble,  0, 0 },
          { "dphi",                &fDPhi,            kDouble,  0, 0 },
          { "start",               &fStartPos,        kDouble,  0, 0 },
          { "nstrips",             &fNStrip,          kInt,     0, 0 },
          { "n_divstrips",         &fNDivStrip,       kInt,     0, 0 },
          { "div_start",           &fDivStart,        kInt,     0, 0 },
          { "deconmode",           &fDeconMode,       kInt,     0, 0 },
          { "mapping",             &mapping,          kTString, 0, 1 }, // not using it right now
          { "chanmap",             &fChanMap,         kIntV,    0, 1 }, // not using it right now
          { "pedestal",            &fPed,             kFloatV,  0, 1 }, // not using it right now
          { 0 }
        };

    status = LoadDB( file, date, request, fPrefix );
  
    assert(fNStrip>0 && fResolution>0 && fMaxClusterSize>0 && fADCMin>=0);
    assert(fSplitFrac>0 && fMaxHits>0 && fMaxTimeSample>0 && fADCSigma>=0);
    assert(fDoNoise>=0 && fCheckPulseShape>=0 && fDoHisto>=0 && fStripPitch>0);
    assert(fStripAngle>=0. && fStripAngle<=360.);
    assert(fStartPos != -100. && fDz != -100. && fDPhi != -100.);
  
    if (status == kOK){
      // Parse the detector map of the data channels
      if( FillDetMap( *detmap, /*use default for now*/0, here ) <= 0 ) status = kInitError;
    }
    delete detmap;
  }catch(...) {
      delete detmap;
      fclose(file);
      throw;
  }
   // Finished reading the database
  fclose(file);
  if( status != kOK ) return status;
  
  fNChan = fNStrip +fNDivStrip;
  
  //change it from deg to rad here once and for all
  fStripAngle *= TMath::DegToRad();
  fSinStripAngle = TMath::Sin(fStripAngle);
  fCosStripAngle = TMath::Cos(fStripAngle);
  
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
  
  if( fNChan <= 0 or fNChan > kMaxNChan ) {
    Error( Here(here), "Invalid number of channels: %d. Must be > 0 and < %d. "
           "Fix database or recompile code with a larger limit.",
           fNChan, kMaxNChan );
    return kInitError;
  }

  Int_t nchan = fDetMap->GetTotNumChan();
  if( nchan != fNChan ) {
    Error( Here(here), "Number of detector map channels (%d) "
           "disagrees with number of strips (%d)", nchan, fNStrip );
    return kInitError;
  }
  
  SafeDelete(fADCraw);
  SafeDelete(fADC);
  SafeDelete(fHitTime);
  SafeDelete(fADCcor);
  SafeDelete(fGoodHit);
  
#ifdef MCDATA
  delete [] fMCHitInfo; fMCHitInfo = 0;
#endif

  // Allocate arrays. The only reason that these are parallel C-arrays is
  // that the global variable system still doesn't support arrays/vectors
  // of structures/objects.
  // Out of memory exceptions from here are caught in Tracker.cxx.
  fADCraw = new Float_t[fNChan];
  fADC = new Float_t[fNChan];
  fHitTime = new Float_t[fNChan];
  fADCcor = new Float_t[fNChan];
  fGoodHit = new Byte_t[fNChan];
  fSigStrips[0].reserve(fNStrip);
  fSigStrips[1].reserve(fNDivStrip);
  fStripsSeen.resize(fNChan);

#ifdef MCDATA
  assert( dynamic_cast<SoLIDTrackerSystem*>(GetMainDetector()) );
  if(dynamic_cast<SoLIDTrackerSystem*>( GetMainDetector() )->TestBit(SoLIDTrackerSystem::kMCData)){
    fMCHitInfo = new Podd::MCHitInfo[fNChan];
    fMCHitList.reserve(fNChan);
  }
#endif

  TString::ECaseCompare cmp = TString::kIgnoreCase;
  if( !mapping.IsNull() ) {
  cout<<"not using mapping for now, should not get here"<<endl;
    if( mapping.Length() >= 3 and
        TString("one-to-one").BeginsWith(mapping,cmp) )
      fMapType = kOneToOne;
    else if( mapping.Length() >=3 and
             TString("reverse").BeginsWith(mapping,cmp) )
      fMapType = kReverse;
    else if( mapping.Length() >= 5 and
             mapping.BeginsWith(TString("gassiplex-adapter"),cmp) ) {
      if( fNChan > 240 ) {
        Error( Here(here), "Gassiplex adapter mapping allows at most 240 "
               "strips, but %d configured. Fix database.", fNChan );
        return kInitError;
      }
      if( fNChan < 240 ) {
        Warning( Here(here), "Gassiplex adapter mapping expects 240 "
                 "strips, but %d configured. Database may be misconfigured "
                 "(or you know what you are doing).", fNChan );
      }
      if( mapping.BeginsWith(TString("gassiplex-adapter-2"),cmp) ) {
        fMapType = kGassiplexAdapter2;
      } else {
        fMapType = kGassiplexAdapter1;
      }
    }
    else if( TString("table").CompareTo(mapping,cmp) ) {
      if( fChanMap.empty() ) {
        Error( Here(here), "Channel mapping table requested, but no map "
               "defined. Specify chanmap in database." );
        return kInitError;
      }
      if( fChanMap.size() != static_cast<UInt_t>(fNChan) ) {
        Error( Here(here), "Number of channel map entries (%u) must equal "
               "number of strips (%d). Fix database.",
               static_cast<unsigned int>(fChanMap.size()), fNChan );
        return kInitError;
      }
      // check if entries in channel map are within range
      for( Vint_t::const_iterator it = fChanMap.begin();
           it != fChanMap.end(); ++it ) {
        if( (*it) < 0 or (*it) >= fNChan ) {
          Error( Here(here), "Illegal chanmap entry: %d. Must be >= 0 and "
                 "< %d. Fix database.", (*it), fNChan );
          return kInitError;
        }
      }
      fMapType = kTable;
    } else {
      Error( Here(here), "Unknown channel mapping type %s. Fix database.",
             mapping.Data() );
      return kInitError;
    }

  } else
    fChanMap.clear();

  if( !fPed.empty() and fPed.size() != static_cast<UInt_t>(fNChan) ) {
  cout<<"not using pedestal yet, should not get here"<<endl;
    Error( Here(here), "Size of pedestal array (%u) must equal "
           "number of strips (%d). Fix database.",
           static_cast<unsigned int>(fPed.size()), fNChan );
    return kInitError;
  }

   // Sanity checks on fMaxTimeSample
  static const UInt_t max_maxsamp = 32; // arbitrary sanity limit
  if( fMaxTimeSample == 0 )
    fMaxTimeSample = 1;
  else if( fMaxTimeSample > (Int_t)max_maxsamp ) {
    Warning( Here(here), "Illegal maximum number of samples: %u. "
             "Adjusted to maximum allowed = %u.", fMaxTimeSample, max_maxsamp );
    fMaxTimeSample = max_maxsamp;
  }
  if( fMaxTimeSample == 1 )
    fCheckPulseShape = 0;

  if( fADCSigma < 0.0 ) {
    Warning( Here(here), "Negative adc.sigma = %lf makes no sense. Adjusted "
             "to positive.", fADCSigma );
    fADCSigma = TMath::Abs( fADCSigma );
  }
  if( fADCSigma < 0.001 ) {
    Warning( Here(here), "adc.sigma = %lf is extremely small. "
             "Double-check database.", fADCSigma );
  }

  UpdateOffset();
  
  SoLIDGEMChamber *parent = dynamic_cast<SoLIDGEMChamber*>(GetParent());
  SoLIDGEMTracker *grand_parent = dynamic_cast<SoLIDGEMTracker*>(parent->GetParent());
  
  fKillCrossTalk = grand_parent->KillCrossTalk();
  fCrossTalkThres = grand_parent->GetCrossTalkThres();
  fCrossStripApart = grand_parent->GetCrossStripApart();

  fIsInit = true;
  return kOK;

}
//__________________________________________________________________________________________
Int_t SoLIDGEMReadOut::DefineVariables( EMode mode )
{
   // initialize global variables
  if( mode == kDefine && fIsSetup ) return kOK;
  fIsSetup = ( mode == kDefine );

  // Register variables in global list

   RVarDef vars[] = {
    { "nrawstrips",     "nstrips with decoder data",        "fNRawStrips" },
    { "nhitstrips",     "nstrips > 0",                      "fNHitStrips" },
    { "nstrips",        "Num strips with hits > adc_min",   "GetNSigStrips()" },
    { "hitocc",         "strips > 0 / n_all_strips",        "fHitOcc" },
    { "occupancy",      "nstrips / n_all_strips",           "fOccupancy" },
    { "strip.adcraw",   "Raw strip ADC sum",                "fADCraw" },
    { "strip.adc",      "Deconvoluted strip ADC sum",       "fADC" },
    { "strip.adc_c",    "Pedestal-sub strip ADC sum",       "fADCcor" },
    { "strip.time",     "Leading time of strip signal (ns)","fHitTime" },
    { "strip.good",     "Good pulse shape on strip",        "fGoodHit" },
    { "nhits",          "Num hits (clusters of strips)",    "GetNhits()" },
    { "noise",          "Noise level (avg below adc.min)",  "fDNoise" },
    { 0 }
  };
  Int_t ret = DefineVarsFromList( vars, mode );

  if( ret != kOK )
    return ret;

  #ifdef MCDATA
  if( !dynamic_cast<SoLIDTrackerSystem*>(GetMainDetector())->TestBit(SoLIDTrackerSystem::kMCData) ) {
#endif
    // Non-Monte Carlo hit data
    RVarDef nonmcvars[] = {
      { "hit.pos",  "Hit centroid (m)",      "fHits.SoLIDRawHit.fPos" },
      { "hit.adc",  "Hit ADC sum",           "fHits.SoLIDRawHit.fADCsum" },
      { "hit.size", "Num strips ",           "fHits.SoLIDRawHit.fSize" },
      { "hit.type", "Hit analysis result",   "fHits.SoLIDRawHit.fType" },
      { 0 }
    };
    ret = DefineVarsFromList( nonmcvars, mode );
#ifdef MCDATA
  } else {
    // Monte Carlo hit data includes the truth information
    // For safety, we make sure that all hit variables are referenced with
    // respect to the MCGEMHit class and not just GEMHit - the memory layout
    // of classes under multiple inheritance might be implemetation-dependent
    RVarDef mcvars[] = {
      { "hit.pos",   "Hit centroid (m)",      "fHits.SoLIDMCRawHit.fPos" },
      { "hit.adc",   "Hit ADC sum",           "fHits.SoLIDMCRawHit.fADCsum" },
      { "hit.size",  "Num strips ",           "fHits.SoLIDMCRawHit.fSize" },
      { "hit.type",  "Hit analysis result",   "fHits.SoLIDMCRawHit.fType" },
      { "hit.mctrk", "MC track number",       "fHits.SoLIDMCRawHit.fMCTrack" },
      { "hit.mcpos", "MC track position (m)", "fHits.SoLIDMCRawHit.fMCPos" },
      { "hit.mctime","MC track time (s)",     "fHits.SoLIDMCRawHit.fMCTime" },
      { "hit.numbg", "MC num backgr hits",    "fHits.SoLIDMCRawHit.fContam" },
      { 0 }
    };
    ret = DefineVarsFromList( mcvars, mode );
  }
#endif
  return ret;

}
//__________________________________________________________________________________________
StripData_t SoLIDGEMReadOut::ChargeDep( const vector<Float_t>& amp )
{
   // Deconvolute signal given by samples in 'amp', return approximate integral.
  // Currently analyzes exactly 3 samples.
  // From Kalyan Allada
  // NIM A326, 112 (1993)

  //FIXME: from database, proper value for Tp
  const Float_t delta_t = 25.0; // time interval between samples (ns)
  const Float_t Tp      = 56.0; // RC filter time constant (ns)

  assert( amp.size() >= 3 );

  //Float_t adcsum;

  //adcsum = amp[0]+amp[1]+amp[2];


  // Weight factors calculated based on the response of the silicon microstrip
  // detector:
  // v(t) = (delta_t/Tp)*exp(-delta_t/Tp)
  // Need to update this for GEM detector response(?):
  // v(t) = A*(1-exp(-(t-t0)/tau1))*exp(-(t-t0)/tau2)
  // where A is the amplitude, t0 the begin of the rise, tau1 the time
  // parameter for the rising edge and tau2 the for the falling edge.

  Float_t x = delta_t/Tp;

  Float_t w1 = TMath::Exp(x-1)/x;
  Float_t w2 = -2*TMath::Exp(-1)/x;
  Float_t w3 = TMath::Exp(-x-1)/x;

  // Deconvoluted signal samples, assuming measurements of zero before the
  // leading edge
  Float_t sig[3] = { amp[0]*w1,
                     amp[1]*w1+amp[0]*w2,
                     amp[2]*w1+amp[1]*w2+amp[0]*w3 };
  Float_t adc;

  //this adc value may not be very reliable, we assume that the time sample befre
  //the trigger start time is 0, but for background, this is not true, and this
  //deconvolution algorithm can give rather strange behavior

  adc    = (sig[0]+sig[1]+sig[2]);

  //adc = amp[2];
  Float_t max = amp[2];
  if (amp[0] > amp[2] && amp[0] > amp[1]) max = amp[0];
  if (amp[1] > amp[0] && amp[1] > amp[2]) max = amp[1];

  Float_t time   = 0;     // TODO

  Bool_t pass;
  // Calculate ratios for 3 samples and check for bad signals
#ifdef PVDIS
  if( amp[2] > 0 ) {
    Double_t r1 = amp[0] < fADCSigma ? 0 : amp[0]/(amp[2]);
    Double_t r2 = amp[1] < fADCSigma ? 0 : amp[1]/(amp[2]);
    pass = (r1 < 1.1 and r2 < 1.3 and (r1) < r2);
    } else
     pass = false;
#endif

#ifdef SIDIS //actually for JPsi
  adc = max;
  pass = true;
  if (amp[0] > amp[1] && amp[0] > amp[2] && amp[1] > amp[2]) pass = false;
#endif
  //pass = true;
  //if (amp[0] > amp[1] && amp[1] > amp[2]) pass = false;
  return StripData_t(max, adc, time, pass);
}
//_____________________________________________________________________________________
StripData_t SoLIDGEMReadOut::ChipChargeDep( const vector<Float_t>& amp )
{
  const Float_t delta_t = 25.0; // time interval between samples (ns)
  const Float_t Tp      = 56.0; // RC filter time constant (ns)

  assert( amp.size() >= 5 );
  
  vector<Float_t> s;
  for (UInt_t i=0; i<amp.size(); i++){
    if (amp[i] < 3.*fADCSigma) s.push_back(0);
    else s.push_back(amp[i]);
  }

  //Float_t adcraw = s[3];

  Float_t x = delta_t/Tp;

  Float_t w1 = TMath::Exp(x-1)/x;
  Float_t w2 = -2*TMath::Exp(-1)/x;
  Float_t w3 = TMath::Exp(-x-1)/x;

  /* Float_t sig[3] = { amp[2]*w1 + amp[1]*w2 + amp[0]*w3,
                      amp[3]*w1 + amp[2]*w2 + amp[1]*w3,
                      amp[4]*w1 + amp[3]*w2 + amp[2]*w3   };*/
  Float_t sig[3] = { amp[1]*w1,
                     amp[2]*w1+amp[1]*w2,
                     amp[3]*w1+amp[2]*w2+amp[1]*w3 };
  
  Float_t adc = sig[2];
  
  Float_t time   = 0;     // TODO

  Bool_t pass = true;

  return StripData_t( adc, adc,time,pass );
}
//______________________________________________________________________________________
StripData_t SoLIDGEMReadOut::SAMPAChargeDep( const vector<Float_t>& amp )
{
    assert( amp.size() >= 5 );
    
    Float_t adc = 0.;
    for (unsigned int i=1; i<7; i++) adc += amp[i];
    adc /= 6.;
    
    Bool_t pass = false;
    
     if (  ((amp[3] > amp[2] && amp[3] > amp[4]) || 
            (amp[4] > amp[3] && amp[4] > amp[5]) || 
            (amp[5] > amp[4] && amp[5] > amp[6]) ||
            (amp[3] > 950 && amp[4] > 950 && amp[5] > 950) ) ) pass = true;
            
     if (! (amp[3] > amp[1] || amp[4] > amp[1] || amp[5] > amp[1])) pass = false;
    
    Float_t time = 0.;
    
    return StripData_t(adc, adc, time, pass);
}
//______________________________________________________________________________________
StripData_t SoLIDGEMReadOut::VMMChargeDep( const vector<Float_t>& amp )
{
    assert( amp.size() >= 1 );
    Float_t adc = amp[0];
    Bool_t pass = true; 
    Float_t time = 0.;
    return StripData_t(adc, adc, time, pass);
}
//______________________________________________________________________________________
Int_t SoLIDGEMReadOut::MapChannel( Int_t idx ) const
{
  // Map hardware channel number to logical strip number based on mapping
  // prescription from database

  assert( idx >= 0 and idx < fNChan );

  Int_t ret = 0;
  switch( fMapType ) {
  case kOneToOne:
    ret = idx;
    break;
  case kReverse:
    ret = fNChan-idx-1;
    break;
  case kGassiplexAdapter1:
    assert( idx < 240 );
    if( idx == 0 )
      ret = 1;
    else if( idx == 239 )
      ret = 238;
    else if( idx % 2 ) // odd
      ret = idx + 2;
    else               // even
      ret = idx - 2;
    break;
  case kGassiplexAdapter2:
    assert( idx < 240 );
    if( idx == 1 )
      ret = 0;
    else if( idx == 238 )
      ret = 239;
    else if( idx % 2 ) // odd
      ret = idx - 2;
    else               // even
      ret = idx + 2;
    break;
  case kTable:
    // Use the mapping lookup table
    assert( fChanMap.size() == static_cast<Vint_t::size_type>(fNChan) );
    ret = fChanMap[idx];
    break;
  }
  assert( ret >= 0 and ret < fNChan );
  return ret;
}
//_____________________________________________________________________________
void SoLIDGEMReadOut::AddStrip( Int_t istrip )
{
  // Record a hit on the given strip number in internal arrays.
  // Utility function used by Decode.

  Float_t adc = fADCraw[istrip];
  if( adc >= fADCMin ) ++fNHitStrips;
    
  if( fGoodHit[istrip] and adc >= fADCMin ) {
    if (istrip < fNStrip) fSigStrips[0].push_back(istrip);
    else fSigStrips[1].push_back(istrip);
  }
#ifdef TESTCODE
  if( fDoHisto ) {
    fHitMap->Fill(istrip);
    fADCMap->Fill(istrip, adc);
  }
#endif
}
//_____________________________________________________________________________
void SoLIDGEMReadOut::UpdateOffset()
{
  //by adding this offset when doing the hit clustering, later on when we combine
  //the 2 hit coordinates, the hit should already in the lab frame.

  //For PVDIS the reference point get from the parent chamber class should be (0,0)
  //because in the PVDIS configuration, if you extend the two lines that parallel to
  //the edge of the chamber, they should have an intersection at (0,0), and this means
  //that all the geometric calculation in the ReadGeometry function in the upper class
  //meke sense
  //
  //For the SIDIS, however, it is more involved. The Chambers will be move inward in
  //the radial direction, which means the two lines no longer pass (0,0), but somewhere
  //on the other side, say (a,b). All the calculation in the ReadGeometry function still
  //meke sense, but with respect to this (a,b).


  fCoorOffset = 0;
  TVector3 trackerCenter = dynamic_cast<SoLIDGEMChamber*>(GetParent())->GetCenterPos();

  fCoorOffset = trackerCenter.X()*TMath::Cos(fStripAngle) +
                trackerCenter.Y()*TMath::Sin(fStripAngle);
  cout<<"offset "<<fCoorOffset<<endl;
}
//_______________________________________________________________________________
Double_t SoLIDGEMReadOut::GetChanPos(Int_t ichan)
{
    if (ichan < fNStrip)
    return GetStart() + ichan * GetPitch();
    else return GetStart() + (ichan - fNStrip + fDivStart)*GetPitch();
}
//________________________________________________________________________________
Int_t SoLIDGEMReadOut::GetBaseStripID(Int_t ichan)
{
    if (ichan < fNStrip) return ichan;
    else return ichan - fNStrip + fDivStart;
}
//________________________________________________________________________________
bool SoLIDGEMReadOut::IsStripDivided(Int_t ichan)
{
    if (fNDivStrip == 0) return false;
    int baseID = GetBaseStripID(ichan);
    if (baseID >= fDivStart && baseID < fDivStart + fNDivStrip) return true;
    else return false;
}












