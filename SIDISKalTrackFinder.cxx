//c++
#include <cmath>
//ROOT
#include "TRandom.h"
//SoLIDTracking
#include "SIDISKalTrackFinder.h"
#include "SoLKalTrackSystem.h"
#include "SoLKalTrackSite.h"
#include "SoLKalTrackState.h"

//these should definitely need to go to the database
#define MAXNTRACKS_FAEC 1000
#define MAXNTRACKS_LAEC 1000

ClassImp(SIDISKalTrackFinder)
SIDISKalTrackFinder::SIDISKalTrackFinder(bool isMC, const char* name, int detConf )
:SoLKalTrackFinder(), THaAnalysisObject(name, "SIDIS_Track_Finder"), 
fDetConf(detConf), fIsMC(isMC),fNGoodTrack(0)
{
  Init();
  fGEMTracker.clear();
  fWindowHits.clear();
  fWindowHits.reserve(MAXWINDOWHIT);

  for (int i=0; i<2; i++) {
   fSeedEfficiency[i] = false;
   fMcTrackEfficiency[i] = false;
  }
  fElectronVertexZ = -1e6;
}
//__________________________________________________________________________
SIDISKalTrackFinder::~SIDISKalTrackFinder()
{
  Clear();
  fGoodHits.clear();
  delete fCoarseTracks;
}
//___________________________________________________________________________
void SIDISKalTrackFinder::Clear( Option_t* opt )
{
  
  if (fCoarseTracks->GetEntries() != 0) 
   fCoarseTracks->Delete();
  
  fCoarseTracks->Clear(opt);
  
  fCaloHits = nullptr;
  fNSeeds = 0;
  fNGoodTrack = 0;
 
  for (int i=0; i<2; i++) {
   fSeedEfficiency[i] = false;
   fMcTrackEfficiency[i] = false;
  }

 
  map< Int_t, vector<SoLIDGEMHit*> >::iterator it;
  for (it = fGoodHits.begin(); it != fGoodHits.end(); it++) { (it->second).clear(); }
  fGoodHits.clear();
  
  map< SeedType, vector<DoubletSeed> >::iterator itt;
  for (itt = fSeedPool.begin(); itt != fSeedPool.end(); itt++) { (itt->second).clear(); }
  
  fElectronVertexZ = -1e6;
}
//______________________________________________________________________________
Int_t SIDISKalTrackFinder::ReadDatabase (const TDatime& date)
{
    FILE* file = OpenFile (date);
    if (!file) return kFileError;
    try{
    const DBRequest request[] =
        {
          { "target_center",               &fTargetCenter,                   kDouble, 0, 0},
          { "target_length",               &fTargetLength,                   kDouble, 0, 0},
          { "chi2_per_ndf_cut",            &fChi2PerNDFCut,                  kDouble, 0, 0},
          { "theta_min",                    fThetaMinCut,                    kDouble, 3, 0},
          { "theta_max",                    fThetaMaxCut,                    kDouble, 3, 0},
          { "momentum_min",                 fMomMinCut,                      kDouble, 3, 0},
          { "momentum_max",                 fMomMaxCut,                      kDouble, 3, 0},
          { "cell_edge_cut",                fCellEdgeCut,                    kDouble, 2, 0},
          { "coarse_cell_edge_cut",         fCoarseCellEdgeCut,              kDouble, 3, 0},
          { "backplane_rlimit_min",        &fRlimitMin[0],                   kDouble, 3, 0},
          { "backplane_rlimit_max",        &fRlimitMax[0],                   kDouble, 3, 0},
          { "midplane_rlimit_min",         &fRlimitMin[1],                   kDouble, 3, 0},
          { "midplane_rlimit_max",         &fRlimitMax[1],                   kDouble, 3, 0},
          { "frontplane_rlimit_min",       &fRlimitMin[2],                   kDouble, 3, 0},
          { "frontplane_rlimit_max",       &fRlimitMax[2],                   kDouble, 3, 0},
          { "midback_dr_min",              &fDeltaRMin[0],                   kDouble, 3, 0},
          { "midback_dr_max",              &fDeltaRMax[0],                   kDouble, 3, 0},
          { "frontmid_dr_min",             &fDeltaRMin[1],                   kDouble, 3, 0},
          { "frontmid_dr_max",             &fDeltaRMax[1],                   kDouble, 3, 0},
          { "frontback_dr_min",            &fDeltaRMin[2],                   kDouble, 3, 0},
          { "frontback_dr_max",            &fDeltaRMax[2],                   kDouble, 3, 0},
          { "midback_dphi_min",            &fDeltaPhiMin[0],                 kDouble, 3, 0},
          { "midback_dphi_max",            &fDeltaPhiMax[0],                 kDouble, 3, 0},
          { "frontmid_dphi_min",           &fDeltaPhiMin[1],                 kDouble, 3, 0},
          { "frontmid_dphi_max",           &fDeltaPhiMax[1],                 kDouble, 3, 0},
          { "frontback_dphi_min",          &fDeltaPhiMin[2],                 kDouble, 3, 0},
          { "frontback_dphi_max",          &fDeltaPhiMax[2],                 kDouble, 3, 0},
          { "backec_dr_min",               &fDeltaRMinBC[0],                 kDouble, 6, 0},
          { "backec_dr_max",               &fDeltaRMaxBC[0],                 kDouble, 6, 0},
          { "backec_dphi_min",             &fDeltaPhiMinBC[0],               kDouble, 6, 0},
          { "backec_dphi_max",             &fDeltaPhiMaxBC[0],               kDouble, 6, 0},
          { "midec_dr_min",                &fDeltaRMinBC[1],                 kDouble, 6, 0},
          { "midec_dr_max",                &fDeltaRMaxBC[1],                 kDouble, 6, 0},
          { "midec_dphi_min",              &fDeltaPhiMinBC[1],               kDouble, 6, 0},
          { "midec_dphi_max",              &fDeltaPhiMaxBC[1],               kDouble, 6, 0},
          { "coarse_ec_pos_cut",            fCoarseECPosCut,                 kDouble, 2, 0},
          { "coarse_spd_r_cut",             fCoarseSPDRCut,                  kDouble, 4, 0},
          { "coarse_spd_phi_cut",           fCoarseSPDPhiCut,                kDouble, 4, 0},
          { "ec_pos_cut",                  &fECPosCut,                       kDouble, 0, 0},
          { "ec_energy_match",             &fECEnergyMatch,                  kInt,    3, 0},
          { 0 }
        };
        Int_t err = LoadDB (file, date, request, fPrefix);
        fclose(file);
        if (err)
        return kInitError;
    }  catch(...) {
        fclose(file);
        throw;
    }

    for (int i=0; i<3; i++) cout<<fThetaMinCut[i]<<" ";
    cout<<endl;
    for (int i=0; i<3; i++) cout<<fThetaMaxCut[i]<<" ";
    cout<<endl;
    for (int i=0; i<3; i++) cout<<fThetaMaxCut[i]<<" ";
    cout<<endl;
    for (int i=0; i<3; i++) cout<<fMomMinCut[i]<<" ";
    cout<<endl;
    for (int i=0; i<3; i++) cout<<fMomMaxCut[i]<<" ";
    cout<<endl;

    return kOK;
}
//___________________________________________________________________________
void SIDISKalTrackFinder::ProcessHits(TClonesArray* theTracks)
{
  if (fGEMTracker.size() == 0) return;
  fNSeeds = 0;
  assert(fCaloHits == nullptr);
  fCaloHits = fECal->GetCaloHits();
  //cout<<"starting electron tracking"<<endl;
  //cout<<fCoarseTracks->GetEntries()<<" "<<fCoarseTracks->GetLast()+1<<endl; 
  fHadronFlag = false;
  //forward angle seed finding
  FindDoubletSeed(4, 5, kFAEC);
  FindDoubletSeed(3, 4, kFAEC);
  FindDoubletSeed(3, 5, kFAEC);
  MergeSeed();
  
  map< SeedType, vector<DoubletSeed> >::iterator itt;
  for (itt = fSeedPool.begin(); itt != fSeedPool.end(); itt++) { (itt->second).clear(); }
  
  //large angle seed finding
  FindDoubletSeed(2, 3, kLAEC);
  FindDoubletSeed(1, 2, kLAEC);
  FindDoubletSeed(1, 3, kLAEC);
  MergeSeed();
  
  
  TrackFollow();
  FindandAddVertex();
  ECalFinalMatch();
  FinalSelection(theTracks);
  //cout<<"done electron tracking"<<endl; 
  //cout<<fCoarseTracks->GetEntries()<<" "<<fCoarseTracks->GetLast()+1<<endl;
  
  
  fHadronFlag = true;
  
  if (fCoarseTracks->GetEntries() != 0) 
  fCoarseTracks->Delete();
  fCoarseTracks->Clear();
  for (itt = fSeedPool.begin(); itt != fSeedPool.end(); itt++) { (itt->second).clear(); }
  fNSeeds = 0;
  //cout<<fCoarseTracks->GetEntries()<<" "<<fCoarseTracks->GetLast()+1<<endl;
  

  //forward angle seed finding
  FindDoubletSeed(4, 5, kFAEC);
  FindDoubletSeed(3, 4, kFAEC);
  FindDoubletSeed(3, 5, kFAEC);
  MergeSeed();
  //cout<<fCoarseTracks->GetEntries()<<" "<<fCoarseTracks->GetLast()+1<<endl;
  TrackFollow();
  FindandAddVertex();
  FinalSelection(theTracks);
  
  fEventNum++;
}

//___________________________________________________________________________________________________________________
void SIDISKalTrackFinder::FindDoubletSeed(Int_t planej, Int_t planek, ECType type)
{
  assert(planek > planej);
  
  double philimit[2] = {0};
  double rlimit[2][2] = {{0, 0}, {0, 0}};
  double deltar[2] = {0};
  double dphi, dr;
  SeedType seedType = kMidBack;
  int countSeed = 0;
  
  int ttype = type;
  if (fHadronFlag) ttype += 1;
  
  double charge = 0;
  if (type == kFAEC){
    if (planek == 5 && planej == 4){
      seedType = kMidBack;
      rlimit[1][0] = fRlimitMin[1][ttype]; rlimit[1][1] = fRlimitMax[1][ttype];
      rlimit[0][0] = fRlimitMin[0][ttype]; rlimit[0][1] = fRlimitMax[0][ttype];
      deltar[0]=fDeltaRMin[0][ttype]; deltar[1]=fDeltaRMax[0][ttype]; 
      philimit[0]=fDeltaPhiMin[0][ttype]; philimit[1]=fDeltaPhiMax[0][ttype];
    }
    else if(planek == 4 && planej == 3){
      seedType = kFrontMid;
      rlimit[1][0] = fRlimitMin[2][ttype]; rlimit[1][1] = fRlimitMax[2][ttype]; 
      rlimit[0][0] = fRlimitMin[1][ttype]; rlimit[0][1] = fRlimitMax[1][ttype];
      deltar[0]=fDeltaRMin[1][ttype]; deltar[1]=fDeltaRMax[1][ttype]; 
      philimit[0]=fDeltaPhiMin[1][ttype]; philimit[1]=fDeltaPhiMax[1][ttype];
    }
    else if (planek == 5 && planej == 3){
      seedType = kFrontBack;
      rlimit[1][0] = fRlimitMin[2][ttype]; rlimit[1][1] = fRlimitMax[2][ttype];
      rlimit[0][0] = fRlimitMin[0][ttype]; rlimit[0][1] = fRlimitMax[0][ttype];
      deltar[0]=fDeltaRMin[2][ttype]; deltar[1]=fDeltaRMax[2][ttype]; 
      philimit[0]=fDeltaPhiMin[2][ttype]; philimit[1]=fDeltaPhiMax[2][ttype];
    } 
  }
  else{ 
    if (planek == 3 && planej == 2){
      seedType = kMidBack;
      rlimit[1][0] = fRlimitMin[1][type]; rlimit[1][1] = fRlimitMax[1][type];
      rlimit[0][0] = fRlimitMin[0][type]; rlimit[0][1] = fRlimitMax[0][type];
      deltar[0]=fDeltaRMin[0][type]; deltar[1]=fDeltaRMax[0][type]; 
      philimit[0]=fDeltaPhiMin[0][type]; philimit[1]=fDeltaPhiMax[0][type];
    }
    else if (planek == 2 && planej == 1){
      seedType = kFrontMid;
      rlimit[1][0] = fRlimitMin[2][type]; rlimit[1][1] = fRlimitMax[2][type]; 
      rlimit[0][0] = fRlimitMin[1][type]; rlimit[0][1] = fRlimitMax[1][type];
      deltar[0]=fDeltaRMin[1][type]; deltar[1]=fDeltaRMax[1][type]; 
      philimit[0]=fDeltaPhiMin[1][type]; philimit[1]=fDeltaPhiMax[1][type];
    }
    else if (planek == 3 && planej == 1){
      seedType = kFrontBack;
      rlimit[1][0] = fRlimitMin[2][type]; rlimit[1][1] = fRlimitMax[2][type];
      rlimit[0][0] = fRlimitMin[0][type]; rlimit[0][1] = fRlimitMax[0][type];
      deltar[0]=fDeltaRMin[2][type]; deltar[1]=fDeltaRMax[2][type]; 
      philimit[0]=fDeltaPhiMin[2][type]; philimit[1]=fDeltaPhiMax[2][type];
    }
    
  }
  /*cout<<planek<<" "<<planej<<endl;
  cout<<rlimit[0][0]<<" "<<rlimit[0][1]<<" "<<rlimit[1][0]<<" "<<rlimit[1][1]<<endl;
  cout<<deltar[0]<<" "<<deltar[1]<<endl;
  cout<<philimit[0]<<" "<<philimit[1]<<endl;*/
  
  for (int k=0; k<fGEMTracker[planek]->GetNChamber(); k++){
    TSeqCollection* planekHitArray = fGEMTracker[planek]->GetChamber(k)->GetHits();
  
    for (int nhitk = 0; nhitk < planekHitArray->GetLast()+1; nhitk++){
      SoLIDGEMHit *hitk = (SoLIDGEMHit*)planekHitArray->At(nhitk);

      if (hitk->GetR() < rlimit[0][0]) continue;
      if (hitk->GetR() > rlimit[0][1]) break; // check if the hit is within r range
      
      Double_t matchEnergy = 0.;
      if (!TriggerCheck(hitk, type, matchEnergy)) continue;
      
      vector<Int_t> jChamberList;
      GetHitChamberList(jChamberList, k, 3);
      
      for (int j=0; j<(int)jChamberList.size(); j++){
        TSeqCollection* planejHitArray = fGEMTracker[planej]->GetChamber(jChamberList.at(j))->GetHits();
      
        for (int nhitj = 0; nhitj < planejHitArray->GetLast()+1; nhitj++){
          SoLIDGEMHit *hitj = (SoLIDGEMHit*)planejHitArray->At(nhitj);
          
          //if the number seeds already exceed the limit, terminate the seed finding process
          if (countSeed > MAXNSEEDS) return;

          if (hitj->IsUsed()) continue;
          if (hitj->GetR()<rlimit[1][0]) continue; 
          if (hitj->GetR()>rlimit[1][1]) break;
          
          dr = CalDeltaR(hitk->GetR(), hitj->GetR());
          if (dr > deltar[1]) continue;
          if (dr < deltar[0]) break;
          
          charge = 0;
	        dphi = CalDeltaPhi(hitj->GetPhi(), hitk->GetPhi());
	        if(((dphi >philimit[0]&& dphi <philimit[1])||(dphi < -1*philimit[0]&& dphi > -1*philimit[1]))){
	         
	          if (dphi >philimit[0]&& dphi <philimit[1]){
	              charge = 1;
	          }else{
	              charge = -1;
	          }
          }
          else continue;
          
          assert(charge != 0); //should never happen
          
          
          //using correction function to calculate initial momentum and angles of the particle at plane k
          double initTheta = 0;
          double initMom   = 0;
          double initPhi   = 0;
          if (!CalInitParForPair(hitj, hitk, charge, initMom, initTheta, initPhi, type)) continue;
          
          if (type == kFAEC && (initTheta > fThetaMaxCut[ttype] || initTheta < fThetaMinCut[ttype])) continue;
          if (type == kLAEC && (initTheta > fThetaMaxCut[kLAEC] || initTheta < fThetaMinCut[kLAEC])) continue;
          if (type == kFAEC && (initMom > fMomMaxCut[ttype] || initMom < fMomMinCut[ttype])) continue;
          if (type == kLAEC && (initMom > fMomMaxCut[kLAEC] || initMom < fMomMinCut[kLAEC])) continue;
          
          Double_t e_reso = fECal->GetEReso()/sqrt(matchEnergy);
          if (fECEnergyMatch[ttype] && (matchEnergy - initMom)/matchEnergy > 5.*e_reso) continue;

          
          TVector3 initDir(cos(initPhi), sin(initPhi), 1./tan(initTheta));
          initDir = initDir.Unit();
		      TVector3 initMomentum = initMom*initDir;
		      TVector3 initPosition(hitk->GetX(), hitk->GetY(), hitk->GetZ());
		      TVector3 finalMomentum;
          TVector3 finalPosition;
          Double_t stepSize = 1.;
          
          Bool_t isSeed = false;
          Double_t toZ = fECal->GetECZ(type);
		    if (type == kFAEC){
        
              fFieldStepper->PropagationClassicalRK4(initMomentum, initPosition, toZ, 
                                                 charge, stepSize, finalMomentum, finalPosition);
              for (UInt_t ec_count=0; ec_count<fCaloHits->GetLast()+1; ec_count++){
	            if ( ((SoLIDCaloHit*)fCaloHits->At(ec_count))->fECID != kFAEC) continue; //not FAEC hit
	            Int_t dettype = ((SoLIDCaloHit*)fCaloHits->At(ec_count))->fSubDetID;
	            if (fHadronFlag) {
	              if (dettype == 0) continue;
	            }
	            else{
	              if (dettype != 0) continue;
	            }
	            Double_t projR   = sqrt( pow(finalPosition.X(), 2) + pow(finalPosition.Y(), 2) );
	            Double_t projPhi = atan2(finalPosition.Y(), finalPosition.X());
	            Double_t spdR    = sqrt(pow(((SoLIDCaloHit*)fCaloHits->At(ec_count))->fXPos, 2) 
	                                  + pow(((SoLIDCaloHit*)fCaloHits->At(ec_count))->fYPos, 2));
	            Double_t spdPhi  = atan2(((SoLIDCaloHit*)fCaloHits->At(ec_count))->fYPos, ((SoLIDCaloHit*)fCaloHits->At(ec_count))->fXPos);
	            Double_t rDiff   = fabs(spdR - projR);
	            Double_t phiDiff = fabs(TVector2::Phi_mpi_pi(spdPhi - projPhi));
	            if (!fHadronFlag){
	              if (sqrt( pow(finalPosition.X() - ((SoLIDCaloHit*)fCaloHits->At(ec_count))->fXPos, 2) +  
	                  pow(finalPosition.Y() - ((SoLIDCaloHit*)fCaloHits->At(ec_count))->fYPos, 2) ) < fCoarseECPosCut[type] ) isSeed = true;
	            }
	            else{
	              if (rDiff < fCoarseSPDRCut[dettype-1] && phiDiff < fCoarseSPDPhiCut[dettype-1]) isSeed = true;
	            } 
		      }
		    }
		    else if (type == kLAEC){
              fFieldStepper->PropagationClassicalRK4(initMomentum, initPosition, toZ, 
                                               charge, stepSize, finalMomentum, finalPosition);
              for (UInt_t ec_count=0; ec_count<fCaloHits->GetLast()+1; ec_count++){
	            if ( ((SoLIDCaloHit*)fCaloHits->At(ec_count))->fECID != kLAEC) continue; //not FAEC hit
	            if (sqrt( pow(finalPosition.X() - ((SoLIDCaloHit*)fCaloHits->At(ec_count))->fXPos, 2) +  
	                  pow(finalPosition.Y() - ((SoLIDCaloHit*)fCaloHits->At(ec_count))->fYPos, 2) ) < fCoarseECPosCut[type] ) isSeed = true;
	          } 
		    }
		      if (!isSeed) continue;
          
          
          fFieldStepper->PropagationClassicalRK4(initMomentum, initPosition, 
                                             fTargetCenter, charge, stepSize, finalMomentum, finalPosition);
      
          double tx = finalMomentum.X()/finalMomentum.Z();
          double ty = finalMomentum.Y()/finalMomentum.Z();
		      double ReconZ = fTargetCenter + (1./(pow(tx,2) + pow(ty,2)))*
                       (tx*(fBPMX-finalPosition.X()) + ty*(fBPMY-finalPosition.Y()) );
      
          if (type == kFAEC && fabs(ReconZ - fTargetCenter) > (fTargetLength/2. + fCoarseCellEdgeCut[ttype] ) ) continue;
          if (type == kLAEC && fabs(ReconZ - fTargetCenter) > (fTargetLength/2. + fCoarseCellEdgeCut[kLAEC] ) ) continue;
          
          
          //so the hit pairs has passed all the cuts, now we can save it into a container and waiting for merge
          countSeed++;
          map< SeedType, vector<DoubletSeed> >::iterator it = fSeedPool.find(seedType);
    
          if (it != fSeedPool.end()){
            (it->second).push_back(DoubletSeed(seedType, hitj, hitk, initMom, initTheta, initPhi, charge, type));
          }
          else{
            cout<<"should never happen, fSeedPool should be init in constructor"<<endl;
            vector<DoubletSeed> thisVector;
            thisVector.push_back(DoubletSeed(seedType, hitj, hitk, initMom, initTheta, initPhi, charge, type));
            fSeedPool.insert(std::pair< SeedType, vector<DoubletSeed> >(seedType, thisVector));
          }
#ifdef MCDATA
          if (dynamic_cast<SoLIDMCGEMHit*>(hitj)->IsSignalHit() == 1 && dynamic_cast<SoLIDMCGEMHit*>(hitk)->IsSignalHit() == 1)
          fSeedEfficiency[0] = true;

          if (dynamic_cast<SoLIDMCGEMHit*>(hitj)->IsSignalHit() == 2 && dynamic_cast<SoLIDMCGEMHit*>(hitk)->IsSignalHit() == 2)
          fSeedEfficiency[1] = true;
#endif
        }
      }
      
    }
  }
  
}
//___________________________________________________________________________________________________________________
void SIDISKalTrackFinder::MergeSeed()
{ 
  //here we will merge the doublet seed into a triplet seed, for which the three type of doublet seed must match at the 
  //common plane. Once a triplet seed is form, its corresponding doublet seeds will be deactivated
  
  for (unsigned int i=0; i<fSeedPool[kMidBack].size(); i++){
    for (unsigned int j=0; j<fSeedPool[kFrontMid].size(); j++){
      //see if there is a common point in between the two doublet seeds
      if (fSeedPool[kMidBack].at(i).hita == fSeedPool[kFrontMid].at(j).hitb){
        //now check the third type of seed
        for (unsigned int k=0; k<fSeedPool[kFrontBack].size(); k++){
          if (fSeedPool[kFrontBack].at(k).hita == fSeedPool[kFrontMid].at(j).hita && 
            fSeedPool[kFrontBack].at(k).hitb == fSeedPool[kMidBack].at(i).hitb){
              
            fSeedPool[kMidBack].at(i).Deactive();
            fSeedPool[kFrontMid].at(j).Deactive();
            fSeedPool[kFrontBack].at(k).Deactive();
            
            SoLKalTrackSite & initSite =  SiteInitWithSeed(&(fSeedPool[kMidBack].at(i)));
            SoLKalTrackSystem *thisSystem = new ((*fCoarseTracks)[fNSeeds++]) SoLKalTrackSystem();
            //cout<<"n hits in track "<<thisSystem->GetLast()+1<<" "<<thisSystem->GetEntries()<<endl;
            thisSystem->SetCharge(fSeedPool[kMidBack].at(i).charge);
            if (!fHadronFlag){
                thisSystem->SetElectron(kTRUE);
                thisSystem->SetMass(kElectronMass);
            }
            else {
                thisSystem->SetElectron(kFALSE);
                thisSystem->SetMass(kPimMass);
            }
            thisSystem->SetAngleFlag(fSeedPool[kMidBack].at(i).flag);
            thisSystem->SetSeedType(kTriplet);
            thisSystem->SetOwner();
            thisSystem->Add(&initSite);
      
            //remember finding tracks always go backward
            SoLKalTrackSite& backSite = *new SoLKalTrackSite(fSeedPool[kMidBack].at(i).hitb, kMdim, kSdim, kMdim*fChi2PerNDFCut);
            if (!(thisSystem->AddAndFilter(backSite))) { thisSystem->SetTrackStatus(false); delete &backSite; }
      
            SoLKalTrackSite& midSite = *new SoLKalTrackSite(fSeedPool[kMidBack].at(i).hita, kMdim, kSdim, kMdim*fChi2PerNDFCut);
            if (!(thisSystem->AddAndFilter(midSite)))  { thisSystem->SetTrackStatus(false); delete &midSite; }
      
            SoLKalTrackSite& frontSite = *new SoLKalTrackSite(fSeedPool[kFrontBack].at(k).hita, kMdim, kSdim, kMdim*fChi2PerNDFCut);
            if (!(thisSystem->AddAndFilter(frontSite))) { thisSystem->SetTrackStatus(false); delete &frontSite; }
		         
          }
        }
      }
    }
  }
  //end of triplet seed matching and begin the remaining doublet seed init
  map< SeedType, vector<DoubletSeed> >::iterator it;
  for (it = fSeedPool.begin(); it != fSeedPool.end(); it++){
    vector<DoubletSeed> & thisVector = (it->second);
     
    for (unsigned int i=0; i<thisVector.size(); i++){
      if (!thisVector.at(i).isActive) continue;
      SoLKalTrackSite & initSite =  SiteInitWithSeed(&(thisVector.at(i)));
      SoLKalTrackSystem *thisSystem = new ((*fCoarseTracks)[fNSeeds++]) SoLKalTrackSystem();
      //cout<<"n hits in track "<<thisSystem->GetLast()+1<<" "<<thisSystem->GetEntries()<<endl;
      thisSystem->SetCharge(thisVector.at(i).charge);
      if (!fHadronFlag){
        thisSystem->SetElectron(kTRUE);
        thisSystem->SetMass(kElectronMass);
      }
      else {
        thisSystem->SetElectron(kFALSE);
        thisSystem->SetMass(kPimMass);
      }
      thisSystem->SetAngleFlag(thisVector.at(i).flag);
      thisSystem->SetSeedType(thisVector.at(i).type);
      thisSystem->SetOwner();
      thisSystem->Add(&initSite);
      //We assume that the doublet seed has already missed one hit (otherwise it is suppose to be part
      //of a triplet seed and thus be set as inactived already)
      thisSystem->AddMissingHits();
    
      SoLKalTrackSite& backSite = *new SoLKalTrackSite(thisVector.at(i).hitb, kMdim, kSdim, kMdim*fChi2PerNDFCut);
      if (!(thisSystem->AddAndFilter(backSite))) { thisSystem->SetTrackStatus(false); delete &backSite; }
      
      SoLKalTrackSite& midSite = *new SoLKalTrackSite(thisVector.at(i).hita, kMdim, kSdim,  kMdim*fChi2PerNDFCut);
      if (!(thisSystem->AddAndFilter(midSite)))  { thisSystem->SetTrackStatus(false); delete &midSite; }
    }
    
  }
  
}
//___________________________________________________________________________________________________________________
void SIDISKalTrackFinder::TrackFollow()
{
  //this function is responsible for propagating the seed track toward the next tracker, find suitable hits
  //the process stop until the track reach the first tracker upstream (track searching always go backward)
  for (Int_t i=0; i<fCoarseTracks->GetLast()+1; i++){
    SoLKalTrackSystem* thisSystem = (SoLKalTrackSystem*)(fCoarseTracks->At(i));
    thisSystem->CheckTrackStatus();
    if (!thisSystem->GetTrackStatus()) continue; //skip the bad tracks
    //--------------------for test----------------------//
      /*SoLKalTrackSystem *newSystem = (SoLKalTrackSystem*)thisSystem->Clone();
      cout<<"chi2: "<<newSystem->GetChi2()<<" "<<thisSystem->GetChi2()<<endl;
      cout<<"mass: "<<newSystem->GetMass()<<" "<<thisSystem->GetMass()<<endl;
      cout<<"state x: "<<(newSystem->GetCurSite()).GetCurState()(0, 0)<<" "<<(thisSystem->GetCurSite()).GetCurState()(0, 0)<<endl;
      cout<<"state x: "<<(newSystem->GetCurSite()).GetCurState()(1, 0)<<" "<<(thisSystem->GetCurSite()).GetCurState()(1, 0)<<endl;
      cout<<"state x: "<<(newSystem->GetCurSite()).GetCurState().GetZ0()<<" "<<(thisSystem->GetCurSite()).GetCurState().GetZ0()<<endl;*/
    //--------------------------------------------------//
    
    thisSystem->SetCurInstancePtr(thisSystem);
    Int_t currentTracker = ((thisSystem->GetCurSite()).GetHit())->GetTrackerID();
    
    Int_t lastTracker = 0;
    
    //seed from type kMidBack will skip the front seed plane. We assume for this type of seed, the hit on the
    //front seed plane is missing, (otherwise the seed should be absorbed into the triplet seed)
    
    if (thisSystem->GetSeedType() == kMidBack) currentTracker--;
    
    while (currentTracker > lastTracker){
      currentTracker--;
      
      thisSystem->CheckTrackStatus();    
      if (!thisSystem->GetTrackStatus()) break; //skip the bad tracks
      
      SoLKalTrackState currentState = (thisSystem->GetCurSite()).GetCurState();
      currentState.InitPredictSV();
      SoLKalTrackState *predictState = currentState.PredictSVatNextZ(fGEMTracker[currentTracker]->GetZ());
      
      bool flag = (thisSystem->GetNHits() >= 3);
      
      int size = GetHitsInWindow(currentTracker, (*predictState)(kIdxX0, 0), (predictState->GetCovMat())(kIdxX0, kIdxX0),
                                (*predictState)(kIdxY0, 0), (predictState->GetCovMat())(kIdxY0, kIdxY0), flag); 
                                
      
      if (size <= 0){
        //when there are too many hits in a small window (usually should not happen), or
        //when there is no hit found in the window, if the track has not missed a hit so far
        //we will keep the track, otherwise it is a bad track (miss too many hits)
        
        //the only exception will be the 0th tracker for a forward angle track, for which it 
        //is not necessary to have a hit, and in that case, we don't count missing hit
        if (thisSystem->GetAngleFlag() == kFAEC && currentTracker == 0) continue;
        
        thisSystem->AddMissingHits();
        
        /*if (thisSystem->GetNMissingHits() > 1 ){
          thisSystem->SetTrackStatus(kFALSE);
        }*/
      }
      else if (size == 1){
        SoLKalTrackSite &newSite = *new SoLKalTrackSite(fWindowHits.at(0), kMdim, kSdim, kMdim*fChi2PerNDFCut);
        newSite.Add(predictState);
        if (newSite.Filter()){
          thisSystem->Add(&newSite);
          thisSystem->IncreaseChi2(newSite.GetDeltaChi2());
          currentState.ClearAttemptSV();
        }
        else{
          thisSystem->AddMissingHits();
          delete &newSite;
        }
          
      }
      else{
        //find the cloest one for now, should use concurrent tracking in the future
        SoLKalTrackSite &newSite = *new SoLKalTrackSite(FindCloestHitInWindow((*predictState)(kIdxX0, 0), 
                                    (*predictState)(kIdxY0, 0)), kMdim, kSdim,  kMdim*fChi2PerNDFCut);
        newSite.Add(predictState);
        if (newSite.Filter()){
          thisSystem->Add(&newSite);
          thisSystem->IncreaseChi2(newSite.GetDeltaChi2());
          currentState.ClearAttemptSV();
        }
        else{
          thisSystem->AddMissingHits();
          delete &newSite;
        }
         
      }
     
    } 
    
    //now that we have all the hits selected, we can look at the chi2 per ndf and charge asymmetry to 
    //get rid of some potential bad tracks, before doing other things
    if (!thisSystem->GetTrackStatus()) continue;
    bool allMC[2] = {true, true};
    for (Int_t j=1; j!=thisSystem->GetLast()+1;j++){
      SoLIDGEMHit* thisHit = (SoLIDGEMHit*)((SoLKalTrackSite*)thisSystem->At(j))->GetHit();
      if (dynamic_cast<SoLIDMCGEMHit*>(thisHit)->IsSignalHit() != 1) allMC[0] = false;
      if (dynamic_cast<SoLIDMCGEMHit*>(thisHit)->IsSignalHit() != 2) allMC[1] = false;
    }
    if (allMC[0]) fMcTrackEfficiency[0] = true;
    if (allMC[1]) fMcTrackEfficiency[1] = true;
    
    
    if (thisSystem->GetChi2perNDF() > fChi2PerNDFCut) {
      thisSystem->SetTrackStatus(kFALSE);
      continue;
    }
    if (!CheckChargeAsy(thisSystem)){
      thisSystem->SetTrackStatus(kFALSE);
      continue;
    }
    //------------check MC track efficiency--------------//
    /*if (!thisSystem->GetTrackStatus()) continue;
    bool allMC[2] = {true, true};
    for (Int_t j=1; j!=thisSystem->GetLast()+1;j++){
      SoLIDGEMHit* thisHit = (SoLIDGEMHit*)((SoLKalTrackSite*)thisSystem->At(j))->GetHit();
      if (dynamic_cast<SoLIDMCGEMHit*>(thisHit)->IsSignalHit() != 1) allMC[0] = false;
      if (dynamic_cast<SoLIDMCGEMHit*>(thisHit)->IsSignalHit() != 2) allMC[1] = false;
    }
    if (allMC[0]) fMcTrackEfficiency[0] = true;
    if (allMC[1]) fMcTrackEfficiency[1] = true;*/
    //---------------------------------------------------//
  }
}
//___________________________________________________________________________________________________________________
void SIDISKalTrackFinder::FindandAddVertex()
{
   //cout<<fHadronFlag<<endl;
   fFieldStepper->UseFineStep();
   for (Int_t i=0; i<fCoarseTracks->GetLast()+1; i++){
      SoLKalTrackSystem* thisSystem = (SoLKalTrackSystem*)(fCoarseTracks->At(i));
      thisSystem->SetCurInstancePtr(thisSystem);
      //cout<<"track status "<<thisSystem->GetTrackStatus()<<endl; 
      thisSystem->CheckTrackStatus();
      if (thisSystem->GetTrackStatus() == kFALSE) continue; //skip bad tracks
      
      SoLKalTrackState currentState = (thisSystem->GetCurSite()).GetCurState();
      //currentState.InitPredictSV();
      
      //SoLKalTrackState *predictState = currentState.PredictSVatNextZ(fTargetCenter);
      //Double_t vertexz = FindVertexZ(predictState);
      
      //now we do the iterations here
      TVector3 initPosition(currentState(kIdxX0, 0), currentState(kIdxY0, 0), currentState.GetZ0());
      TVector3 initMomentum(0, 0, 0);
      TVector3 finalPosition(0, 0, 0);
      TVector3 finalMomentum(0, 0, 0);
      currentState.CalcMomVec(initMomentum);
      double vertexz = fTargetCenter;
      double step = 0.02;
      double charge = thisSystem->GetCharge();
      
      for (int iter = 0; iter < 3; iter++){
        fFieldStepper->PropagationClassicalRK4(initMomentum, initPosition, vertexz, 
                                               charge, step, finalMomentum, finalPosition);
        double tx = finalMomentum.X()/finalMomentum.Z();
        double ty = finalMomentum.Y()/finalMomentum.Z();
        vertexz   = vertexz + (1./(pow(tx,2) + pow(ty,2)))*
                           (tx*(fBPMX-finalPosition.X()) + ty*(fBPMY-finalPosition.Y()) );
      }
      
      if (thisSystem->GetAngleFlag() == kFAEC && fabs(vertexz - fTargetCenter) > (fTargetLength/2. + fCellEdgeCut[kFAEC]) ){
        thisSystem->SetTrackStatus(kFALSE);
        continue;
      }
      else if (thisSystem->GetAngleFlag() == kLAEC && fabs(vertexz - fTargetCenter) > (fTargetLength/2. + fCellEdgeCut[kLAEC]) ){
        thisSystem->SetTrackStatus(kFALSE);
        continue;
      }
      if (fHadronFlag && thisSystem->GetAngleFlag() == kFAEC && fElectronVertexZ > -1e5){
        //cout<<fElectronVertexZ<<" "<<vertexz<<" vertexz info"<<endl;
        if (fabs(vertexz - fElectronVertexZ) > 0.1) {
          thisSystem->SetTrackStatus(kFALSE);
          continue;
        }
      }
      
      
      //propagate the state vector to the interaction vertex that just found
      //not sure if this is the best way to add vertex
      currentState.InitPredictSV();
      SoLKalTrackState *predictState = currentState.PredictSVatNextZ(vertexz);
      
      //make a site at the interaction vertex to add to the fitting
      SoLKalTrackSite &vertexSite = *new SoLKalTrackSite(kMdim, kSdim,  kGiga);
      vertexSite.SetMeasurement(fBPMX, fBPMY);
      if (fDetConf == 3) vertexSite.SetHitResolution(1e-3, 1e-3);
      else vertexSite.SetHitResolution(3e-4, 3e-4);
      
      vertexSite.Add(predictState);
      if (vertexSite.Filter()){
        //calculate vertex variables and set info to the track system
        Double_t temp_tx =  vertexSite.GetCurState()(kIdxTX, 0);
        Double_t temp_ty =  vertexSite.GetCurState()(kIdxTY, 0);
        Double_t temp_qp =  vertexSite.GetCurState()(kIdxQP, 0);
        
        TVector3 vertex_vdir;
        vertex_vdir.SetZ( 1./(TMath::Sqrt(temp_tx*temp_tx + temp_ty*temp_ty + 1. )) );
        vertex_vdir.SetX(temp_tx * vertex_vdir.Z());
        vertex_vdir.SetY(temp_ty * vertex_vdir.Z());
        vertex_vdir = vertex_vdir.Unit();
        //cout<<fHadronFlag<<" "<<i<<" "<<thisSystem->GetCharge()/temp_qp<<endl; 
        thisSystem->SetMomentum(thisSystem->GetCharge()/temp_qp);
        thisSystem->SetTheta(acos(1./sqrt(1. + pow( (vertex_vdir.X()/vertex_vdir.Z()), 2) 
                             + pow((vertex_vdir.Y()/vertex_vdir.Z()) , 2))));
                             
        thisSystem->SetVertexZ(vertexz);
        thisSystem->SetPhi(atan2(vertex_vdir.Y(), vertex_vdir.X()));
      }
      else{
        thisSystem->SetTrackStatus(kFALSE); 
      }
      int detType = thisSystem->GetAngleFlag();
      if (detType == kFAEC && fHadronFlag) detType += 1;
      if (thisSystem->GetMomentum() > fMomMaxCut[detType] || 
          thisSystem->GetMomentum() < fMomMinCut[detType]) thisSystem->SetTrackStatus(kFALSE);
      
      currentState.ClearAttemptSV();
      delete &vertexSite;
   }
   fFieldStepper->UseDefaultStep();
}
//___________________________________________________________________________________________________________________
void SIDISKalTrackFinder::FinalSelection(TClonesArray *theTracks)
{
  fCoarseTracks->Sort();
  
  /*Int_t countTrack = 0;
  
  for (Int_t i=0; i<fCoarseTracks->GetLast()+1; i++){
    SoLKalTrackSystem *thisSystem = (SoLKalTrackSystem*)(fCoarseTracks->At(i));
    if (thisSystem->GetTrackStatus()) countTrack++; 
  }
  cout<<countTrack<<endl;*/
  
  for (Int_t i=0; i<fCoarseTracks->GetLast()+1; i++){
  
    SoLKalTrackSystem *thisSystem = (SoLKalTrackSystem*)(fCoarseTracks->At(i));
    thisSystem->SetCurInstancePtr(thisSystem);
    if (!thisSystem->GetTrackStatus()) continue;
    Int_t flag = 0;
    
    //start from 1 because the 0th is the dummy site that we used to initialize Kalman Filter
    //TODO remember not to add the last one since later it will be the BPM, not GEM hit
    for (Int_t j=1;j!=thisSystem->GetLast()+1;j++){
      SoLIDGEMHit *thishit = (SoLIDGEMHit*)(static_cast<SoLKalTrackSite*>(thisSystem->At(j))->GetHit());
      Int_t layer = thishit->GetTrackerID();
      
      map< Int_t, vector<SoLIDGEMHit*> >::iterator it = fGoodHits.find(layer);
      
      if (it != fGoodHits.end()){
        for (UInt_t n = 0; n<(it->second).size(); n++){
          if ((thishit->GetX() == ((it->second).at(n))->GetX()) && 
              (thishit->GetY() == ((it->second).at(n))->GetY())) { flag = 1; }
        }
      }
      
    }
    if (flag == 0){
      SoLIDTrack* newtrack = 0;
      if (fIsMC){
#ifdef MCDATA
        newtrack = new ((*theTracks)[fNGoodTrack++]) SoLIDMCTrack();
#endif
      }
      else{
        newtrack = new ((*theTracks)[fNGoodTrack++]) SoLIDTrack();
      }
    CopyTrack(newtrack, thisSystem);
    if (!fHadronFlag && fElectronVertexZ < -1e5) fElectronVertexZ = thisSystem->GetVertexZ();
    }
    
  }
  
}
//___________________________________________________________________________________________________________________
void SIDISKalTrackFinder::ECalFinalMatch()
{
  for (Int_t i=0; i<fCoarseTracks->GetLast()+1; i++){
    SoLKalTrackSystem *thisSystem = (SoLKalTrackSystem*)(fCoarseTracks->At(i));
    
    if ( !(thisSystem->GetTrackStatus()) ) continue; //skip bad tracks
    Double_t ecalZ = fECal->GetECZ((ECType)thisSystem->GetAngleFlag());
    
    //using Kalman Filter smoother to smooth the track back to the first measurement site
    //so that we don't need to propagate and fit back again
    thisSystem->SetCurInstancePtr(thisSystem);
    thisSystem->SmoothBackTo(1);
    
    SoLKalTrackState currentState = (thisSystem->GetCurSite()).GetCurState();
    currentState.InitPredictSV();
    
    SoLKalTrackState *predictState = currentState.PredictSVatNextZ(ecalZ);
    
    thisSystem->SetTrackStatus(kFALSE);
    for (Int_t ec_count=0; ec_count<fCaloHits->GetLast()+1; ec_count++){
        SoLIDCaloHit* thisHit = ((SoLIDCaloHit*)fCaloHits->At(ec_count));
    
	    if (thisHit->fECID != thisSystem->GetAngleFlag()) continue;
	    if ( fabs(thisHit->fXPos - (*predictState)(kIdxX0, 0)) <fECPosCut &&
	         fabs(thisHit->fYPos - (*predictState)(kIdxY0, 0)) <fECPosCut  ){
	      
	      //for large angle, require also that the momentum of the track needs to match the 
	      //cluster energy. This is difficult to do for forward angle since we detect both hadron
	      //and electron there and there is a long distance betwee the FAEC and the last GEM, during
	      //which there could be significant energy loss but we don't have other tracking detectors
	      //and field integral to measure it
	      Int_t type = thisSystem->GetAngleFlag();
	      Double_t momentum = thisSystem->GetCharge() / (*predictState)(kIdxQP, 0);
	      Double_t e_reso = fECal->GetEReso()/sqrt(thisHit->fEdp);
	      if (fECEnergyMatch[type]){  
	        if (type == kLAEC){
	            if ( ((thisHit->fEdp) - momentum)/thisHit->fEdp > 3.*e_reso ) continue;
	            if ( fabs( (momentum - thisHit->fEdp) / momentum) > 0.5) continue;
	        }else{
	            if ( ((thisHit->fEdp) - momentum)/thisHit->fEdp > 3.*e_reso ) continue;
	        }
	      }
	      thisSystem->SetTrackStatus(kTRUE);
	      thisSystem->fDeltaECX = thisHit->fXPos;// - (*predictState)(kIdxX0, 0);
	      thisSystem->fDeltaECY = thisHit->fYPos;// - (*predictState)(kIdxY0, 0);
	      thisSystem->fDeltaECE = thisHit->fEdp;    
	    }
	  }
	  
	  thisSystem->SetSitePtrToLastSite();
  }
}
//___________________________________________________________________________________________________________________

//assistant functions below
inline SoLKalTrackSite & SIDISKalTrackFinder::SiteInitWithSeed(DoubletSeed* thisSeed)
{
  TVector3 initDir(cos(thisSeed->initPhi), sin(thisSeed->initPhi), 
                   1./tan(thisSeed->initTheta));
  initDir = initDir.Unit();
            
  //-----------prepare seeds for Kalman Filter track finding------------//
  SoLKalMatrix svd(kSdim,1);
  svd.Zero();
  svd(kIdxX0,0) = (thisSeed->hitb)->GetX();
  svd(kIdxY0,0) = (thisSeed->hitb)->GetY();
  svd(kIdxTX,0) = initDir.X()/initDir.Z();
  svd(kIdxTY,0) = initDir.Y()/initDir.Z();
  svd(kIdxQP,0) = thisSeed->charge/thisSeed->initMom;
  
  SoLKalMatrix C(kSdim,kSdim);
  C.Zero();
  /*for (int index=0; index<kSdim; index++) {
    C(index,index) = 1.;   // dummy error matrix
  }*/
  Double_t phi = (thisSeed->hitb)->GetPhi();
  Double_t dr = 6.e-4;
  Double_t drphi = 6.4e-5;
  
  Double_t dx = sqrt( pow( cos(phi)*dr, 2) + pow( sin(phi)*drphi, 2) );
  Double_t dy = sqrt( pow( sin(phi)*dr, 2) + pow( cos(phi)*drphi, 2) );
  
  C(kIdxX0, kIdxX0) = 10*pow(dx, 2);
  C(kIdxY0, kIdxY0) = 10*pow(dy, 2);
  C(kIdxTX, kIdxTX) = 0.001; 
  C(kIdxTY, kIdxTY) = 0.001;
  C(kIdxQP, kIdxQP) = 0.01;
      
  SoLKalTrackSite& initSite = *new SoLKalTrackSite(thisSeed->hitb, kMdim, kSdim, kMdim*fChi2PerNDFCut);
      
  initSite.Add(new SoLKalTrackState(svd, C, initSite, SoLKalTrackSite::kPredicted));
  initSite.Add(new SoLKalTrackState(svd, C, initSite, SoLKalTrackSite::kFiltered));
  initSite.SetHitResolution(kGiga, kGiga); //give it a very large resolution (100m) since it is a virtual site
  
  return initSite;
}
//___________________________________________________________________________________________________________________
inline bool SIDISKalTrackFinder::TriggerCheck(SoLIDGEMHit *theHit, ECType type, Double_t& matchEnergy)
{
  if (type == kLAEC){
    for (Int_t ec_count=0; ec_count<fCaloHits->GetLast()+1; ec_count++){
        SoLIDCaloHit* thisHit = (SoLIDCaloHit*)(fCaloHits->At(ec_count));
    
	    if (thisHit->fECID != kLAEC) continue; //not LAEC hit
	    Double_t ecHitPhi = TMath::ATan2(thisHit->fYPos, thisHit->fXPos);
	    Double_t ecHitR = TMath::Sqrt( TMath::Power(thisHit->fXPos, 2) + 
	    			  	   TMath::Power(thisHit->fYPos, 2) );
	    Double_t tmpDeltaPhi = CalDeltaPhi(ecHitPhi, theHit->GetPhi());
	    Double_t tmpDeltaR   = CalDeltaR(ecHitR, theHit->GetR());
	    if (theHit->GetTrackerID()==2 && (tmpDeltaPhi < fDeltaPhiMaxBC[1][0] && tmpDeltaPhi > fDeltaPhiMinBC[1][0]) && 
	       (tmpDeltaR < fDeltaRMaxBC[1][0] && tmpDeltaR > fDeltaRMinBC[1][0])){
            matchEnergy = thisHit->fEdp;
            return true;
        }else if (theHit->GetTrackerID()==3 && (tmpDeltaPhi<fDeltaPhiMaxBC[0][0] && tmpDeltaPhi> fDeltaPhiMinBC[0][0]) && 
                 (tmpDeltaR<fDeltaRMaxBC[0][0] && tmpDeltaR > fDeltaRMinBC[0][0]) ){
            matchEnergy = thisHit->fEdp;
            return true;
        }
	  }
  }
  else if (type == kFAEC){
    for (Int_t ec_count=0; ec_count<fCaloHits->GetLast()+1; ec_count++)
	    {
	      SoLIDCaloHit* thisHit = (SoLIDCaloHit*)(fCaloHits->At(ec_count));
	      if (thisHit->fECID != kFAEC) continue; //not FAEC hit
	      if (fHadronFlag) {
	        if (thisHit->fSubDetID == 0) continue;
	      }
	      else{
	        if (thisHit->fSubDetID != 0) continue;
	      }
	      Double_t ecHitPhi = TMath::ATan2(thisHit->fYPos, thisHit->fXPos);
	      Double_t ecHitR = TMath::Sqrt( TMath::Power(thisHit->fXPos, 2) + 
				  	   TMath::Power(thisHit->fYPos, 2) );
	      Double_t tmpDeltaPhi = CalDeltaPhi(ecHitPhi, theHit->GetPhi());
	      Double_t tmpDeltaR   = CalDeltaR(ecHitR, theHit->GetR());
	      Int_t    dettype     = thisHit->fSubDetID+1;
	      if (theHit->GetTrackerID()==4 && (tmpDeltaPhi < fDeltaPhiMaxBC[1][dettype] && tmpDeltaPhi > fDeltaPhiMinBC[1][dettype]) && 
	         (tmpDeltaR < fDeltaRMaxBC[1][dettype] && tmpDeltaR > fDeltaRMinBC[1][dettype])){
	        matchEnergy = thisHit->fEdp;
	        return true;
	      }else if (theHit->GetTrackerID()==5 && (tmpDeltaPhi<fDeltaPhiMaxBC[0][dettype] && tmpDeltaPhi> fDeltaPhiMinBC[0][dettype]) && 
	               (tmpDeltaR<fDeltaRMaxBC[0][dettype] && tmpDeltaR > fDeltaRMinBC[0][dettype]) ){
	        matchEnergy = thisHit->fEdp;
	        return true;
	      }
	    }
  }
  return false;
}
//___________________________________________________________________________________________________________________
inline SoLIDGEMHit* SIDISKalTrackFinder::FindCloestHitInWindow(double &x, double &y){
  double minD = kGiga;
  SoLIDGEMHit *minHit = nullptr;
  for (unsigned int i=0; i<fWindowHits.size(); i++){
    double r = sqrt(pow(x - fWindowHits.at(i)->GetX(), 2) + pow(y - fWindowHits.at(i)->GetY(), 2));
    if (r < minD) {
      minHit = fWindowHits.at(i);
      minD = r;
    }
  }
  return minHit;
}
//___________________________________________________________________________________________________________________
inline int SIDISKalTrackFinder::GetHitsInWindow(int plane, double x, double wx, double y, double wy, bool flag)
{
  assert(plane >= 0);
  fWindowHits.clear();
  vector<Int_t> ChamberList;
  
  double thisR = sqrt(x*x + y*y);
  
  
  GetHitChamberList(ChamberList, GetChamIDFromPos(x, y, plane), 1);
  for (int i=0; i<(int)ChamberList.size(); i++){
    TSeqCollection* HitArray = fGEMTracker[plane]->GetChamber(ChamberList.at(i))->GetHits();
    for (int nhit = 0; nhit < HitArray->GetLast()+1; nhit++){ 
      SoLIDGEMHit *hit = (SoLIDGEMHit*)HitArray->At(nhit);
      
      if (hit->IsUsed()) continue;
      if (hit->GetR() < thisR - 0.03) continue;
      if (hit->GetR() > thisR + 0.03) break;     
      
      bool condition;
      if (!flag) condition = sqrt( pow(hit->GetX() - x, 2) + pow(hit->GetY() - y, 2) ) < 0.015 ;
      else condition = ( fabs(hit->GetX() - x) < 10.*sqrt(wx) && fabs(hit->GetY() - y) < 10.*sqrt(wy) ); 
      
      
      
      if (condition){
        fWindowHits.push_back(hit);
        if (fWindowHits.size() > MAXWINDOWHIT) return -1; //too many hits to be considered       
      }        
    }      
  }
  
  return fWindowHits.size();
}
//____________________________________________________________________________________________________________________
inline Double_t SIDISKalTrackFinder::FindVertexZ(SoLKalTrackState* thisState)
{
  //cout<<thisState->GetZ0()<<" "<<fTargetCenter<<" "<<(*thisState)(kIdxX0, 0)<<" "<<(*thisState)(kIdxY0,0)<<endl;
  //assert(fabs(thisState->GetZ0() - fTargetCenter) < 0.01);
  Double_t vertexz = fTargetCenter+ (1./(pow((*thisState)(kIdxTX, 0),2) + pow((*thisState)(kIdxTY, 0),2)))*
  ((*thisState)(kIdxTX,0)*(fBPMX-(*thisState)(kIdxX0,0)) + (*thisState)(kIdxTY,0)*(fBPMY-(*thisState)(kIdxY0, 0)));
  
  return vertexz;
}
//____________________________________________________________________________________________________________________
inline void SIDISKalTrackFinder::CopyTrack(SoLIDTrack* soltrack, SoLKalTrackSystem* kaltrack)
{
  soltrack->SetStatus(kaltrack->GetTrackStatus());
  soltrack->SetCoarseFitStatus(kTRUE);
  soltrack->SetAngleFlag(kaltrack->GetAngleFlag());
  soltrack->SetCharge(kaltrack->GetCharge());
  soltrack->SetNDF(kaltrack->GetNDF());
  soltrack->SetCoarseChi2(kaltrack->GetChi2perNDF());
  soltrack->SetMomentum(kaltrack->GetMomentum());
  soltrack->SetVertexZ(kaltrack->GetVertexZ());
  soltrack->SetPhi(kaltrack->GetPhi());
  soltrack->SetTheta(kaltrack->GetTheta());
  soltrack->SetECX(kaltrack->fDeltaECX);
  soltrack->SetECY(kaltrack->fDeltaECY);
  soltrack->SetECE(kaltrack->fDeltaECE);
  soltrack->SetElectron(kaltrack->IsElectron());
  //cout<<fNGoodTrack<<" "<<kaltrack->GetMomentum()<<" "<<kaltrack->GetTheta()<<endl;
  
  //start from 1 because the 0th is the dummy site that we used to initialize Kalman Filter
  //TODO remember not to add the last one since later it will be the BPM, not GEM hit
  
  for (Int_t j=1; j!=kaltrack->GetLast()+1;j++){
    SoLIDGEMHit* thishit = 0;
    Int_t layer = 0;
    
    thishit = (SoLIDGEMHit*)(static_cast<SoLKalTrackSite*>(kaltrack->At(j))->GetPredInfoHit());
    thishit->SetUsed();
    layer = thishit->GetTrackerID();
    
    map< Int_t, vector<SoLIDGEMHit*> >::iterator it = fGoodHits.find(layer);
    
    if (it != fGoodHits.end()){
      (it->second).push_back(thishit);
    }
    else{
      vector<SoLIDGEMHit*> thisVector;
      thisVector.push_back(thishit);
      fGoodHits.insert(std::pair<Int_t, vector<SoLIDGEMHit*> >(layer, thisVector));
    }
      
    assert(thishit != 0);
    soltrack->AddHit(thishit);
    //cout<<layer<<" hit info "<<thishit->GetX()<<" "<<thishit->GetY()<<endl;
  }
  
  //for back track info
  SoLKalTrackState & thisState = dynamic_cast<SoLKalTrackSite*>(kaltrack->At(kaltrack->GetLast()))->GetState(SoLKalTrackSite::kFiltered);
  Double_t &x  = thisState(kIdxX0, 0);
  Double_t &y  = thisState(kIdxY0, 0);
  Double_t &tx = thisState(kIdxTX, 0);
  Double_t &ty = thisState(kIdxTY, 0);
  soltrack->SetBackX(x);
  soltrack->SetBackY(y);
  
  TVector3 vdir;
  vdir.SetZ( 1./(TMath::Sqrt(tx*tx + ty*ty + 1. )) );
  vdir.SetX(tx * vdir.Z());
  vdir.SetY(ty * vdir.Z());
  vdir = vdir.Unit();
  
  soltrack->SetBackTheta(acos(1./sqrt(1. + pow( (vdir.X()/vdir.Z()), 2)
                             + pow((vdir.Y()/vdir.Z()) , 2))));

  soltrack->SetBackPhi(atan2(vdir.Y(), vdir.X()));
}
//______________________________________________________________________________________________________________________
inline Bool_t SIDISKalTrackFinder::CheckChargeAsy(SoLKalTrackSystem* theSystem)
{
  Int_t countCharge = 0;
  
  for (Int_t i=1; i<theSystem->GetLast()+1; i++){
    SoLIDGEMHit* theHit = (SoLIDMCGEMHit*)(static_cast<SoLKalTrackSite*>(theSystem->At(i))->GetHit());
    if (fabs( (theHit->GetQU() - theHit->GetQV())/(theHit->GetQU() + theHit->GetQV()) ) < 0.6) countCharge++;
  }
  if (countCharge >= 3){
    return kTRUE;
  }
  else{
    return kFALSE;
  }
}
//_______________________________________________________________________________________________________________________
inline void SIDISKalTrackFinder::GetHitChamberList(vector<Int_t> &theList, Int_t thisChamber, Int_t size)
{
  //choose 7 chambers around the current hit chamber, in order to reduce the search region 
  theList.resize(2*size+1);
  thisChamber -= size;
  for (unsigned int i=0; i<theList.size(); i++){
    theList[i] = thisChamber;
    thisChamber++;
    if (theList[i] < 0) theList[i] += 30;
    if (theList[i] >29) theList[i] -= 30;
  }
}
//_______________________________________________________________________________________________________________________
inline Int_t SIDISKalTrackFinder::GetChamIDFromPos(Double_t &x, Double_t &y, Int_t TrackerID)
{
  double phi = atan2(y, x);
  for (int i=0; i<fGEMTracker[TrackerID]->GetNChamber(); i++){
    SoLIDGEMChamber* thisChamber = fGEMTracker[TrackerID]->GetChamber(i);
    double dphi = phi - thisChamber->GetPhiInLab();
    dphi = TVector2::Phi_mpi_pi(dphi);
    if ( dphi < thisChamber->GetPhiCover()/2. && dphi > -1*thisChamber->GetPhiCover()/2.) 
    return i;
  }
  return -1; 
}
//_______________________________________________________________________________________________________________________
inline double SIDISKalTrackFinder::PredictR(Int_t &plane, SoLIDGEMHit* hit1, SoLIDGEMHit* hit2)
{
  double z = fGEMTracker[plane]->GetZ();
  return (hit1->GetR() - hit2->GetR())/(hit1->GetZ() - hit2->GetZ())*(z - hit1->GetZ()) 
         + hit1->GetR();
}
//________________________________________________________________________________________________________________________
inline Bool_t SIDISKalTrackFinder::CalInitParForPair(SoLIDGEMHit* hita, SoLIDGEMHit* hitb, Double_t &charge,
                                                 Double_t& mom, Double_t& theta, Double_t& phi, ECType& type)
{
    double deltaR = hitb->GetR() - hita->GetR();
    double deltaPhi = TVector2::Phi_mpi_pi(hitb->GetPhi() - hita->GetPhi());
    deltaPhi = fabs(deltaPhi);

    double deltaZ = hitb->GetZ() - hita->GetZ();
    double deltaY = hitb->GetY() - hita->GetY();
    double deltaX = hitb->GetX() - hita->GetX();

    theta = atan(deltaR/deltaZ);
    phi = atan2(deltaY, deltaX);
  
    if (fDetConf == 3){ //SIDIS NH3
        if (type == kFAEC && hita->GetTrackerID() == 4 && hitb->GetTrackerID() == 5){
            if (deltaPhi < -3.54398e-04 + 1e-6) return kFALSE;
            mom = 1.82273e-01/(deltaPhi - -3.54398e-04);
            theta += -3.15393e-05 + -0.00194545*deltaPhi + 1.92086*deltaPhi*deltaPhi;
            phi += 0.000503764 + 0.843494*deltaPhi + 0.613864*deltaPhi*deltaPhi;
        }
        else if (type == kFAEC && hita->GetTrackerID() == 3 && hitb->GetTrackerID() == 4){
            if (deltaPhi < -1.44365e-04 + 1e-6) return kFALSE;
            mom = 1.52902e-01/(deltaPhi - -1.44365e-04);
            theta += -6.27864e-05 + 0.00150745*deltaPhi + 1.87546*deltaPhi*deltaPhi;
            phi +=  0.00283427+ 0.925668*deltaPhi + 1.29002*deltaPhi*deltaPhi;
        }
        else if (type == kFAEC && hita->GetTrackerID() == 3 && hitb->GetTrackerID() == 5){
            if (deltaPhi < -2.67014e-04 + 1e-6) return kFALSE;
            mom = 3.40491e-01/(deltaPhi - -2.67014e-04);
            theta +=  -0.000104712+ 0.00152549*deltaPhi + 0.451974*deltaPhi*deltaPhi;
            phi += 0.00410935 + 0.832924*deltaPhi + 0.607294*deltaPhi*deltaPhi;
        }
        else if (type == kLAEC && hita->GetTrackerID() == 2 && hitb->GetTrackerID() == 3){
            if (deltaPhi < -1.18851e-03 + 1e-6) return kFALSE;
            mom = 1.12414e-01/(deltaPhi - -1.18851e-03);
            theta +=  0.000239174+ -0.0296545*deltaPhi + 5.32953*deltaPhi*deltaPhi;
            phi += 0.00472525+ 0.88319*deltaPhi + 5.14816*deltaPhi*deltaPhi;
        }
        else if (type == kLAEC && hita->GetTrackerID() == 1 && hitb->GetTrackerID() == 2){
            if (deltaPhi < -9.32190e-04 + 1e-6) return kFALSE;
            mom = 6.42127e-02/(deltaPhi - -9.32190e-04);
            theta += 4.88206e-05 + -0.0184821*deltaPhi + 9.93942*deltaPhi*deltaPhi;
            phi += 0.00560415+ 0.608656*deltaPhi + 20.6661*deltaPhi*deltaPhi;
        }
        else if (type == kLAEC && hita->GetTrackerID() == 1 && hitb->GetTrackerID() == 3){
            if (deltaPhi < -2.10386e-03 + 1e-6) return kFALSE;
            mom = 1.76726e-01/(deltaPhi - -2.10386e-03);
            theta += 2.66966e-05 + -0.00712552*deltaPhi + 1.88514*deltaPhi*deltaPhi;
            phi +=0.0103316 + 0.798952*deltaPhi + 4.67825*deltaPhi*deltaPhi;
        }
    }
    else{ //SIDIS He3 and JPsi
        if (type == kFAEC && hita->GetTrackerID() == 4 && hitb->GetTrackerID() == 5){
            if (deltaPhi < -0.000260251 + 1e-6) return kFALSE;
            mom = 0.184566/(deltaPhi - -0.000260251);
            theta += 0.00117211 + -0.0519514*deltaPhi + 2.25303*deltaPhi*deltaPhi;
            phi +=(-1.*charge)*( 0.000617591 + 0.827572*deltaPhi + 0.53233*deltaPhi*deltaPhi );
        }
        else if (type == kFAEC && hita->GetTrackerID() == 3 && hitb->GetTrackerID() == 4){
            if (deltaPhi < -0.000456619 + 1e-6) return kFALSE;
            mom = 0.156196/(deltaPhi - -0.000456619);
            theta += 0.000614522 + -0.0317596*deltaPhi + 2.14209*deltaPhi*deltaPhi;
            phi += (-1*charge)*(-0.000117883 + 1.0794*deltaPhi + -0.178718*deltaPhi*deltaPhi);
        }
        else if (type == kFAEC && hita->GetTrackerID() == 3 && hitb->GetTrackerID() == 5){
            if (deltaPhi < -0.00068094 + 1e-6) return kFALSE;
            mom = 0.340493/(deltaPhi - -0.00068094);
            theta += 0.000992593 + -0.0239743*deltaPhi + 0.549664*deltaPhi*deltaPhi;
            phi += (-1*charge)*(0.000672179 + 0.908331*deltaPhi + 0.166936*deltaPhi*deltaPhi);
        }
        else if (type == kLAEC && hita->GetTrackerID() == 2 && hitb->GetTrackerID() == 3){
            if (deltaPhi < -0.000396581 + 1e-6) return kFALSE;
            mom = 0.109309/(deltaPhi - -0.000396581);
            theta += 3.61953e-05 + -1.06941e-05*deltaPhi + 4.23056*deltaPhi*deltaPhi;
            phi += (-1*charge)*(-0.000283361 + 1.2515*deltaPhi + -0.637376*deltaPhi*deltaPhi);
        }
        else if (type == kLAEC && hita->GetTrackerID() == 1 && hitb->GetTrackerID() == 2){
            if (deltaPhi < -0.000371282 + 1e-6) return kFALSE;
            mom = 0.0610608/(deltaPhi - -0.000371282);
            theta += 5.30086e-05 + -0.00477822*deltaPhi + 8.29537*deltaPhi*deltaPhi;
            phi += (-1*charge)*(-0.00018356 + 1.39156*deltaPhi + -1.29865*deltaPhi*deltaPhi);
        }
        else if (type == kLAEC && hita->GetTrackerID() == 1 && hitb->GetTrackerID() == 3){
            if (deltaPhi < -0.000743132 + 1e-6) return kFALSE;
            mom = 0.17019/(deltaPhi - -0.000743132);
            theta += -5.87987e-06 + 0.00195822*deltaPhi + 1.58111*deltaPhi*deltaPhi;
            phi += (-1*charge)*(-0.000508897 + 1.30589*deltaPhi + -0.468172*deltaPhi*deltaPhi);
        }
    }
  
    phi = TVector2::Phi_mpi_pi(phi);
    return kTRUE;
  
}
//________________________________________________________________________________________________________________________
inline Int_t SIDISKalTrackFinder::BinarySearchForR(TSeqCollection* array, Double_t &lowr)
{
  //search for the first index that is above lowr, if not found, return the length of the array
  //the array needs to be sorted in increasing r order before use
  
  Int_t low = 0;
  Int_t high = array->GetEntries();
  
  while (low != high){
    Int_t mid = (low + high)/2;
    SoLIDGEMHit* hit = (SoLIDGEMHit*)array->At(mid);
    if ( hit->GetR() <= lowr ){
      
      low = mid + 1;
    }
    else{
      high = mid;
    }
  }
  
  assert(low == high);
  return low;
}
//__________________________________________________________________________
inline double SIDISKalTrackFinder::CalDeltaPhi(const double & phi1, const double & phi2)
{
  double deltaPhi = phi1 - phi2;
  return TVector2::Phi_mpi_pi(deltaPhi);
}
//___________________________________________________________________________
inline double SIDISKalTrackFinder::CalDeltaR(const double & r1, const double & r2)
{
  return r1 - r2;
}


