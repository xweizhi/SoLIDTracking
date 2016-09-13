//ROOT
#include "TVector2.h"
#include "TMath.h"
//SoLIDTracking
#include "ProgressiveTracking.h"
#include "SoLIDUtility.h"

ProgressiveTracking::ProgressiveTracking(Int_t ntracker, Bool_t isMC) : 
fNTracker(ntracker), fDoMC(isMC),
fNElectron(1), fNHadron(0), fNTrack(0), fIsIterBackward(kTRUE), fHasCaloHit(kTRUE)
{

  if (fDoMC){
#ifdef MCDATA
    fCoarseTracks = new TClonesArray("SoLIDMCTrack", 200);
    for (Int_t i=0; i<fNTracker; i++){
      fGoodHits[i] = new TClonesArray("SoLIDMCGEMHit", 20);
      fNGoodHits[i] = 0;
#endif
    }
  }
  else{
    fCoarseTracks = new TClonesArray("SoLIDTrack", 200);
    for (Int_t i=0; i<fNTracker; i++){
      fGoodHits[i] = new TClonesArray("SoLIDGEMHit", 20);
      fNGoodHits[i] = 0;
    }
  }
}
//____________________________________________________________________________
ProgressiveTracking::~ProgressiveTracking()
{
  Clear();
  delete fCoarseTracks;
  for (Int_t i=0; i<fNTracker; i++){
    delete fGoodHits[i];
  }
}
//____________________________________________________________________________
void ProgressiveTracking::ProcessHits(map<Int_t, vector<TSeqCollection*> > *theHitMap, 
				      TClonesArray* theTracks)
{
  //temperary for testing SoLIDTrack class
  /*SoLIDGEMHit *hitSave[6] = {0};
  
  for (UInt_t i=0; i< (*theHitMap).size()-2; i++){
    vector<TSeqCollection*> hitsOnTracker = (*theHitMap).find(i)->second;
    for (UInt_t j=0; j<hitsOnTracker.size(); j++){
      TSeqCollection* hitarray = hitsOnTracker.at(j);
      TIterator* hitit = hitarray->MakeIterator();
      SoLIDGEMHit *thisHit = 0;
      while ( ( thisHit = static_cast<SoLIDGEMHit*>(hitit->Next()) ) ){
        if (dynamic_cast<SoLIDMCGEMHit*>(thisHit)->IsSignalHit()){
          hitSave[i] = thisHit;
          break;
        }
      }
    }
  }
  Int_t hitCount = 0;
  for (UInt_t i=0; i<(*theHitMap).size(); i++){
    if (hitSave[i] != 0){
      hitCount++;
    }
  }
  
  if (hitCount>0){
    new ( (*theTracks)[0] ) SoLIDMCTrack();
    
    for (UInt_t i=0; i<(*theHitMap).size(); i++){
      if (hitSave[i] != 0){
        dynamic_cast<SoLIDMCTrack*>(theTracks->At(0))->AddHit(hitSave[i]);
      }
    }
  }*/
  //--------------------------------------
  FindTrack(0, 0, theHitMap);
  FindTrack(0, 1, theHitMap);
  FindTrack(0, 2, theHitMap);
  FindTrack(0, 3, theHitMap);
  FindTrack(0, 4, theHitMap);
  //FindTrack(1, 0, theHitMap);
  //FindTrack(1, 1, theHitMap);
  //FindTrack(1, 2, theHitMap);
  //FindTrack(1, 3, theHitMap);
  //FindTrack(1, 4, theHitMap);
  //FindTrack(1, 5, theHitMap);
  //cout<<fNTrack<<endl;
  CheckTracks();
  CombineTrackRoad(theTracks);
  //cout<<"after combine tracks: "<<theTracks->GetLast()+1<<endl;
}
//_______________________________________________________________________________
Int_t ProgressiveTracking::ReadDataBase()
{
  //TODO: move all the hard-coding to here
  return 0; 
}
//_______________________________________________________________________________
void ProgressiveTracking::Clear( Option_t* opt )
{
  fCoarseTracks->Clear(opt);
  fCaloHits.clear();
  fNTrack = 0;
  for (Int_t i=0; i<fNTracker; i++){
    fGoodHits[i]->Clear(opt);
    fNGoodHits[i] = 0;
  }
}
//________________________________________________________________________________
void ProgressiveTracking::FindTrack(Int_t angleflag, Int_t type, 
                                    map<Int_t, vector<TSeqCollection*> > *theHitMap)
{
  Double_t philimit[4][2];
  Double_t rlimit[5][2];
  Double_t deltar[4][2];
  Double_t dr, dphi;

  Int_t layer_array[5],nlayer=0;
  Int_t Is_Seed = 1;
  if (angleflag == kLAEC && type == 0){
    // large angle 1->2->3->4
    layer_array[0] = 3;layer_array[1] = 2;layer_array[2] = 1;layer_array[3] = 0;nlayer=4;
    rlimit[3][0] = 0.434; rlimit[3][1] = 0.860; rlimit[2][0] = 0.51; rlimit[2][1] = 0.98;
    rlimit[1][0] = 0.60; rlimit[1][1] = 1.12; rlimit[0][0] = 0.74; rlimit[0][1] = 1.35;
    rlimit[4][0] = 0.; rlimit[4][1] = 0.;
    philimit[2][0]=0.003; philimit[2][1]=0.05; philimit[1][0]=0.006;philimit[1][1]=0.07;
    philimit[0][0]=0.012; philimit[0][1]=0.127; philimit[3][0]=0;philimit[3][1]=0;
    deltar[2][0]=0.058; deltar[2][1]=0.122; deltar[1][0]=0.076; deltar[1][1]=0.151;
    deltar[0][0]=0.121; deltar[0][1]=0.242; deltar[3][0]=0.;deltar[3][1]=100.;
  }else if (angleflag == kLAEC && type == 1){
    // large angle 1->3->4
    layer_array[0] = 3;layer_array[1] = 2;layer_array[2] = 0;nlayer=3;
    rlimit[2][0] = 0.434; rlimit[2][1] = 0.860; rlimit[1][0] = 0.6; rlimit[1][1] = 1.12;
    rlimit[0][0] = 0.74; rlimit[0][1] = 1.35; rlimit[3][0] = 0.; rlimit[3][1] = 0.;
    rlimit[4][0] = 0.; rlimit[4][1] = 0.;
    philimit[1][0]=0.01; philimit[1][1]=0.12; philimit[0][0]=0.012;philimit[0][1]=0.127;
    philimit[2][0]=0; philimit[2][1]=0; philimit[3][0]=0; philimit[3][1]=0;
    deltar[1][0]=0.139; deltar[1][1]=0.269; deltar[0][0]=0.121; deltar[0][1]=0.242;
    deltar[2][0]=0.; deltar[2][1]=100.; deltar[3][0]=0.;deltar[3][1]=100.;
  }else if (angleflag == kLAEC && type == 2){
    // large angle 1->2->4
    rlimit[2][0] = 0.434; rlimit[2][1] = 0.860; rlimit[1][0] = 0.51; rlimit[1][1] = 0.98;
    rlimit[0][0] = 0.74; rlimit[0][1] = 1.35; rlimit[3][0] = 0.; rlimit[3][1] = 0.;
    rlimit[4][0] = 0.; rlimit[4][1] = 0.;
    layer_array[0] = 3;layer_array[1] = 1;layer_array[2] = 0;nlayer=3;
    philimit[1][0]=0.003; philimit[1][1]=0.05; philimit[0][0]=0.018;philimit[0][1]=0.2;
    philimit[2][0]=0; philimit[2][1]=0; philimit[3][0]=0; philimit[3][1]=0;
    deltar[1][0]=0.058; deltar[1][1]=0.122; deltar[0][0]=0.2; deltar[0][1]=0.39;
    deltar[2][0]=0.; deltar[2][1]=100.; deltar[3][0]=0.;deltar[3][1]=100.;
  }else if (angleflag == kLAEC && type == 3){
    // large angle 1->2->3
    rlimit[2][0] = 0.434; rlimit[2][1] = 0.860; rlimit[1][0] = 0.51; rlimit[1][1] = 0.98;
    rlimit[0][0] = 0.6; rlimit[0][1] = 1.12; rlimit[3][0] = 0.; rlimit[3][1] = 0.;
    rlimit[4][0] = 0.; rlimit[4][1] = 0.;
    layer_array[0] = 2;layer_array[1] = 1;layer_array[2] = 0;nlayer=3;
    philimit[1][0]=0.003; philimit[1][1]=0.05; philimit[0][0]=0.006;philimit[0][1]=0.07;
    philimit[2][0]=0; philimit[2][1]=0; philimit[3][0]=0; philimit[3][1]=0;
    deltar[1][0]=0.058; deltar[1][1]=0.122; deltar[0][0]=0.076; deltar[0][1]=0.151;
    deltar[2][0]=0.; deltar[2][1]=100.; deltar[3][0]=0.;deltar[3][1]=100.;
  }else if (angleflag == kLAEC && type == 4){
    // large angle 2->3->4
    layer_array[0] = 3;layer_array[1] = 2;layer_array[2] = 1;nlayer=3;
    //rlimit[0]=51.1;rlimit[1]=98;
    rlimit[2][0] = 0.51; rlimit[2][1] = 0.98; rlimit[1][0] = 0.6; rlimit[1][1] = 1.12;
    rlimit[0][0] = 0.74; rlimit[0][1] = 1.35; rlimit[3][0] = 0.; rlimit[3][1] = 0.;
    rlimit[4][0] = 0.; rlimit[4][1] = 0.;
    philimit[1][0]=0.006; philimit[1][1]=0.07; philimit[0][0]=0.012;philimit[0][1]=0.127;
    philimit[2][0]=0; philimit[2][1]=0; philimit[3][0]=0; philimit[3][1]=0;
    deltar[1][0]=0.076; deltar[1][1]=0.151; deltar[0][0]=0.121; deltar[0][1]=0.242;
    deltar[2][0]=0.; deltar[2][1]=100; deltar[3][0]=0.;deltar[3][1]=100.;
  }else if (angleflag ==kFAEC && type ==0){
    //small angle 1->2->3->4->5
    rlimit[4][0] = 0.21; rlimit[4][1] = 0.57; rlimit[3][0] = 0.27; rlimit[3][1] = 0.66;
    rlimit[2][0] = 0.33; rlimit[2][1] = 0.78; rlimit[1][0] = 0.43; rlimit[1][1] = 0.95;
    rlimit[0][0] = 0.54; rlimit[0][1] = 1.19;
    layer_array[0] = 5;layer_array[1] = 4;layer_array[2] = 3;layer_array[3] = 2;layer_array[4]=1; nlayer=5;
    philimit[3][0]=0.005; philimit[3][1]=0.066; philimit[2][0]=0.009;philimit[2][1]=0.114;
    philimit[1][0]=0.017; philimit[1][1]=0.173; philimit[0][0]=0.022; philimit[0][1]=0.207;
    deltar[3][0]=0.033; deltar[3][1]=0.088; deltar[2][0]=0.056; deltar[2][1]=0.141;
    deltar[1][0]=0.081; deltar[1][1]=0.201; deltar[0][0]=0.084;deltar[0][1]=0.241;
  }else if (angleflag ==kFAEC && type ==1){
    //small angle 1->3->4->5
    rlimit[3][0] = 0.21; rlimit[3][1] = 0.57; rlimit[2][0] = 0.33; rlimit[2][1] = 0.78;
    rlimit[1][0] = 0.43; rlimit[1][1] = 0.95; rlimit[0][0] = 0.54; rlimit[0][1] = 1.19;
    rlimit[4][0] = 0.; rlimit[4][1] = 0.;
    layer_array[0] = 5;layer_array[1] = 4;layer_array[2] = 3;layer_array[3] = 1;nlayer=4;
    philimit[2][0]=0.015; philimit[2][1]=0.175; philimit[1][0]=0.017;philimit[1][1]=0.173;
    philimit[0][0]=0.022; philimit[0][1]=0.207; philimit[3][0]=0; philimit[3][1]=0; 
    deltar[2][0]=0.091; deltar[2][1]=0.226; deltar[1][0]=0.081; deltar[1][1]=0.201;
    deltar[0][0]=0.084; deltar[0][1]=0.241; deltar[3][0]=0.;deltar[3][1]=100.;
  }else if (angleflag ==kFAEC && type ==2){
    //small angle 1->2->4->5
    rlimit[3][0] = 0.21; rlimit[3][1] = 0.57; rlimit[2][0] = 0.27; rlimit[2][1] = 0.66;
    rlimit[1][0] = 0.43; rlimit[1][1] = 0.95; rlimit[0][0] = 0.54; rlimit[0][1] = 1.19;
    rlimit[4][0] = 0.; rlimit[4][1] = 0.;
    layer_array[0] = 5;layer_array[1] = 4;layer_array[2] = 2;layer_array[3] = 1;nlayer=4;
    philimit[2][0]=0.005; philimit[2][1]=0.066; philimit[1][0]=0.03;philimit[1][1]=0.286;
    philimit[0][0]=0.022; philimit[0][1]=0.207; philimit[3][0]=0; philimit[3][1]=0;
    deltar[2][0]=0.033; deltar[2][1]=0.088; deltar[1][0]=0.141; deltar[1][1]=0.345;
    deltar[0][0]=0.084; deltar[0][1]=0.241; deltar[3][0]=0.;deltar[3][1]=100.;
  }else if (angleflag ==kFAEC && type ==3){
    //small angle 1->2->3->5
    rlimit[3][0] = 0.21; rlimit[3][1] = 0.57; rlimit[2][0] = 0.27; rlimit[2][1] = 0.66;
    rlimit[1][0] = 0.33; rlimit[1][1] = 0.78; rlimit[0][0] = 0.54; rlimit[0][1] = 1.19;
    rlimit[4][0] = 0.; rlimit[4][1] = 0.;
    layer_array[0] = 5;layer_array[1] = 3;layer_array[2] = 2;layer_array[3] = 1;nlayer=4;
    philimit[2][0]=0.005; philimit[2][1]=0.066; philimit[1][0]=0.009;philimit[1][1]=0.114;
    philimit[0][0]=0.042; philimit[0][1]=0.378; philimit[3][0]=0; philimit[3][1]=0;
    deltar[2][0]=0.033; deltar[2][1]=0.088; deltar[1][0]=0.056; deltar[1][1]=0.141;
    deltar[0][0]=0.157; deltar[0][1]=0.439; deltar[3][0]=0.;deltar[3][1]=100.;
  }else if (angleflag ==kFAEC && type ==4){
    //small angle 1->2->3->4
    rlimit[3][0] = 0.21; rlimit[3][1] = 0.57; rlimit[2][0] = 0.27; rlimit[2][1] = 0.66;
    rlimit[1][0] = 0.33; rlimit[1][1] = 0.78; rlimit[0][0] = 0.43; rlimit[0][1] = 0.95;
    rlimit[4][0] = 0.; rlimit[4][1] = 0.;
    layer_array[0] = 4;layer_array[1] = 3;layer_array[2] = 2;layer_array[3] = 1;nlayer=4;
    philimit[2][0]=0.005; philimit[2][1]=0.066; philimit[1][0]=0.009;philimit[1][1]=0.114;
    philimit[0][0]=0.017; philimit[0][1]=0.173; philimit[3][0]=0; philimit[3][1]=0;
    deltar[2][0]=0.033; deltar[2][1]=0.088; deltar[1][0]=0.056; deltar[1][1]=0.141;
    deltar[0][0]=0.081; deltar[0][1]=0.201; deltar[3][0]=0.;deltar[3][1]=100.;
  }else if (angleflag ==kFAEC && type ==5){
    //small angle 2->3->4->5
    rlimit[3][0] = 0.27; rlimit[3][1] = 0.66; rlimit[2][0] = 0.33; rlimit[2][1] = 0.78;
    rlimit[1][0] = 0.43; rlimit[1][1] = 0.95; rlimit[0][0] = 0.54; rlimit[0][1] = 1.19;
    rlimit[4][0] = 0.; rlimit[4][1] = 0.;
    layer_array[0] = 5;layer_array[1] = 4;layer_array[2] = 3;layer_array[3] = 2;nlayer=4;
    philimit[2][0]=0.009; philimit[2][1]=0.114; philimit[1][0]=0.017;philimit[1][1]=0.173;
    philimit[0][0]=0.022; philimit[0][1]=0.207; philimit[3][0]=0; philimit[3][1]=0;
    deltar[2][0]=0.056; deltar[2][1]=0.141; deltar[1][0]=0.081; deltar[1][1]=0.201;
    deltar[0][0]=0.084; deltar[0][1]=0.241; deltar[3][0]=0.;deltar[3][1]=100.;
  }else{
    // large angle 1->2->3->4
    layer_array[0] = 3;layer_array[1] = 2;layer_array[2] = 1;layer_array[3] = 0;nlayer=4;
    rlimit[3][0] = 0.434; rlimit[3][1] = 0.860; rlimit[2][0] = 0.51; rlimit[2][1] = 0.98;
    rlimit[1][0] = 0.60; rlimit[1][1] = 1.12; rlimit[0][0] = 0.74; rlimit[0][1] = 1.35;
    rlimit[4][0] = 0.; rlimit[4][1] = 0.;
    philimit[2][0]=0.003; philimit[2][1]=0.05; philimit[1][0]=0.006;philimit[1][1]=0.07;
    philimit[0][0]=0.012; philimit[0][1]=0.127; philimit[3][0]=0;philimit[3][1]=0;
    deltar[2][0]=0.058; deltar[2][1]=0.122; deltar[1][0]=0.076; deltar[1][1]=0.151;
    deltar[0][0]=0.121; deltar[0][1]=0.242; deltar[3][0]=0.;deltar[3][1]=100.;
  }
  //----------------------------------------------------------------//
  // first layer
  vector<TSeqCollection*> hitArray1 = (*theHitMap).find(layer_array[0])->second;
  for (UInt_t chamber1 = 0; chamber1<hitArray1.size(); chamber1++){
  //cout<<hitArray1[chamber1]->GetLast()+1<<endl;
    for (Int_t layer1 = 0; layer1!=hitArray1[chamber1]->GetLast()+1;layer1++){
      //cout<<hitArray1[chamber1]->GetLast()+1<<endl;
      SoLIDGEMHit *hit1 = (SoLIDGEMHit*)hitArray1[chamber1]->At(layer1);
      //for large angle check the hit with EC
      Is_Seed = 1;
      if (angleflag == kLAEC )
	    {
	       for (UInt_t ec_count=0; ec_count<fCaloHits.size(); ec_count++)
	      {
	      if (fCaloHits[ec_count].fECID != kLAEC) continue; //not LAEC hit
	      Double_t ecHitPhi = TMath::ATan2(fCaloHits[ec_count].fYPos, fCaloHits[ec_count].fXPos);
	      Double_t ecHitR = TMath::Sqrt( TMath::Power(fCaloHits[ec_count].fXPos, 2) + 
	                                     TMath::Power(fCaloHits[ec_count].fYPos, 2) );
	      Double_t tmpDeltaPhi = CalDeltaPhi(ecHitPhi, hit1->GetPhi());
	      Double_t tmpDeltaR   = ecHitR - hit1->GetR();
	      if (hit1->GetTrackerID()==2 && (tmpDeltaPhi < 0.15 && tmpDeltaPhi > -0.05) && (tmpDeltaR < 0.286 && tmpDeltaR > 0.095))
	      {
	      Is_Seed = 1;
      	}
	      else if (hit1->GetTrackerID()==3 && (tmpDeltaPhi<0.06 && tmpDeltaPhi>-0.06) && (tmpDeltaR<0.054 && tmpDeltaR > -0.039) )
	      {
	      Is_Seed = 1;
      	}
	      else
	      {
      	Is_Seed = 0;
       	}
	      }
	    }else if (angleflag == kFAEC){
	       for (UInt_t ec_count=0; ec_count<fCaloHits.size(); ec_count++)
	      {
	      if (fCaloHits[ec_count].fECID != kFAEC) continue; //not FAEC hit
	      Double_t ecHitPhi = TMath::ATan2(fCaloHits[ec_count].fYPos, fCaloHits[ec_count].fXPos);
	      Double_t ecHitR = TMath::Sqrt( TMath::Power(fCaloHits[ec_count].fXPos, 2) + 
	                                     TMath::Power(fCaloHits[ec_count].fYPos, 2) );
	      Double_t tmpDeltaPhi = CalDeltaPhi(ecHitPhi, hit1->GetPhi());
	      Double_t tmpDeltaR   = ecHitR - hit1->GetR();
	      if (hit1->GetTrackerID()==4 && (tmpDeltaPhi < 1.1 && tmpDeltaPhi > 0.04) && (tmpDeltaR < 1.11 && tmpDeltaR > 0.425))
	      {
	      Is_Seed = 1;
      	}
	      else if (hit1->GetTrackerID()==5 && (tmpDeltaPhi<0.9 && tmpDeltaPhi>0.02) && (tmpDeltaR<0.88 && tmpDeltaR > 0.31) )
	      {
	      Is_Seed = 1;
      	}
	      else
	      {
      	Is_Seed = 0;
       	}
	      }
	    }
      if (Is_Seed == 0) continue;
      //checking if the hits is within a good r range
      if (hit1->GetR()<rlimit[0][0] || hit1->GetR() > rlimit[0][1]) continue;
      SoLIDTrack *newtrack = 0;
      if (fDoMC){
#ifdef MCDATA
        newtrack = new ((*fCoarseTracks)[fNTrack++]) SoLIDMCTrack();
#endif
      }
      else{
        newtrack = new ((*fCoarseTracks)[fNTrack++]) SoLIDTrack();
      }
      //second layer
      vector<TSeqCollection*> hitArray2 = (*theHitMap).find(layer_array[1])->second;
      for (UInt_t chamber2 = 0; chamber2<hitArray2.size(); chamber2++) {
	for (Int_t layer2 = 0; layer2!=hitArray2[chamber2]->GetLast()+1;layer2++){
	  SoLIDGEMHit *hit2 = (SoLIDGEMHit*)hitArray2[chamber2]->At(layer2);
	  if (hit2->GetR()<rlimit[1][0] || hit2->GetR() > rlimit[1][1]) continue;
     
	  //crude check
	  Double_t mom_min = 0.6,mom_max = 11; Double_t charge = 0;
	  Double_t theta_min,theta_max;

	  Double_t mom_min_save[4],mom_max_save[4];
	  Double_t theta_min_save[4],theta_max_save[4];

	  if (angleflag==0){
	    theta_min = 14;theta_max = 29;
	  }else{
	    theta_min = 6.;theta_max = 17.5;
	  }
      
	  mom_min_save[0] = mom_min;
	  mom_max_save[0] = mom_max;
	  theta_min_save[0] = theta_min;
	  theta_max_save[0] = theta_max;
	  dr = hit1->GetR()-hit2->GetR();
	  dphi = CalDeltaPhi(hit2->GetPhi(), hit1->GetPhi());
	  if(dr<deltar[0][1] && dr>deltar[0][0] && hit1->GetR()>rlimit[0][0] && hit1->GetR() < rlimit[0][1] && ((dphi >philimit[0][0]&& dphi <philimit[0][1])||(dphi < -1*philimit[0][0]&& dphi > -1*philimit[0][1]))
	     && FindThetaRange(hit2->GetR(),layer_array[1],hit1->GetR(),layer_array[0],&theta_min,&theta_max,angleflag) && 
	     FindMomRange(layer_array[1],hit2->GetPhi(),layer_array[0], hit1->GetPhi(),&mom_min,&mom_max,angleflag)){
	    // check r range
	    // judge the track charge
	    if (dphi >philimit[0][0]&& dphi <philimit[0][1]){
	      charge = 1;
	    }else{
	      charge = -1;
	    }
	    // third layer
	    vector<TSeqCollection*> hitArray3 = (*theHitMap).find(layer_array[2])->second;
	    for (UInt_t chamber3 = 0; chamber3<hitArray3.size(); chamber3++) {
	      for (Int_t layer3 = 0; layer3!=hitArray3[chamber3]->GetLast()+1;layer3++){
		SoLIDGEMHit *hit3 = (SoLIDGEMHit*)hitArray3[chamber3]->At(layer3);
		if (hit3->GetR()<rlimit[2][0] || hit3->GetR() > rlimit[2][1]) continue;
		mom_min_save[1] = mom_min;
		mom_max_save[1] = mom_max;
		theta_min_save[1] = theta_min;
		theta_max_save[1] = theta_max;
		dr = hit2->GetR()-hit3->GetR();
		dphi = CalDeltaPhi(hit3->GetPhi(), hit2->GetPhi());
		if(dr>deltar[1][0] && dr<deltar[1][1] && ((dphi >philimit[1][0]&& dphi <philimit[1][1] && charge==1)||(dphi < -1*philimit[1][0]&& dphi > -1*philimit[1][1]&& charge==-1))
		   &&FindThetaRange(hit3->GetR(),layer_array[2],hit2->GetR(),layer_array[1],&theta_min,&theta_max,angleflag) && 
		   FindMomRange(layer_array[2],hit3->GetPhi(),layer_array[1], hit2->GetPhi(),&mom_min,&mom_max,angleflag)) {
		  // check r range
		  if (nlayer>3){
	       
		    // fourth layer
		    vector<TSeqCollection*> hitArray4 = (*theHitMap).find(layer_array[3])->second;
		    for (UInt_t chamber4 = 0; chamber4<hitArray4.size(); chamber4++) {
		      for (Int_t layer4 = 0; layer4!=hitArray4[chamber4]->GetLast()+1;layer4++){
			SoLIDGEMHit *hit4 = (SoLIDGEMHit*)hitArray4[chamber4]->At(layer4);
			if (hit4->GetR()<rlimit[3][0] || hit4->GetR() > rlimit[3][1]) continue;
		
			mom_min_save[2] = mom_min;
			mom_max_save[2] = mom_max;
			theta_min_save[2] = theta_min;
			theta_max_save[2] = theta_max;
			dr = hit3->GetR()-hit4->GetR();
			dphi = CalDeltaPhi(hit4->GetPhi(), hit3->GetPhi());
			if(dr>deltar[2][0] && dr<deltar[2][1] && ((dphi >philimit[2][0]&& dphi <philimit[2][1] && charge==1)||(dphi < -1*philimit[2][0]&& dphi > -1*philimit[2][1]&& charge==-1))
			   && FindThetaRange(hit4->GetR(),layer_array[3],hit3->GetR(),layer_array[2],&theta_min,&theta_max,angleflag) && 
			   FindMomRange(layer_array[3],hit4->GetPhi(),layer_array[2], hit3->GetPhi(),&mom_min,&mom_max,angleflag)){
			  if (nlayer>4){
		    
			    // fifth layer
			    vector<TSeqCollection*> hitArray5 = (*theHitMap).find(layer_array[4])->second;
			    for (UInt_t chamber5 = 0; chamber5<hitArray5.size(); chamber5++) {
			      for (Int_t layer5 = 0; layer5!=hitArray5[chamber5]->GetLast()+1;layer5++){
				SoLIDGEMHit *hit5 = (SoLIDGEMHit*)hitArray5[chamber5]->At(layer5);
				if (hit5->GetR()<rlimit[4][0] || hit5->GetR() > rlimit[4][1]) continue;
		      
				mom_min_save[3] = mom_min;
				mom_max_save[3] = mom_max;
				theta_min_save[3] = theta_min;
				theta_max_save[3] = theta_max;
				dr = hit4->GetR() - hit5->GetR();
				dphi = CalDeltaPhi(hit5->GetPhi(), hit4->GetPhi());
				if(dr>deltar[3][0] && dr<deltar[3][1] && ((dphi >philimit[3][0]&& dphi <philimit[3][1] && charge==1)||(dphi < -1*philimit[3][0]&& dphi > -1*philimit[3][1]&& charge==-1))
				   &&FindThetaRange(hit5->GetR(),layer_array[4],hit4->GetR(),layer_array[3],&theta_min,&theta_max,angleflag) && 
				   FindMomRange(layer_array[4],hit5->GetPhi(),layer_array[3], hit4->GetPhi(),&mom_min,&mom_max,angleflag)){
				  newtrack->AddHit(hit5);
				  newtrack->AddHit(hit4);
				  newtrack->AddHit(hit3);
				  newtrack->AddHit(hit2);
				  newtrack->AddHit(hit1);
				  newtrack->SetCharge(charge);
				  newtrack->SetAngleFlag(angleflag);
				  newtrack->SetMomMin(mom_min);
				  newtrack->SetMomMax(mom_max);
				  newtrack->SetThetaMin(theta_min);
				  newtrack->SetThetaMax(theta_max);
				  if (fDoMC){
#ifdef MCDATA
				    newtrack = new ((*fCoarseTracks)[fNTrack++]) SoLIDMCTrack();
#endif
				  }
				  else{
				    newtrack = new ((*fCoarseTracks)[fNTrack++]) SoLIDTrack();  
				  }
				  newtrack->Clear();
				}else{
				  newtrack->Clear();
				  mom_min = mom_min_save[3];
				  mom_max = mom_max_save[3];
				  theta_min = theta_min_save[3];
				  theta_max = theta_max_save[3];
				}
			      }//loop on all hits on chamber on fifth tracker 
			    }//loop on all chamber on fifth tracker
			  }else{
			    if (nlayer==4){
			      newtrack->AddHit(hit4);
			      newtrack->AddHit(hit3);
			      newtrack->AddHit(hit2);
			      newtrack->AddHit(hit1);
			      newtrack->SetCharge(charge);
			      newtrack->SetAngleFlag(angleflag);
			      newtrack->SetMomMin(mom_min);
			      newtrack->SetMomMax(mom_max);
			      newtrack->SetThetaMin(theta_min);
			      newtrack->SetThetaMax(theta_max);
			      if (fDoMC){
#ifdef MCDATA
				    newtrack = new ((*fCoarseTracks)[fNTrack++]) SoLIDMCTrack();
#endif
				  }
				  else{
				    newtrack = new ((*fCoarseTracks)[fNTrack++]) SoLIDTrack();  
				  }
			      newtrack->Clear();
			    }
			  }
			}else{
			  mom_min = mom_min_save[2];
			  mom_max = mom_max_save[2];
			  theta_min = theta_min_save[2];
			  theta_max = theta_max_save[2];
			}
		      }//loop on all his on chamber on fourth tracker
		    }//loop on all chamber on fourth tracker
		  }else{
		    if (nlayer==3){
		      newtrack->AddHit(hit3);
		      newtrack->AddHit(hit2);
		      newtrack->AddHit(hit1);
		      newtrack->SetMomMin(mom_min);
		      newtrack->SetMomMax(mom_max);
		      newtrack->SetThetaMin(theta_min);
		      newtrack->SetThetaMax(theta_max);
		      newtrack->SetCharge(charge);
		      newtrack->SetAngleFlag(angleflag);
		      if (fDoMC){
#ifdef MCDATA
				    newtrack = new ((*fCoarseTracks)[fNTrack++]) SoLIDMCTrack();
#endif
				  }
				  else{
				    newtrack = new ((*fCoarseTracks)[fNTrack++]) SoLIDTrack();  
				  }
		      newtrack->Clear();
		    }
		  }
		}else{
		  mom_min = mom_min_save[1];
		  mom_max = mom_max_save[1];
		  theta_min = theta_min_save[1];
		  theta_max = theta_max_save[1];
		}
	      }//loop on all hits on chamber on third tracker
	    }//loop on all chamber on third tracker
	  }else{
	    mom_min = mom_min_save[0];
	    mom_max = mom_max_save[0];
	    theta_min = theta_min_save[0];
	    theta_max = theta_max_save[0];
	  }
	}//loop on all hits on chamber on second tracker
      }//loop on all chamber on second tracker
      newtrack->Clear();
      fCoarseTracks->Remove(newtrack);
      fNTrack--;
    }//loop on all hits on chamber on first tracker
  }//loop on all chamber on first tracker
}
//________________________________________________________________________________
void ProgressiveTracking::CheckTracks()
{
  Int_t nmom;
  Int_t c_count = 0;
  Double_t q_asy = -1;
  Double_t Rmom[20]={0},Xc[20]={0},Yc[20]={0},Rc[20]={0};
  Double_t reconVertexZ = 0.;
  Double_t reconTheta = 0.;
  for (Int_t i=0;i!=fCoarseTracks->GetLast()+1;i++){
    SoLIDTrack *thistrack = (SoLIDTrack*)(fCoarseTracks->At(i));
    if (thistrack->GetAngleFlag()==kLAEC && thistrack->GetNHits()==3){
      nmom=0;
      c_count = 0;
      for (UInt_t j=0; j<thistrack->GetNHits(); j++){
	q_asy = CalChargeAsy(thistrack->GetHit(j)->GetQU(), thistrack->GetHit(j)->GetQV());
	if (TMath::Abs(q_asy)<0.25) {c_count++;}
      }
      if (c_count<3)
	{
	    thistrack->SetCoarseFitStatus(kTRUE); //not really a fit yet
	    thistrack->SetStatus(kFALSE);
	    thistrack->SetCoarseChi2(kINFINITY);
	  continue;
	}
	     
      SoLIDGEMHit *hit1 = (SoLIDGEMHit*)(thistrack->GetHit(0));
      SoLIDGEMHit *hit2 = (SoLIDGEMHit*)(thistrack->GetHit(1));
      SoLIDGEMHit *hit3 = (SoLIDGEMHit*)(thistrack->GetHit(2));
      Int_t flag = 1;
      
      if (flag == 1)
	{
	  CalCircle(hit1->GetX(),hit1->GetY(),hit2->GetX(),hit2->GetY(),hit3->GetX(),hit3->GetY(),Rmom,Xc,Yc);
	  Rc[0] = sqrt(Xc[0]*Xc[0]+Yc[0]*Yc[0]);
	  if (hit1->GetTrackerID()==0 && hit2->GetTrackerID()==1&&hit3->GetTrackerID()==2){
	    Int_t j=0;
	    if ( (Rmom[j]-(-0.0712984+1.01393*Rc[j]))<0.1 && (Rmom[j]-(-0.0712984+1.01393*Rc[j]))>-0.2 
	         && CalVertexZ(hit1, hit2, hit3, &Rmom[j], kLAEC, &reconVertexZ, &reconTheta) ){
	      
	      thistrack->SetStatus(kTRUE);
	      thistrack->SetCoarseFitStatus(kTRUE); //not really a fit yet
	      thistrack->SetCoarseChi2(fabs(Rmom[0]+0.0712984-1.01393*Rc[0]));
	      
	    }else{
	      flag = 0;
	    }
	       
	  }
	  if (hit1->GetTrackerID()==0 && hit2->GetTrackerID()==1&&hit3->GetTrackerID()==3){
	    Int_t j=0;
	    if ( (Rmom[j]-(-0.0719744+1.01388*Rc[j]))<0.1 && (Rmom[j]-(-0.0719744+1.01388*Rc[j]))>-0.2 
	         && CalVertexZ(hit1, hit2, hit3, &Rmom[j], kLAEC, &reconVertexZ, &reconTheta) ){
	      thistrack->SetStatus(kTRUE);
	      thistrack->SetCoarseFitStatus(kTRUE); //not really a fit yet
	      thistrack->SetCoarseChi2(fabs(Rmom[0]+0.0719744-1.01388*Rc[0]));
	    }else{
	      flag = 0;
	    }
	  }
	  if (hit1->GetTrackerID()==0 && hit2->GetTrackerID()==2&&hit3->GetTrackerID()==3){
	    Int_t j=0;
	    if ( (Rmom[j]-(-0.0737444+1.01419*Rc[j]))<0.1 && (Rmom[j]-(-0.0737444+1.01419*Rc[j]))>-0.2 
	         && CalVertexZ(hit1, hit2, hit3, &Rmom[j], kLAEC, &reconVertexZ, &reconTheta) ){
	      thistrack->SetStatus(kTRUE);
	      thistrack->SetCoarseFitStatus(kTRUE); //not really a fit yet
	      thistrack->SetCoarseChi2(fabs(Rmom[0]+0.0737444-1.01419*Rc[0]));
	    }else{
	      flag = 0;
	    }
	  }
	  if (hit1->GetTrackerID()==1 && hit2->GetTrackerID()==2&&hit3->GetTrackerID()==3){
	    Int_t j=0;
	    if ( (Rmom[j]-(-0.0765712+1.01482*Rc[j]))<0.1 && (Rmom[j]-(-0.0765712+1.01482*Rc[j]))>-0.2 
	         && CalVertexZ(hit1, hit2, hit3, &Rmom[j], kLAEC, &reconVertexZ, &reconTheta) ){
	      thistrack->SetStatus(kTRUE);
	      thistrack->SetCoarseFitStatus(kTRUE); //not really a fit yet
	      thistrack->SetCoarseChi2(fabs(Rmom[0]+0.0765712-1.01482*Rc[0]));
	    }else{
	      flag = 0;
	    }
	  }
	}
      
      if (flag == 1){
      thistrack->SetVertexZ(reconVertexZ);
      thistrack->SetTheta(reconTheta);
      }else{
      thistrack->SetCoarseFitStatus(kTRUE); //not really a fit yet
	    thistrack->SetStatus(kFALSE);
	    thistrack->SetCoarseChi2(kINFINITY);
      }
    }else if (thistrack->GetAngleFlag()==kLAEC && thistrack->GetNHits()==4){
      nmom  = 0;
      c_count = 0;
      for (UInt_t j=0; j<thistrack->GetNHits(); j++){
	      q_asy = CalChargeAsy(thistrack->GetHit(j)->GetQU(), thistrack->GetHit(j)->GetQV());
	      if (TMath::Abs(q_asy)<0.25) {c_count++;}
      }
      if (c_count<3)
	    {
	     thistrack->SetCoarseFitStatus(kTRUE); //not really a fit yet
	    thistrack->SetStatus(kFALSE);
	    thistrack->SetCoarseChi2(kINFINITY);
	     continue;
	    }
      SoLIDGEMHit *hit1 = (SoLIDGEMHit*)(thistrack->GetHit(0));
      SoLIDGEMHit *hit2 = (SoLIDGEMHit*)(thistrack->GetHit(1));
      SoLIDGEMHit *hit3 = (SoLIDGEMHit*)(thistrack->GetHit(2));
      SoLIDGEMHit *hit4 = (SoLIDGEMHit*)(thistrack->GetHit(3));
      
      Int_t flag = 1;
      if (flag == 1)
	{
	  CalCircle(hit1->GetX(),hit1->GetY(),hit2->GetX(),hit2->GetY(),hit3->GetX(),hit3->GetY(),&Rmom[nmom],&Xc[nmom],&Yc[nmom]);
	  Rc[nmom] = sqrt(Xc[nmom]*Xc[nmom]+Yc[nmom]*Yc[nmom]);
	  nmom++;
	  CalCircle(hit1->GetX(),hit1->GetY(),hit2->GetX(),hit2->GetY(),hit4->GetX(),hit4->GetY(),&Rmom[nmom],&Xc[nmom],&Yc[nmom]);
	  Rc[nmom] = sqrt(Xc[nmom]*Xc[nmom]+Yc[nmom]*Yc[nmom]);
	  nmom++;
	  CalCircle(hit1->GetX(),hit1->GetY(),hit3->GetX(),hit3->GetY(),hit4->GetX(),hit4->GetY(),&Rmom[nmom],&Xc[nmom],&Yc[nmom]);
	  Rc[nmom] = sqrt(Xc[nmom]*Xc[nmom]+Yc[nmom]*Yc[nmom]);
	  nmom++;
	  CalCircle(hit2->GetX(),hit2->GetY(),hit3->GetX(),hit3->GetY(),hit4->GetX(),hit4->GetY(),&Rmom[nmom],&Xc[nmom],&Yc[nmom]);
	  Rc[nmom] = sqrt(Xc[nmom]*Xc[nmom]+Yc[nmom]*Yc[nmom]);
	  nmom++;
	     
	  //compare single
	  //012
	  Double_t vertexzAvg = 0.;
	  Double_t thetaAvg = 0.;
	  Int_t countAvg = 0;
	  Int_t j=0;
	  if ( (Rmom[j]-(-0.0712984+1.01393*Rc[j]))<0.1 && (Rmom[j]-(-0.0712984+1.01393*Rc[j]))>-0.2
	        && CalVertexZ(hit1, hit2, hit3, &Rmom[j], kLAEC, &reconVertexZ, &reconTheta) ){
	        vertexzAvg += reconVertexZ;
	        thetaAvg += reconTheta;
	        countAvg++;
	  }else{
	    flag = 0;
	  }
	  //013
	  j=1;
	  if ( (Rmom[j]-(-0.0719744+1.01388*Rc[j]))<0.1 && (Rmom[j]-(-0.0719744+1.01388*Rc[j]))>-0.2 
	        && CalVertexZ(hit1, hit2, hit4, &Rmom[j], kLAEC, &reconVertexZ, &reconTheta) ){
	        vertexzAvg += reconVertexZ;
	        thetaAvg += reconTheta;
	        countAvg++;
	  }else{
	    flag = 0;
	  }
	  //023
	  j=2;
	  if ( (Rmom[j]-(-0.0737444+1.01419*Rc[j]))<0.1 && (Rmom[j]-(-0.0737444+1.01419*Rc[j]))>-0.2 
	      && CalVertexZ(hit1, hit3, hit4, &Rmom[j], kLAEC, &reconVertexZ, &reconTheta) ){
	      vertexzAvg += reconVertexZ;
	      thetaAvg += reconTheta;
	      countAvg++;
	  }else{
	    flag = 0;
	  }
	  //123
	  j=3;
	  if ( (Rmom[j]-(-0.0765712+1.01482*Rc[j]))<0.1 && (Rmom[j]-(-0.0765712+1.01482*Rc[j]))>-0.2 
	       && CalVertexZ(hit2, hit3, hit4, &Rmom[j], kLAEC, &reconVertexZ, &reconTheta) ){
	       vertexzAvg += reconVertexZ;
	       thetaAvg += reconTheta;
	       countAvg++;
	  }else{
	    flag = 0;
	  }
	      
	  //compare between 2
	  if ((Rmom[3]-Rmom[2]-(-0.000907009+1.02995*(Rc[3]-Rc[2])))>-0.08&&(Rmom[3]-Rmom[2]-(-0.000907009+1.02995*(Rc[3]-Rc[2])))<0.08 && sqrt((Xc[3]-Xc[2])*(Xc[3]-Xc[2])+(Yc[3]-Yc[2])*(Yc[3]-Yc[2]))<0.8*sqrt(2)){
	  }else{
	    flag = 0;
	  }
	     
	  if (sqrt((Xc[2]-Xc[1])*(Xc[2]-Xc[1])+(Yc[2]-Yc[1])*(Yc[2]-Yc[1]))<0.8*sqrt(2)&& (Rmom[2]-Rmom[1]-(-0.000611265+1.02364*(Rc[2]-Rc[1])))>-0.08&&(Rmom[2]-Rmom[1]-(-0.000611265+1.02364*(Rc[2]-Rc[1])))<0.08){
	  }else{
	    flag = 0;
	  }
	  if (flag == 1) {
	    assert(countAvg == 4);
	    reconTheta = thetaAvg/4.;
	    reconVertexZ = vertexzAvg/4.;
	  }
	}
	 
	 
      if (flag==1){
	     thistrack->SetStatus(kTRUE);
	     thistrack->SetCoarseChi2((sqrt((Xc[1]-Xc[0])*(Xc[1]-Xc[0])+(Yc[1]-Yc[0])*(Yc[1]-Yc[0]))+sqrt((Xc[2]-Xc[1])*(Xc[2]-Xc[1])+(Yc[2]-Yc[1])*(Yc[2]-Yc[1]))+sqrt((Xc[3]-Xc[2])*(Xc[3]-Xc[2])+(Yc[3]-Yc[2])*(Yc[3]-Yc[2]))+sqrt((Xc[3]-Xc[0])*(Xc[3]-Xc[0])+(Yc[3]-Yc[0])*(Yc[3]-Yc[0])))/4.);
	     thistrack->SetCoarseFitStatus(kTRUE);
	     thistrack->SetVertexZ(reconVertexZ);
	     thistrack->SetTheta(reconTheta);
      }else{
	     thistrack->SetCoarseFitStatus(kTRUE); //not really a fit yet
	    thistrack->SetStatus(kFALSE);
	    thistrack->SetCoarseChi2(kINFINITY);
      }
    }else if (thistrack->GetAngleFlag()==kFAEC && thistrack->GetNHits()==4){
      c_count = 0;
       for (UInt_t j=0; j<thistrack->GetNHits(); j++){
	       q_asy = CalChargeAsy(thistrack->GetHit(j)->GetQU(), thistrack->GetHit(j)->GetQV());
	       if (TMath::Abs(q_asy)<0.25) {c_count++;}
      }

      if (c_count<3)
	    {
	      thistrack->SetCoarseFitStatus(kTRUE); //not really a fit yet
	      thistrack->SetStatus(kFALSE);
	      thistrack->SetCoarseChi2(kINFINITY);
	      continue;
	    }
      SoLIDGEMHit *hit1 = (SoLIDGEMHit*)(thistrack->GetHit(0));
      SoLIDGEMHit *hit2 = (SoLIDGEMHit*)(thistrack->GetHit(1));
      SoLIDGEMHit *hit3 = (SoLIDGEMHit*)(thistrack->GetHit(2));
      SoLIDGEMHit *hit4 = (SoLIDGEMHit*)(thistrack->GetHit(3));
      Int_t flag = 1;
      
 
      if (flag == 1)
	{
	  nmom  = 0;
	  CalCircle(hit1->GetX(),hit1->GetY(),hit2->GetX(),hit2->GetY(),hit3->GetX(),hit3->GetY(),&Rmom[nmom],&Xc[nmom],&Yc[nmom]);
	  Rc[nmom] = sqrt(Xc[nmom]*Xc[nmom]+Yc[nmom]*Yc[nmom]);
	  nmom++;
	  CalCircle(hit1->GetX(),hit1->GetY(),hit2->GetX(),hit2->GetY(),hit4->GetX(),hit4->GetY(),&Rmom[nmom],&Xc[nmom],&Yc[nmom]);
	  Rc[nmom] = sqrt(Xc[nmom]*Xc[nmom]+Yc[nmom]*Yc[nmom]);
	  nmom++;
	  CalCircle(hit1->GetX(),hit1->GetY(),hit3->GetX(),hit3->GetY(),hit4->GetX(),hit4->GetY(),&Rmom[nmom],&Xc[nmom],&Yc[nmom]);
	  Rc[nmom] = sqrt(Xc[nmom]*Xc[nmom]+Yc[nmom]*Yc[nmom]);
	  nmom++;
	  CalCircle(hit2->GetX(),hit2->GetY(),hit3->GetX(),hit3->GetY(),hit4->GetX(),hit4->GetY(),&Rmom[nmom],&Xc[nmom],&Yc[nmom]);
	  Rc[nmom] = sqrt(Xc[nmom]*Xc[nmom]+Yc[nmom]*Yc[nmom]);
	  nmom++;
	  Double_t vertexzAvg = 0.;
	  Double_t thetaAvg = 0.;
	  Int_t countAvg = 0;
	  //1234
	  if (hit1->GetTrackerID()==1 && hit2->GetTrackerID()==2 && hit3->GetTrackerID()==3 && hit4->GetTrackerID()==4){
	    //123
	    if ( Rmom[0]-(-0.0374431+1.0118*Rc[0])< 0.06 && Rmom[0]-(-0.0374431+1.0118*Rc[0])>-0.1
	       && CalVertexZ(hit1, hit2, hit3, &Rmom[0], kFAEC, &reconVertexZ, &reconTheta) ){
	       vertexzAvg += reconVertexZ;
	      thetaAvg += reconTheta;
	      countAvg++;
	    }else{
	      flag = 0;
	    }
	    //124
	    if ( Rmom[1]-(-0.0345228+1.01045*Rc[1])< 0.06 && Rmom[1]-(-0.0345228+1.01045*Rc[1])>-0.1
	       && CalVertexZ(hit1, hit2, hit4, &Rmom[1], kFAEC, &reconVertexZ, &reconTheta) ){
	       vertexzAvg += reconVertexZ;
	      thetaAvg += reconTheta;
	      countAvg++;
	    }else{
	      flag = 0;
	    }
	    //134
	    if ( Rmom[2]-(-0.032753+1.0098*Rc[2])< 0.06 && Rmom[2]-(-0.032753+1.0098*Rc[2])>-0.1
	       && CalVertexZ(hit1, hit3, hit4, &Rmom[2], kFAEC, &reconVertexZ, &reconTheta) ){
	       vertexzAvg += reconVertexZ;
	      thetaAvg += reconTheta;
	      countAvg++;
	    }else{
	      flag = 0;
	    }
	    //234
	    if (Rmom[3]-(-0.0316562+1.00954*Rc[3])<0.06 && Rmom[3]-(-0.0316562+1.00954*Rc[3])>-0.1
	       && CalVertexZ(hit2, hit3, hit4, &Rmom[3], kFAEC, &reconVertexZ, &reconTheta) ){
	       vertexzAvg += reconVertexZ;
	      thetaAvg += reconTheta;
	      countAvg++;
	    }else{
	      flag = 0;
	    }
	    //cross-check
	    if (sqrt((Xc[3]-Xc[2])*(Xc[3]-Xc[2])+(Yc[3]-Yc[2])*(Yc[3]-Yc[2]))<0.4*sqrt(2) && Rmom[3]-Rmom[2]-(-0.000290564+1.04145*(Rc[3]-Rc[2]))<0.05 && Rmom[3]-Rmom[2]-(0.000290564+1.04145*(Rc[3]-Rc[2]))>-0.05){
	    }else{
	      flag = 0;
	    }
	  }
	   
	  //1235
	  if (hit1->GetTrackerID()==1 && hit2->GetTrackerID()==2 && hit3->GetTrackerID()==3 && hit4->GetTrackerID()==5){
	    //123
	    if ( Rmom[0]-(-0.0374431+1.0118*Rc[0])< 0.06&& Rmom[0]-(-0.0374431+1.0118*Rc[0])>-0.1
	       && CalVertexZ(hit1, hit2, hit3, &Rmom[0], kFAEC, &reconVertexZ, &reconTheta) ){
	        vertexzAvg += reconVertexZ;
	      thetaAvg += reconTheta;
	      countAvg++;
	    }else{
	      flag = 0;
	    }
	    //125
	    if ( Rmom[1]-(-0.0315397+1.00942*Rc[1])< 0.06&& Rmom[1]-(-0.0315397+1.00942*Rc[1])>-0.1
	       && CalVertexZ(hit1, hit2, hit4, &Rmom[1], kFAEC, &reconVertexZ, &reconTheta) ){
	       vertexzAvg += reconVertexZ;
	      thetaAvg += reconTheta;
	      countAvg++;
	    }else{
	      flag = 0;
	    }
	    //135
	    if ( Rmom[2]-(-0.028547+1.00847*Rc[2])<0.06 && Rmom[2]-(-0.028547+1.00847*Rc[2])>-0.1
	       && CalVertexZ(hit1, hit3, hit4, &Rmom[2], kFAEC, &reconVertexZ, &reconTheta) ){
	       vertexzAvg += reconVertexZ;
	      thetaAvg += reconTheta;
	      countAvg++;
	    }else{
	      flag = 0;
	    }
	    //235
	    if ( Rmom[3]-(-0.0259752+1.00773*Rc[3])< 0.06&& Rmom[3]-(-0.0259752+1.00773*Rc[3])>-0.1
	       && CalVertexZ(hit2, hit3, hit4, &Rmom[3], kFAEC, &reconVertexZ, &reconTheta) ){
	       vertexzAvg += reconVertexZ;
	      thetaAvg += reconTheta;
	      countAvg++;
	    }else{
	      flag = 0;
	    }
	    if (sqrt((Xc[3]-Xc[2])*(Xc[3]-Xc[2])+(Yc[3]-Yc[2])*(Yc[3]-Yc[2]))<0.3*sqrt(2) && Rmom[3]-Rmom[2]-(0.000678468+1.0465*(Rc[3]-Rc[2]))<0.05 && Rmom[3]-Rmom[2]-(0.000678468+1.0465*(Rc[3]-Rc[2]))>-0.05){
	    }else{
	      flag = 0;
	    }
	  }
	   
	  //1245
	  if (hit1->GetTrackerID()==1 && hit2->GetTrackerID()==2 && hit3->GetTrackerID()==4 && hit4->GetTrackerID()==5){
	    //124
	    if ( Rmom[0]-(-0.0345228+1.01045*Rc[0])< 0.06&& Rmom[0]-(-0.0345228+1.01045*Rc[0])>-0.1
	       && CalVertexZ(hit1, hit2, hit3, &Rmom[0], kFAEC, &reconVertexZ, &reconTheta) ){
	       vertexzAvg += reconVertexZ;
	      thetaAvg += reconTheta;
	      countAvg++;
	    }else{
	      flag = 0;
	    }
	    //125
	    if ( Rmom[1]-(-0.0315397+1.00942*Rc[1])< 0.06&& Rmom[1]-(-0.0315397+1.00942*Rc[1])>-0.1
	       && CalVertexZ(hit1, hit2, hit4, &Rmom[1], kFAEC, &reconVertexZ, &reconTheta) ){
	       vertexzAvg += reconVertexZ;
	      thetaAvg += reconTheta;
	      countAvg++;
	    }else{
	      flag = 0;
	    }
	    //145
	    if ( Rmom[2]-(-0.0227897+1.00688*Rc[2])<0.06 && Rmom[2]-(-0.0227897+1.00688*Rc[2])>-0.15
	       && CalVertexZ(hit1, hit3, hit4, &Rmom[2], kFAEC, &reconVertexZ, &reconTheta) ){
	       vertexzAvg += reconVertexZ;
	      thetaAvg += reconTheta;
	      countAvg++;
	    }else{
	      flag = 0;
	    }
	    //245
	    if ( Rmom[3]-(-0.0190085+1.00589*Rc[3])<0.06 && Rmom[3]-(-0.0190085+1.00589*Rc[3])>-0.15
	       && CalVertexZ(hit2, hit3, hit4, &Rmom[3], kFAEC, &reconVertexZ, &reconTheta) ){
	       vertexzAvg += reconVertexZ;
	      thetaAvg += reconTheta;
	      countAvg++;
	    }else{
	      flag = 0;
	    }
	    if (sqrt((Xc[3]-Xc[2])*(Xc[3]-Xc[2])+(Yc[3]-Yc[2])*(Yc[3]-Yc[2]))<0.16*sqrt(2) && Rmom[3]-Rmom[2]-(-0.00134498+1.04092*(Rc[3]-Rc[2]))< 0.05&& Rmom[3]-Rmom[2]-(-0.00134498+1.04092*(Rc[3]-Rc[2]))>-0.05){
	    }else{
	      flag = 0;
	    }
	  }
	   
	   
	  //1345
	  if (hit1->GetTrackerID()==1 && hit2->GetTrackerID()==3 && hit3->GetTrackerID()==4 && hit4->GetTrackerID()==5){
	    //134
	    if ( Rmom[0]-(-0.032753+1.0098*Rc[0])< 0.06&& Rmom[0]-(-0.032753+1.0098*Rc[0])>-0.1
	       && CalVertexZ(hit1, hit2, hit3, &Rmom[0], kFAEC, &reconVertexZ, &reconTheta) ){
	       vertexzAvg += reconVertexZ;
	      thetaAvg += reconTheta;
	      countAvg++;
	    }else{
	      flag = 0;
	    }
	    //135
	    if ( Rmom[1]-(-0.028547+1.00847*Rc[1])<0.06 && Rmom[1]-(-0.028547+1.00847*Rc[1])>-0.1
	       && CalVertexZ(hit1, hit2, hit4, &Rmom[1], kFAEC, &reconVertexZ, &reconTheta) ){
	       vertexzAvg += reconVertexZ;
	      thetaAvg += reconTheta;
	      countAvg++;
	    }else{
	      flag = 0;
	    }
	    //145
	    if ( Rmom[2]-(-0.0227897+1.00688*Rc[2])<0.06 && Rmom[2]-(-0.0227897+1.00688*Rc[2])>-0.15
	       && CalVertexZ(hit1, hit3, hit4, &Rmom[2], kFAEC, &reconVertexZ, &reconTheta) ){
	       vertexzAvg += reconVertexZ;
	      thetaAvg += reconTheta;
	      countAvg++;
	    }else{
	      flag = 0;
	    }
	    //345
	    if ( Rmom[3]-(-0.0116765+1.00414*Rc[3])<0.07 && Rmom[3]-(-0.0116765+1.00414*Rc[3])>-0.2
	       && CalVertexZ(hit2, hit3, hit4, &Rmom[3], kFAEC, &reconVertexZ, &reconTheta) ){
	       vertexzAvg += reconVertexZ;
	      thetaAvg += reconTheta;
	      countAvg++;
	    }else{
	      flag = 0;
	    }
	    if (sqrt((Xc[2]-Xc[1])*(Xc[2]-Xc[1])+(Yc[2]-Yc[1])*(Yc[2]-Yc[1]))<0.4*sqrt(2) && Rmom[2]-Rmom[1]-(-0.00156074+1.039*(Rc[2]-Rc[1]))< 0.05 && Rmom[2]-Rmom[1]-(-0.00156074+1.039*(Rc[2]-Rc[1]))>-0.1){
	    }else{
	      flag = 0;
	    }
	  }
	   
	  //2345
	  if (hit1->GetTrackerID()==2 && hit2->GetTrackerID()==3 && hit3->GetTrackerID()==4 && hit4->GetTrackerID()==5){
	    //234
	    if ( Rmom[0]-(-0.0316562+1.00954*Rc[0])<0.06 && Rmom[0]-(-0.0316562+1.00954*Rc[0])>-0.1
	       && CalVertexZ(hit1, hit2, hit3, &Rmom[0], kFAEC, &reconVertexZ, &reconTheta) ){
	       vertexzAvg += reconVertexZ;
	      thetaAvg += reconTheta;
	      countAvg++;
	    }else{
	      flag = 0;
	    }
	    //235
	    if ( Rmom[1]-(-0.0259752+1.00773*Rc[1])< 0.06&& Rmom[1]-(-0.0259752+1.00773*Rc[1])>-0.1
	       && CalVertexZ(hit1, hit2, hit4, &Rmom[1], kFAEC, &reconVertexZ, &reconTheta) ){
	       vertexzAvg += reconVertexZ;
	      thetaAvg += reconTheta;
	      countAvg++;
	    }else{
	      flag = 0;
	    }
	    //245
	    if ( Rmom[2]-(-0.0190085+1.00589*Rc[2])<0.06 && Rmom[2]-(-0.0190085+1.00589*Rc[2])>-0.15
	       && CalVertexZ(hit1, hit3, hit4, &Rmom[2], kFAEC, &reconVertexZ, &reconTheta) ){
	       vertexzAvg += reconVertexZ;
	      thetaAvg += reconTheta;
	      countAvg++;
	    }else{
	      flag = 0;
	    }
	    //345
	    if ( Rmom[3]-(-0.0116765+1.00414*Rc[3])<0.07 && Rmom[3]-(-0.0116765+1.00414*Rc[3])>-0.2
	       && CalVertexZ(hit2, hit3, hit4, &Rmom[3], kFAEC, &reconVertexZ, &reconTheta) ){
	       vertexzAvg += reconVertexZ;
	      thetaAvg += reconTheta;
	      countAvg++;
	    }else{
	      flag = 0;
	    }
	    if (sqrt((Xc[2]-Xc[1])*(Xc[2]-Xc[1])+(Yc[2]-Yc[1])*(Yc[2]-Yc[1]))<0.4*sqrt(2) && Rmom[2]-Rmom[1]-(-0.00196282+1.04375*(Rc[2]-Rc[1]))<0.06 && Rmom[2]-Rmom[1]-(-0.00196282+1.04375*(Rc[2]-Rc[1]))>-0.08){
	    }else{
	      flag = 0;
	    }
	  }
	  if (flag == 1) {
	    assert(countAvg == 4);
	    reconTheta = thetaAvg/4.;
	    reconVertexZ = vertexzAvg/4.;
	  }
	  
	}      
       
       
      if (flag==1){
	      thistrack->SetStatus(kTRUE);
	      thistrack->SetCoarseChi2((sqrt((Xc[1]-Xc[0])*(Xc[1]-Xc[0])+(Yc[1]-Yc[0])*(Yc[1]-Yc[0]))+sqrt((Xc[2]-Xc[1])*(Xc[2]-Xc[1])+(Yc[2]-Yc[1])*(Yc[2]-Yc[1]))+sqrt((Xc[3]-Xc[2])*(Xc[3]-Xc[2])+(Yc[3]-Yc[2])*(Yc[3]-Yc[2]))+sqrt((Xc[3]-Xc[0])*(Xc[3]-Xc[0])+(Yc[3]-Yc[0])*(Yc[3]-Yc[0])))/4.);
	      thistrack->SetCoarseFitStatus(kTRUE);
	      thistrack->SetVertexZ(reconVertexZ);
	      thistrack->SetTheta(reconTheta);
      }else{
        thistrack->SetCoarseFitStatus(kTRUE);
	      thistrack->SetStatus(kFALSE);
	      thistrack->SetCoarseChi2(kINFINITY);
      }
    }else if (thistrack->GetAngleFlag()==1 && thistrack->GetNHits()==5){
      c_count = 0;
       for (UInt_t j=0; j<thistrack->GetNHits(); j++){
	       q_asy = CalChargeAsy(thistrack->GetHit(j)->GetQU(), thistrack->GetHit(j)->GetQV());
	       if (TMath::Abs(q_asy)<0.25) {c_count++;}
      }
      if (c_count<3)
	{
	   thistrack->SetCoarseFitStatus(kTRUE); //not really a fit yet
	    thistrack->SetStatus(kFALSE);
	    thistrack->SetCoarseChi2(kINFINITY);
	  continue;
	}
      SoLIDGEMHit *hit1 = (SoLIDGEMHit*)(thistrack->GetHit(0));
      SoLIDGEMHit *hit2 = (SoLIDGEMHit*)(thistrack->GetHit(1));
      SoLIDGEMHit *hit3 = (SoLIDGEMHit*)(thistrack->GetHit(2));
      SoLIDGEMHit *hit4 = (SoLIDGEMHit*)(thistrack->GetHit(3));
      SoLIDGEMHit *hit5 = (SoLIDGEMHit*)(thistrack->GetHit(4));
      Int_t flag = 1;

      if (flag == 1)
	{
	  nmom  = 0;
	  CalCircle(hit1->GetX(),hit1->GetY(),hit2->GetX(),hit2->GetY(),hit3->GetX(),hit3->GetY(),&Rmom[nmom],&Xc[nmom],&Yc[nmom]);
	  Rc[nmom] = sqrt(Xc[nmom]*Xc[nmom]+Yc[nmom]*Yc[nmom]);
	  nmom++;
	  CalCircle(hit1->GetX(),hit1->GetY(),hit2->GetX(),hit2->GetY(),hit4->GetX(),hit4->GetY(),&Rmom[nmom],&Xc[nmom],&Yc[nmom]);
	  Rc[nmom] = sqrt(Xc[nmom]*Xc[nmom]+Yc[nmom]*Yc[nmom]);
	  nmom++;
	  CalCircle(hit1->GetX(),hit1->GetY(),hit3->GetX(),hit3->GetY(),hit4->GetX(),hit4->GetY(),&Rmom[nmom],&Xc[nmom],&Yc[nmom]);
	  Rc[nmom] = sqrt(Xc[nmom]*Xc[nmom]+Yc[nmom]*Yc[nmom]);
	  nmom++;
	  CalCircle(hit2->GetX(),hit2->GetY(),hit3->GetX(),hit3->GetY(),hit4->GetX(),hit4->GetY(),&Rmom[nmom],&Xc[nmom],&Yc[nmom]);
	  Rc[nmom] = sqrt(Xc[nmom]*Xc[nmom]+Yc[nmom]*Yc[nmom]);
	  nmom++;
	  CalCircle(hit1->GetX(),hit1->GetY(),hit2->GetX(),hit2->GetY(),hit5->GetX(),hit5->GetY(),&Rmom[nmom],&Xc[nmom],&Yc[nmom]);
	  Rc[nmom] = sqrt(Xc[nmom]*Xc[nmom]+Yc[nmom]*Yc[nmom]);
	  nmom++;
	  CalCircle(hit1->GetX(),hit1->GetY(),hit3->GetX(),hit3->GetY(),hit5->GetX(),hit5->GetY(),&Rmom[nmom],&Xc[nmom],&Yc[nmom]);
	  Rc[nmom] = sqrt(Xc[nmom]*Xc[nmom]+Yc[nmom]*Yc[nmom]);
	  nmom++;
	  CalCircle(hit2->GetX(),hit2->GetY(),hit3->GetX(),hit3->GetY(),hit5->GetX(),hit5->GetY(),&Rmom[nmom],&Xc[nmom],&Yc[nmom]);
	  Rc[nmom] = sqrt(Xc[nmom]*Xc[nmom]+Yc[nmom]*Yc[nmom]);
	  nmom++;
	  CalCircle(hit1->GetX(),hit1->GetY(),hit4->GetX(),hit4->GetY(),hit5->GetX(),hit5->GetY(),&Rmom[nmom],&Xc[nmom],&Yc[nmom]);
	  Rc[nmom] = sqrt(Xc[nmom]*Xc[nmom]+Yc[nmom]*Yc[nmom]);
	  nmom++;
	  CalCircle(hit2->GetX(),hit2->GetY(),hit4->GetX(),hit4->GetY(),hit5->GetX(),hit5->GetY(),&Rmom[nmom],&Xc[nmom],&Yc[nmom]);
	  Rc[nmom] = sqrt(Xc[nmom]*Xc[nmom]+Yc[nmom]*Yc[nmom]);
	  nmom++;
	  CalCircle(hit3->GetX(),hit3->GetY(),hit4->GetX(),hit4->GetY(),hit5->GetX(),hit5->GetY(),&Rmom[nmom],&Xc[nmom],&Yc[nmom]);
	  Rc[nmom] = sqrt(Xc[nmom]*Xc[nmom]+Yc[nmom]*Yc[nmom]);
	  nmom++;
	   
	  Double_t vertexzAvg = 0.;
	  Double_t thetaAvg = 0.;
	  Int_t countAvg = 0;
	  //123
	  if ( Rmom[0]-(-0.0374431+1.0118*Rc[0])< 0.06&& Rmom[0]-(-0.0374431+1.0118*Rc[0])>-0.1
	     && CalVertexZ(hit1, hit2, hit3, &Rmom[0], kFAEC, &reconVertexZ, &reconTheta) ){
	      vertexzAvg += reconVertexZ;
	      thetaAvg += reconTheta;
	      countAvg++;
	  }else{
	    flag = 0;
	  }
	  //124
	  if ( Rmom[1]-(-0.0345228+1.01045*Rc[1])< 0.06&& Rmom[1]-(-0.0345228+1.01045*Rc[1])>-0.1
	     && CalVertexZ(hit1, hit2, hit4, &Rmom[1], kFAEC, &reconVertexZ, &reconTheta) ){
	      vertexzAvg += reconVertexZ;
	      thetaAvg += reconTheta;
	      countAvg++;
	  }else{
	    flag = 0;
	  }
	  //134
	  if ( Rmom[2]-(-0.032753+1.0098*Rc[2])< 0.06&& Rmom[2]-(-0.032753+1.0098*Rc[2])>-0.1
	     && CalVertexZ(hit1, hit3, hit4, &Rmom[2], kFAEC, &reconVertexZ, &reconTheta) ){
	      vertexzAvg += reconVertexZ;
	      thetaAvg += reconTheta;
	      countAvg++;
	  }else{
	    flag = 0;
	  }
	  //234
	  if ( Rmom[3]-(-0.0316562+1.00954*Rc[3])<0.06 && Rmom[3]-(-0.0316562+1.00954*Rc[3])>-0.1
	     && CalVertexZ(hit2, hit3, hit4, &Rmom[3], kFAEC, &reconVertexZ, &reconTheta) ){
	      vertexzAvg += reconVertexZ;
	      thetaAvg += reconTheta;
	      countAvg++;
	  }else{
	    flag = 0;
	  }
	  //125
	  if ( Rmom[4]-(-0.0315397+1.00942*Rc[4])< 0.06&& Rmom[4]-(-0.0315397+1.00942*Rc[4])>-0.1
	     && CalVertexZ(hit1, hit2, hit5, &Rmom[4], kFAEC, &reconVertexZ, &reconTheta) ){
	      vertexzAvg += reconVertexZ;
	      thetaAvg += reconTheta;
	      countAvg++;
	  }else{
	    flag = 0;
	  }
	  //135
	  if ( Rmom[5]-(-0.028547+1.00847*Rc[5])<0.06 && Rmom[5]-(-0.028547+1.00847*Rc[5])>-0.1
	     && CalVertexZ(hit1, hit3, hit5, &Rmom[5], kFAEC, &reconVertexZ, &reconTheta) ){
	      vertexzAvg += reconVertexZ;
	      thetaAvg += reconTheta;
	      countAvg++;
	  }else{
	    flag = 0;
	  }
	  //235
	  if ( Rmom[6]-(-0.0259752+1.00773*Rc[6])< 0.06&& Rmom[6]-(-0.0259752+1.00773*Rc[6])>-0.1
	     && CalVertexZ(hit2, hit3, hit5, &Rmom[6], kFAEC, &reconVertexZ, &reconTheta) ){
	      vertexzAvg += reconVertexZ;
	      thetaAvg += reconTheta;
	      countAvg++;
	  }else{
	    flag = 0;
	  }
	  //145
	  if ( Rmom[7]-(-0.0227897+1.00688*Rc[7])<0.06 && Rmom[7]-(-0.0227897+1.00688*Rc[7])>-0.15
	     && CalVertexZ(hit1, hit4, hit5, &Rmom[7], kFAEC, &reconVertexZ, &reconTheta) ){
	      vertexzAvg += reconVertexZ;
	      thetaAvg += reconTheta;
	      countAvg++;
	  }else{
	    flag = 0;
	  }
	  //245
	  if ( Rmom[8]-(-0.0190085+1.00589*Rc[8])<0.06 && Rmom[8]-(-0.0190085+1.00589*Rc[8])>-0.15
	     && CalVertexZ(hit2, hit4, hit5, &Rmom[8], kFAEC, &reconVertexZ, &reconTheta) ){
	      vertexzAvg += reconVertexZ;
	      thetaAvg += reconTheta;
	      countAvg++;
	  }else{
	    flag = 0;
	  }
	  //345
	  if ( Rmom[9]-(-0.0116765+1.00414*Rc[9])<0.07 && Rmom[9]-(-0.0116765+1.00414*Rc[9])>-0.2
	     && CalVertexZ(hit3, hit4, hit5, &Rmom[9], kFAEC, &reconVertexZ, &reconTheta) ){
	      vertexzAvg += reconVertexZ;
	      thetaAvg += reconTheta;
	      countAvg++;
	  }else{
	    flag = 0;
	  }
	  //cross-check
	  if (sqrt((Xc[8]-Xc[7])*(Xc[8]-Xc[7])+(Yc[8]-Yc[7])*(Yc[8]-Yc[7]))<0.16*sqrt(2) && Rmom[8]-Rmom[7]-(0.00134498+1.04092*(Rc[8]-Rc[7]))<0.05 && Rmom[8]-Rmom[7]-(0.00134498+1.04092*(Rc[8]-Rc[7]))>-0.05){
	  }else{
	    flag = 0;
	  }
	   
	  if (sqrt((Xc[6]-Xc[5])*(Xc[6]-Xc[5])+(Yc[6]-Yc[5])*(Yc[6]-Yc[5]))<0.3*sqrt(2)&& Rmom[6]-Rmom[5]-(-0.000678468+1.0465*(Rc[6]-Rc[5]))<0.05 && Rmom[6]-Rmom[5]-(-0.000678468+1.0465*(Rc[6]-Rc[5]))>-0.05){
	  }else{
	    flag = 0;
	  }
	  if (flag == 1) {
	    assert(countAvg == 10);
	    reconTheta = thetaAvg/10.;
	    reconVertexZ = vertexzAvg/10.;
	  }
	}
      if (flag==1){
	           thistrack->SetStatus(kTRUE);
	           thistrack->SetCoarseChi2((sqrt((Xc[1]-Xc[0])*(Xc[1]-Xc[0])+(Yc[1]-Yc[0])*(Yc[1]-Yc[0]))
			       +sqrt((Xc[2]-Xc[1])*(Xc[2]-Xc[1])+(Yc[2]-Yc[1])*(Yc[2]-Yc[1]))
			       +sqrt((Xc[3]-Xc[2])*(Xc[3]-Xc[2])+(Yc[3]-Yc[2])*(Yc[3]-Yc[2]))
			       +sqrt((Xc[4]-Xc[3])*(Xc[4]-Xc[3])+(Yc[4]-Yc[3])*(Yc[4]-Yc[3]))
			       +sqrt((Xc[5]-Xc[4])*(Xc[5]-Xc[4])+(Yc[5]-Yc[4])*(Yc[5]-Yc[4]))
			       +sqrt((Xc[6]-Xc[5])*(Xc[6]-Xc[5])+(Yc[6]-Yc[5])*(Yc[6]-Yc[5]))
			       +sqrt((Xc[7]-Xc[6])*(Xc[7]-Xc[6])+(Yc[7]-Yc[6])*(Yc[7]-Yc[6]))
			       +sqrt((Xc[8]-Xc[7])*(Xc[8]-Xc[7])+(Yc[8]-Yc[7])*(Yc[8]-Yc[7]))
			       +sqrt((Xc[9]-Xc[8])*(Xc[9]-Xc[8])+(Yc[9]-Yc[8])*(Yc[9]-Yc[8]))
			       +sqrt((Xc[0]-Xc[9])*(Xc[0]-Xc[9])+(Yc[0]-Yc[9])*(Yc[0]-Yc[9]))
			       )/10.);
			       thistrack->SetCoarseFitStatus(kTRUE);
			       thistrack->SetVertexZ(reconVertexZ);
	           thistrack->SetTheta(reconTheta);
      }else{
	     thistrack->SetCoarseFitStatus(kTRUE);
	      thistrack->SetStatus(kFALSE);
	      thistrack->SetCoarseChi2(kINFINITY);
      }

    }
  }
  fCoarseTracks->Sort();
}
//_________________________________________________________________________________
Int_t ProgressiveTracking::FindMomRange(Int_t layer1, Double_t phi0, Int_t layer2, Double_t phi1, 
					Double_t* mom_min,Double_t* mom_max,Int_t angleflag){
  
  Double_t phi2 = phi1-phi0;
  if (phi2<0)
    {
      phi2 = phi2 + 2.*3.1415926;
    }
  if (layer1==0&&layer2==1&&angleflag==kLAEC){
    Double_t par[6]={0.0502148,-0.00143936,-0.003,0.006, 0.3, -0.1};
    Double_t tempmin,tempmax;
    tempmin = TMath::Abs(par[0]/(phi2-par[1]-par[2])-par[4]);
    tempmax = TMath::Abs(par[0]/(phi2-par[1]-par[3])-par[5]);
    if (tempmax > 11.) tempmax = 11.;
    if (tempmin>*mom_min) *mom_min = tempmin;
    if (tempmax<*mom_max) *mom_max = tempmax;
  }else if(layer1==0&&layer2==2&&angleflag==kLAEC){
    Double_t par[6]={0.118713,-0.00310567,-0.005,0.01, 0.3, -0.1};
    Double_t tempmin,tempmax;
    tempmin = TMath::Abs(par[0]/(phi2-par[1]-par[2])-par[4]);
    tempmax = TMath::Abs(par[0]/(phi2-par[1]-par[3])-par[5]);
    if (tempmax > 11.) tempmax = 11.;
    if (tempmin>*mom_min) *mom_min = tempmin;
    if (tempmax<*mom_max) *mom_max = tempmax;
  }else if(layer1==1&&layer2==2&&angleflag==kLAEC){
    Double_t par[6]={0.06932,-0.00184587,-0.003,0.008, 0.3, -0.1};
    Double_t tempmin,tempmax;
    tempmin = TMath::Abs(par[0]/(phi2-par[1]-par[2])-par[4]);
    tempmax = TMath::Abs(par[0]/(phi2-par[1]-par[3])-par[5]);
    if (tempmax > 11.) tempmax = 11.;
    if (tempmin>*mom_min) *mom_min = tempmin;
    if (tempmax<*mom_max) *mom_max = tempmax;
  }else if(layer1==1&&layer2==3&&angleflag==kLAEC){
    Double_t par[6]={0.196285,-0.00521495,-0.006,0.015, 0.3, -0.1};
    Double_t tempmin,tempmax;
    tempmin = TMath::Abs(par[0]/(phi2-par[1]-par[2])-par[4]);
    tempmax = TMath::Abs(par[0]/(phi2-par[1]-par[3])-par[5]);
    if (tempmax > 11.) tempmax = 11.;
    if (tempmin>*mom_min) *mom_min = tempmin;
    if (tempmax<*mom_max) *mom_max = tempmax;
  }else if(layer1==2&&layer2==3&&angleflag==kLAEC){
    Double_t par[6]={0.127511,-0.0034897,-0.005,0.01, 0.3, -0.1};
    Double_t tempmin,tempmax;
    tempmin = TMath::Abs(par[0]/(phi2-par[1]-par[2])-par[4]);
    tempmax = TMath::Abs(par[0]/(phi2-par[1]-par[3])-par[5]);
    if (tempmax > 11.) tempmax = 11.;
    if (tempmin>*mom_min) *mom_min = tempmin;
    if (tempmax<*mom_max) *mom_max = tempmax;
  }else if(layer1==1&&layer2==2&&angleflag==kFAEC){
    Double_t par[6]={0.0714656,-0.00281194,-0.003,0.004, 0.5, -0.1};
    Double_t tempmin,tempmax;
    tempmin = TMath::Abs(par[0]/(phi2-par[1]-par[2])-par[4]);
    tempmax = TMath::Abs(par[0]/(phi2-par[1]-par[3])-par[5]);
    if (tempmax > 11.) tempmax = 11.;
    if (tempmin>*mom_min) *mom_min = tempmin;
    if (tempmax<*mom_max) *mom_max = tempmax;
  }else if(layer1==1&&layer2==3&&angleflag==kFAEC){
    Double_t par[6]={0.198465,-0.00711192,-0.005,0.006, 0.6, -0.1};
    Double_t tempmin,tempmax;
    tempmin = TMath::Abs(par[0]/(phi2-par[1]-par[2])-par[4]);
    tempmax = TMath::Abs(par[0]/(phi2-par[1]-par[3])-par[5]);
    if (tempmax > 11.) tempmax = 11.;
    if (tempmin>*mom_min) *mom_min = tempmin;
    if (tempmax<*mom_max) *mom_max = tempmax;
  }else if(layer1==2&&layer2==3&&angleflag==kFAEC){
    Double_t par[6]={0.12868,-0.00465945,-0.004,0.005, 0.5, -0.1};
    Double_t tempmin,tempmax;
    tempmin = TMath::Abs(par[0]/(phi2-par[1]-par[2])-par[4]);
    tempmax = TMath::Abs(par[0]/(phi2-par[1]-par[3])-par[5]);
    if (tempmax > 11.) tempmax = 11.;
    if (tempmin>*mom_min) *mom_min = tempmin;
    if (tempmax<*mom_max) *mom_max = tempmax;
  }else if(layer1==2&&layer2==4&&angleflag==kFAEC){
    Double_t par[6]={0.326676,-0.0117276,-0.01,0.015, 0.4, -0.2};
    Double_t tempmin,tempmax;
    tempmin = TMath::Abs(par[0]/(phi2-par[1]-par[2])-par[4]);
    tempmax = TMath::Abs(par[0]/(phi2-par[1]-par[3])-par[5]);
    if (tempmax > 11.) tempmax = 11.;
    if (tempmin>*mom_min) *mom_min = tempmin;
    if (tempmax<*mom_max) *mom_max = tempmax;
  }else if(layer1==3&&layer2==4&&angleflag==kFAEC){
    Double_t par[6]={0.198643,-0.00721605,-0.006,0.008, 0.5, -0.1};
    Double_t tempmin,tempmax;
    tempmin = TMath::Abs(par[0]/(phi2-par[1]-par[2])-par[4]);
    tempmax = TMath::Abs(par[0]/(phi2-par[1]-par[3])-par[5]);
    if (tempmax > 11.) tempmax = 11.;
    if (tempmin>*mom_min) *mom_min = tempmin;
    if (tempmax<*mom_max) *mom_max = tempmax;
  }else if(layer1==3&&layer2==5&&angleflag==kFAEC){
    Double_t par[6]={0.442922,-0.0168963,-0.01,0.015, 0.5, -0.1};
    Double_t tempmin,tempmax;
    tempmin = TMath::Abs(par[0]/(phi2-par[1]-par[2])-par[4]);
    tempmax = TMath::Abs(par[0]/(phi2-par[1]-par[3])-par[5]);
    if (tempmax > 11.) tempmax = 11.;
    if (tempmin>*mom_min) *mom_min = tempmin;
    if (tempmax<*mom_max) *mom_max = tempmax;
  }else if(layer1==4&&layer2==5&&angleflag==kFAEC){
    Double_t par[6]={0.24772,-0.0103976,-0.006,0.009, 0.5, -0.1};
    Double_t tempmin,tempmax;
    tempmin = TMath::Abs(par[0]/(phi2-par[1]-par[2])-par[4]);
    tempmax = TMath::Abs(par[0]/(phi2-par[1]-par[3])-par[5]);
    if (tempmax > 11.) tempmax = 11.;
    if (tempmin>*mom_min) *mom_min = tempmin;
    if (tempmax<*mom_max) *mom_max = tempmax;
  }
  if (*mom_min<=*mom_max){
    return 1;
  }else{
    //cout<<layer1<<" "<<layer2<<" "<<*mom_min<<" "<<*mom_max<<" "<<phi2<<endl;
    return 0;
  }

}
//___________________________________________________________________________________
void ProgressiveTracking::CombineTrackRoad(TClonesArray* theTracks)
{
  Int_t nGoodTrack = 0;
  for (Int_t i=0;i!=fCoarseTracks->GetLast()+1;i++){
    SoLIDTrack *thistrack = (SoLIDTrack*)(fCoarseTracks->At(i));
    if (thistrack->GetStatus()){
      Int_t flag = 0;
      //(this part may not be perfect since the flag will be 1 as long as one 
      //of the hits are the same)
      for (Int_t j=0;j!=thistrack->GetNHits();j++){
	     SoLIDGEMHit *thishit = (SoLIDGEMHit*)(thistrack->GetHit(j));
	     Int_t layer = thishit->GetTrackerID();
	for (Int_t k=0;k!=fNGoodHits[layer];k++){
	  if ((thishit->GetX() == ((SoLIDGEMHit*)(fGoodHits[layer]->At(k)))->GetX()) && (thishit->GetY() == ((SoLIDGEMHit*)(fGoodHits[layer]->At(k)))->GetY())){
	    flag = 1;
	  }
       	}
      }
      if (flag == 0){
      SoLIDTrack* newtrack = 0;
      if (fDoMC){
#ifdef MCDATA
				    newtrack = new ((*theTracks)[nGoodTrack++]) SoLIDMCTrack();
#endif
				  }
				  else{
				    newtrack = new ((*theTracks)[nGoodTrack++]) SoLIDTrack();  
				  }
	newtrack->SetStatus(kTRUE);
	newtrack->SetCoarseFitStatus(kTRUE);
	newtrack->SetAngleFlag(thistrack->GetAngleFlag());
	newtrack->SetCharge(thistrack->GetCharge());
	newtrack->SetCoarseChi2(thistrack->GetChi2());
	newtrack->SetMomMin(thistrack->GetMomMin());
	newtrack->SetMomMax(thistrack->GetMomMax());
	newtrack->SetThetaMin(thistrack->GetThetaMin());
	newtrack->SetThetaMax(thistrack->GetThetaMax());
	for (Int_t j=0;j!=thistrack->GetNHits();j++){
	  SoLIDGEMHit* thishit = 0;
    Int_t layer = 0;
	  if (fDoMC){
#ifdef MCDATA
      thishit = (SoLIDMCGEMHit*)(thistrack->GetHit(j));
	    layer = thishit->GetTrackerID();
	    new ((*fGoodHits[layer])[fNGoodHits[layer]++]) SoLIDMCGEMHit(*(dynamic_cast<SoLIDMCGEMHit*>(thishit)));
#endif
	  }
	  else{
	    thishit = (SoLIDGEMHit*)(thistrack->GetHit(j));
	    layer = thishit->GetTrackerID();
	    new ((*fGoodHits[layer])[fNGoodHits[layer]++]) SoLIDGEMHit(*thishit);
	  }
	  assert(thishit != 0);
	  newtrack->AddHit(thishit); 
	}
      }
    }
  }
}
//___________________________________________________________________________________
Int_t ProgressiveTracking::FindThetaRange(Double_t r1, Int_t layer1, Double_t r2, Int_t layer2,
                                          Double_t* theta_min,Double_t* theta_max, Int_t angleflag)
{

  Double_t tempmin,tempmax;
  if (layer1==0&&layer2==1&&angleflag==kLAEC){
    Double_t par[4]={1.69204, 201.147, 1.2, -1.2};
    tempmax = par[0]+par[1]*(r2-r1)+par[2];
    tempmin = par[0]+par[1]*(r2-r1)+par[3];
    if (tempmin>*theta_min) *theta_min = tempmin;
    if (tempmax<*theta_max) *theta_max = tempmax;
  }else if (layer1==0 && layer2==2&&angleflag==kLAEC){
    Double_t par[4]={1.60416, 90.3763, 1.5, -1.};
    tempmax = par[0]+par[1]*(r2-r1)+par[2];
    tempmin = par[0]+par[1]*(r2-r1)+par[3];
    if (tempmin>*theta_min) *theta_min = tempmin;
    if (tempmax<*theta_max) *theta_max = tempmax;
  }else if (layer1==1 && layer2==2&&angleflag==kLAEC){
    Double_t par[4]={1.67377, 162.847, 1.5, -1.2};
    tempmax = par[0]+par[1]*(r2-r1)+par[2];
    tempmin = par[0]+par[1]*(r2-r1)+par[3];
    if (tempmin>*theta_min) *theta_min = tempmin;
    if (tempmax<*theta_max) *theta_max = tempmax;
  }else if (layer1==1 && layer2==3&&angleflag==kLAEC){
    Double_t par[4]={1.6713, 61.8236, 2., -1.};
    tempmax = par[0]+par[1]*(r2-r1)+par[2];
    tempmin = par[0]+par[1]*(r2-r1)+par[3];
    if (tempmin>*theta_min) *theta_min = tempmin;
    if (tempmax<*theta_max) *theta_max = tempmax;
  }else if (layer1==2 && layer2==3&&angleflag==kLAEC){
    Double_t par[4]={1.76328, 99.1469, 2, -1};
    tempmax = par[0]+par[1]*(r2-r1)+par[2];
    tempmin = par[0]+par[1]*(r2-r1)+par[3]; 
    if (tempmin>*theta_min) *theta_min = tempmin;
    if (tempmax<*theta_max) *theta_max = tempmax;
  }else if (layer1==1 && layer2==2&&angleflag==kFAEC){
    Double_t par[4]={0.35128, 177.195, 1.5, -1.2};
    tempmax = par[0]+par[1]*(r2-r1)+par[2];
    tempmin = par[0]+par[1]*(r2-r1)+par[3]; 
    if (tempmin>*theta_min) *theta_min = tempmin;
    if (tempmax<*theta_max) *theta_max = tempmax;
  }else if (layer1==1 && layer2==3&&angleflag==kFAEC){
    Double_t par[4]={0.279777, 67.6577, 1.2, -1.};
    tempmax = par[0]+par[1]*(r2-r1)+par[2];
    tempmin = par[0]+par[1]*(r2-r1)+par[3]; 
    if (tempmin>*theta_min) *theta_min = tempmin;
    if (tempmax<*theta_max) *theta_max = tempmax;
  }else if (layer1==2 && layer2==3&&angleflag==kFAEC){
    Double_t par[4]={0.31703, 108.612, 1.5, -1.};
    tempmax = par[0]+par[1]*(r2-r1)+par[2];
    tempmin = par[0]+par[1]*(r2-r1)+par[3]; 
    if (tempmin>*theta_min) *theta_min = tempmin;
    if (tempmax<*theta_max) *theta_max = tempmax;
  }else if (layer1==2 && layer2==4&&angleflag==kFAEC){
    Double_t par[4]={0.360605, 44.7933, 2., -1.};
    tempmax = par[0]+par[1]*(r2-r1)+par[2];
    tempmin = par[0]+par[1]*(r2-r1)+par[3];
    if (tempmin>*theta_min) *theta_min = tempmin;
    if (tempmax<*theta_max) *theta_max = tempmax;
  }else if (layer1==3 && layer2==4&&angleflag==kFAEC){
    Double_t par[4]={0.448258, 75.8155, 2.2, -1.};
    tempmax = par[0]+par[1]*(r2-r1)+par[2];
    tempmin = par[0]+par[1]*(r2-r1)+par[3];
    if (tempmin>*theta_min) *theta_min = tempmin;
    if (tempmax<*theta_max) *theta_max = tempmax;
  }else if (layer1==3 && layer2==5&&angleflag==kFAEC){
    Double_t par[4]={0.717657, 34.1022, 3., -1.};
    tempmax = par[0]+par[1]*(r2-r1)+par[2];
    tempmin = par[0]+par[1]*(r2-r1)+par[3];
    if (tempmin>*theta_min) *theta_min = tempmin;
    if (tempmax<*theta_max) *theta_max = tempmax;
  }else if (layer1==4 && layer2==5&&angleflag==kFAEC){
    Double_t par[4]={1.02099, 61.4545, 4, -1};
    tempmax = par[0]+par[1]*(r2-r1)+par[2];
    tempmin = par[0]+par[1]*(r2-r1)+par[3];
    if (tempmin>*theta_min) *theta_min = tempmin;
    if (tempmax<*theta_max) *theta_max = tempmax;
  }
  
  if (*theta_min<=*theta_max){
    return 1;
  }else{
    return 0;
  }
}
//________________________________________________________________________________
inline void ProgressiveTracking::CalCircle(Double_t x1,Double_t y1,Double_t x2,Double_t y2,Double_t x3,
                                           Double_t y3, Double_t* R,Double_t* Xc, Double_t* Yc){
  if (y1!=y2&& (x1-x2)*(y2-y3)!=(x2-x3)*(y1-y2)){
    *Xc = ((y2-y3)*(x1*x1-x2*x2+y1*y1-y2*y2)/2.-(y1-y2)*(x2*x2-x3*x3+y2*y2-y3*y3)/2.)/((x1-x2)*(y2-y3)-(x2-x3)*(y1-y2));
    *Yc = ((x1*x1-x2*x2+y1*y1-y2*y2)/2.-(*Xc)*(x1-x2))/(y1-y2);
    *R = sqrt((*Xc-x1)*(*Xc-x1)+(*Yc-y1)*(*Yc-y1));
    
  }else if (y2!=y3&&(x2-x3)*(y3-y1)!=(x3-x1)*(y2-y3)){
    *Xc = ((y3-y1)*(x2*x2-x3*x3+y2*y2-y3*y3)/2.-(y2-y3)*(x3*x3-x1*x1+y3*y3-y1*y1)/2.)/((x2-x3)*(y3-y1)-(x3-x1)*(y2-y3));
    *Yc = ((x2*x2-x3*x3+y2*y2-y3*y3)/2.-(*Xc)*(x2-x3))/(y2-y3);
    *R = sqrt((*Xc-x2)*(*Xc-x2)+(*Yc-y2)*(*Yc-y2));
   
  }else if (y3!=y1&&(x3-x1)*(y1-y2)!=(x1-x2)*(y3-y1)){
    *Xc = ((y1-y2)*(x3*x3-x1*x1+y3*y3-y1*y1)/2.-(y3-y1)*(x1*x1-x2*x2+y1*y1-y2*y2)/2.)/((x3-x1)*(y1-y2)-(x1-x2)*(y3-y1));
    *Yc = ((x3*x3-x1*x1+y3*y3-y1*y1)/2.-(*Xc)*(x3-x1))/(y3-y1);
    *R = sqrt((*Xc-x3)*(*Xc-x3)+(*Yc-y3)*(*Yc-y3));
    
  }else{
    *Xc = 0;
    *Yc = 0;
    *R = 0;
  }
  
}
//________________________________________________________________________________________
inline Double_t ProgressiveTracking::CalDeltaPhi(Double_t phi1, Double_t phi2)
{
  Double_t deltaPhi = phi1 - phi2;
  deltaPhi = TVector2::Phi_mpi_pi(deltaPhi);
  return deltaPhi;
}
//________________________________________________________________________________________
inline Double_t ProgressiveTracking::CalChargeAsy(Double_t qu, Double_t qv)
{
  assert(qu+qv != 0.);
  return (qu-qv)/(qu+qv);
}
//_____________________________________________________________________________________________
Bool_t ProgressiveTracking::CalVertexZ(SoLIDGEMHit* hit1, SoLIDGEMHit* hit2, SoLIDGEMHit *hit3, Double_t *Rmom, Int_t angleFlag,
                                Double_t *reconVertexZ, Double_t *reconTheta)
{
  assert(hit1->GetTrackerID() < hit2->GetTrackerID() && hit2->GetTrackerID() < hit3->GetTrackerID());
  Double_t rotphi[3] = {0.};
  Double_t x1 = hit1->GetX();
  Double_t x2 = hit2->GetX();
  Double_t x3 = hit3->GetX();
  Double_t y1 = hit1->GetY();
  Double_t y2 = hit2->GetY();
  Double_t y3 = hit3->GetY();
  Double_t z1 = hit1->GetZ();
  Double_t z2 = hit2->GetZ();
  Double_t z3 = hit3->GetZ();
  
  rotphi[0] = 2.*TMath::ASin(TMath::Sqrt(pow(x2 - x1, 2)+pow(y2 - y1, 2))/2./(*Rmom));
  rotphi[1] = 2.*TMath::ASin(TMath::Sqrt(pow(x3 - x2, 2)+pow(y3 - y2, 2))/2./(*Rmom));
  rotphi[2] = 2.*TMath::ASin(TMath::Sqrt(pow(x3 - x1, 2)+pow(y3 - y1, 2))/2./(*Rmom));
  *reconTheta = (TMath::ATan((*Rmom)*rotphi[0]/(z2 - z1))+ 
                 TMath::ATan((*Rmom)*rotphi[1]/(z3 - z2))+ TMath::ATan((*Rmom)*rotphi[2]/(z3 - z1)))/3.;
  *reconTheta = (*reconTheta)*180./TMath::Pi();
  Bool_t pass = kTRUE;
  if (angleFlag == kLAEC){
    if (hit1->GetTrackerID() == 0 && hit2->GetTrackerID() == 1 && hit3->GetTrackerID() == 2){
    
       *reconTheta += 0.0031035 + -0.0817294*rotphi[0] + -46.1366*rotphi[0]*rotphi[0];
       *reconVertexZ = -TMath::Sqrt(x1*x1+y1*y1)/TMath::Tan((*reconTheta)/180.*TMath::Pi());
       *reconVertexZ = -1.9853 + 0.742107*(*reconVertexZ) + -0.0711807*(*reconVertexZ)*(*reconVertexZ);
       if ( (*reconVertexZ) > -3.2 || (*reconVertexZ) < -3.8) { pass = kFALSE; }  
    }
    else if (hit1->GetTrackerID() == 0 && hit2->GetTrackerID() == 1 && hit3->GetTrackerID() == 3){
       *reconTheta += -0.00499077 + 0.452268*rotphi[0] + -61.6921*rotphi[0]*rotphi[0];
       *reconVertexZ = -TMath::Sqrt(x1*x1+y1*y1)/TMath::Tan((*reconTheta)/180.*TMath::Pi());
       *reconVertexZ = -1.77836 + 0.975878*(*reconVertexZ) + -0.00546159*(*reconVertexZ)*(*reconVertexZ);
       if ( (*reconVertexZ) > -3.2 || (*reconVertexZ) < -3.8) { pass = kFALSE; }
    }
    else if (hit1->GetTrackerID() == 0 && hit2->GetTrackerID() == 2 && hit3->GetTrackerID() == 3){
      *reconTheta += -0.0152338 + 0.357038*rotphi[0] + -13.6711*rotphi[0]*rotphi[0];
      *reconVertexZ = -TMath::Sqrt(x1*x1+y1*y1)/TMath::Tan((*reconTheta)/180.*TMath::Pi());
      *reconVertexZ = -1.75096 + 1.00334*(*reconVertexZ) + 0.00117971*(*reconVertexZ)*(*reconVertexZ);
      if ( (*reconVertexZ) > -3.2 || (*reconVertexZ) < -3.8) { pass = kFALSE; }  
    }
    else if (hit1->GetTrackerID() == 1 && hit2->GetTrackerID() == 2 && hit3->GetTrackerID() == 3){
      *reconTheta += -0.00676339 + 0.356862*rotphi[0] + -45.2325*rotphi[0]*rotphi[0];
      *reconVertexZ = -TMath::Sqrt(x1*x1+y1*y1)/TMath::Tan((*reconTheta)/180.*TMath::Pi());
      *reconVertexZ = -1.5141 + 0.995579*(*reconVertexZ) + 0.000449439*(*reconVertexZ)*(*reconVertexZ);
      if ( (*reconVertexZ) > -3.2 || (*reconVertexZ) < -3.8) { pass = kFALSE; }
    }
  }
  else if (angleFlag == kFAEC){
    if (hit1->GetTrackerID() == 1 && hit2->GetTrackerID() == 2 && hit3->GetTrackerID() == 3){
      *reconTheta += 0.0108978 + -0.297584*rotphi[0] + -24.4199*rotphi[0]*rotphi[0];
      *reconVertexZ = -TMath::Sqrt(x1*x1+y1*y1)/TMath::Tan((*reconTheta)/180.*TMath::Pi());
      *reconVertexZ = -0.104674 + 2.46256*(*reconVertexZ) + 0.380788*(*reconVertexZ)*(*reconVertexZ);
      if ( (*reconVertexZ) > -3.2 || (*reconVertexZ) < -3.8) { pass = kFALSE; }
    }
    else if (hit1->GetTrackerID() == 1 && hit2->GetTrackerID() == 2 && hit3->GetTrackerID() == 4){
      *reconTheta += 0.00611748 + -0.0187155*rotphi[0] + -32.6257*rotphi[0]*rotphi[0];
      *reconVertexZ = -TMath::Sqrt(x1*x1+y1*y1)/TMath::Tan((*reconTheta)/180.*TMath::Pi());
      *reconVertexZ = -0.31037 + 2.24668*(*reconVertexZ) + 0.324404*(*reconVertexZ)*(*reconVertexZ);
      if ( (*reconVertexZ) > -3.2 || (*reconVertexZ) < -3.8) { pass = kFALSE; }      
    }
    else if (hit1->GetTrackerID() == 1 && hit2->GetTrackerID() == 3 && hit3->GetTrackerID() == 4){
      *reconTheta += -0.00901722 + 0.163468*rotphi[0] + -5.49729*rotphi[0]*rotphi[0];
      *reconVertexZ = -TMath::Sqrt(x1*x1+y1*y1)/TMath::Tan((*reconTheta)/180.*TMath::Pi());
      *reconVertexZ = -0.539831 + 2.0037*(*reconVertexZ) + 0.260459*(*reconVertexZ)*(*reconVertexZ);
      if ( (*reconVertexZ) > -3.2 || (*reconVertexZ) < -3.8) { pass = kFALSE; }
    }
    else if (hit1->GetTrackerID() == 2 && hit2->GetTrackerID() == 3 && hit3->GetTrackerID() == 4){
      *reconTheta += -0.000859183 + 0.0709105*rotphi[0] + -14.2965*rotphi[0]*rotphi[0];
      *reconVertexZ = -TMath::Sqrt(x1*x1+y1*y1)/TMath::Tan((*reconTheta)/180.*TMath::Pi());
      *reconVertexZ = 0.526385 + 2.5496*(*reconVertexZ) + 0.347781*(*reconVertexZ)*(*reconVertexZ);
      if ( (*reconVertexZ) > -3.2 || (*reconVertexZ) < -3.8) { pass = kFALSE; }
    }
    else if (hit1->GetTrackerID() == 1 && hit2->GetTrackerID() == 2 && hit3->GetTrackerID() == 5){
      *reconTheta += -0.000143187 + 0.230999*rotphi[0] + -38.6706*rotphi[0]*rotphi[0];
      *reconVertexZ = -TMath::Sqrt(x1*x1+y1*y1)/TMath::Tan((*reconTheta)/180.*TMath::Pi());
      *reconVertexZ = -0.382876 + 2.17201*(*reconVertexZ) + 0.305258*(*reconVertexZ)*(*reconVertexZ);
      if ( (*reconVertexZ) > -3.2 || (*reconVertexZ) < -3.8) { pass = kFALSE; }
    }
    else if (hit1->GetTrackerID() == 1 && hit2->GetTrackerID() == 3 && hit3->GetTrackerID() == 5){
      *reconTheta += -0.0105459 + 0.202513*rotphi[0] + -6.17457*rotphi[0]*rotphi[0];
      *reconVertexZ = -TMath::Sqrt(x1*x1+y1*y1)/TMath::Tan((*reconTheta)/180.*TMath::Pi());
      *reconVertexZ = -0.649638 + 1.88758*(*reconVertexZ) + 0.229914*(*reconVertexZ)*(*reconVertexZ);
      if ( (*reconVertexZ) > -3.2 || (*reconVertexZ) < -3.8) { pass = kFALSE; }
    }
    else if (hit1->GetTrackerID() == 2 && hit2->GetTrackerID() == 3 && hit3->GetTrackerID() == 5){
      *reconTheta += -0.00514327 + 0.191252*rotphi[0] + -16.2823*rotphi[0]*rotphi[0];
      *reconVertexZ = -TMath::Sqrt(x1*x1+y1*y1)/TMath::Tan((*reconTheta)/180.*TMath::Pi());
      *reconVertexZ = 0.437893 + 2.46523*(*reconVertexZ) + 0.327882*(*reconVertexZ)*(*reconVertexZ);
      if ( (*reconVertexZ) > -3.2 || (*reconVertexZ) < -3.8) { pass = kFALSE; }
    }
    else if (hit1->GetTrackerID() == 2 && hit2->GetTrackerID() == 4 && hit3->GetTrackerID() == 5){
      *reconTheta += -0.00775158 + 0.0747926*rotphi[0] + -2.87822*rotphi[0]*rotphi[0];
      *reconVertexZ = -TMath::Sqrt(x1*x1+y1*y1)/TMath::Tan((*reconTheta)/180.*TMath::Pi());
      *reconVertexZ = 0.371146 + 2.40374*(*reconVertexZ) + 0.313795*(*reconVertexZ)*(*reconVertexZ); 
      if ( (*reconVertexZ) > -3.2 || (*reconVertexZ) < -3.8) { pass = kFALSE; }
    }
    else if (hit1->GetTrackerID() == 2 && hit2->GetTrackerID() == 4 && hit3->GetTrackerID() == 5){
      *reconTheta += -5.48472e-05 + -0.0508633*rotphi[0] + -8.36052*rotphi[0]*rotphi[0];
      *reconVertexZ = -TMath::Sqrt(x1*x1+y1*y1)/TMath::Tan((*reconTheta)/180.*TMath::Pi());
      *reconVertexZ = 2.10465 + 3.06856*(*reconVertexZ) + 0.381855*(*reconVertexZ)*(*reconVertexZ);
      if ( (*reconVertexZ) > -3.1 || (*reconVertexZ) < -3.9) { pass = kFALSE; }
    }
    else if (hit1->GetTrackerID() == 1 && hit2->GetTrackerID() == 4 && hit3->GetTrackerID() == 5){
      *reconTheta += -0.0114525 + 0.0978447*rotphi[0] + -1.79422*rotphi[0]*rotphi[0];
      *reconVertexZ = -TMath::Sqrt(x1*x1+y1*y1)/TMath::Tan((*reconTheta)/180.*TMath::Pi());
      *reconVertexZ = -0.611219 + 1.92689*(*reconVertexZ) + 0.239943*(*reconVertexZ)*(*reconVertexZ);
      if ( (*reconVertexZ) > -3.2 || (*reconVertexZ) < -3.8) { pass = kFALSE; }
    }
  }
  else{
    *reconTheta = 0;
    *reconVertexZ = 0;
    return kFALSE;
  }
  
  if (!pass){
    *reconTheta = 0;
    *reconVertexZ = 0;
    return kFALSE;
  }else{
    return kTRUE;
  }

}










































