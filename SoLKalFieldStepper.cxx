//c++
#include <iostream>
#include <cassert>
#include <cmath>
//ROOT
#include "TMath.h"
//SoLIDTracking
#include "SoLKalFieldStepper.h"
#include "SoLKalTrackSite.h"
#include "SoLIDGEMHit.h"
using namespace std;
SoLKalFieldStepper * SoLKalFieldStepper::fSoLKalFieldStepper = NULL;

//_________________________________________________________________
SoLKalFieldStepper::SoLKalFieldStepper()
{
  fIsMSOn = kTRUE;
  fIsDEDXOn = kTRUE;
  UseDefaultStep();
  fMass = kElectronMass;
  fCharge = -1.;
  fIsElectron = kTRUE;
  fFieldMap = SoLIDFieldMap::GetInstance();
  InitDetMaterial();
}
//_________________________________________________________________
SoLKalFieldStepper::~SoLKalFieldStepper()
{
}
//_________________________________________________________________
void SoLKalFieldStepper::UseDefaultStep()
{
  initialStepSize = 0.2;     //m
  stepSizeDec     = 0.75;
  stepSizeInc     = 1.25;
  maxStepSize     = 0.5;     //m
  minStepSize     = 0.05;    //m
  minPrecision    = 1.e-4;
  maxPrecision    = 1.e-5;
  maxNumSteps     = 1e3;
  maxDist         = 5.e-3;   //m
  minLengthCalcQ  = 1e-2;    //m
}
//_________________________________________________________________
void SoLKalFieldStepper::UseFineStep()
{
  initialStepSize = 0.02;     //m
  stepSizeDec     = 0.75;
  stepSizeInc     = 1.25;
  maxStepSize     = 0.05;     //m
  minStepSize     = 0.01;    //m
  minPrecision    = 1.e-4;
  maxPrecision    = 1.e-5;
  maxNumSteps     = 1e3;
  maxDist         = 1.e-3;   //m
  minLengthCalcQ  = 1e-3;    //m
}
//__________________________________________________________________
//a simple fixed step Runge-Kutta propagation method
void SoLKalFieldStepper::PropagationClassicalRK4(TVector3 &inMom, TVector3 &inPos, 
                         Double_t &finalZ, Double_t &charge, 
                         Double_t &stepSize, TVector3 &fiMom, TVector3 &fiPos)
{
  const Int_t nvar = 6;
  Double_t yt[nvar], yIn[nvar], yOut[nvar];
  Double_t dydx[nvar], dydxt[nvar], dydxm[nvar];
  Double_t h = stepSize; //initial step size
  Double_t hh = h*0.5;
  Double_t h6 = h/6.0;
  Double_t delta_z;
  Double_t delta_z_save;
  Double_t distance2plane = 0.1;//(m) stop propagation when the particle is this close to the target plane
  Bool_t do_loop = kTRUE;
  Double_t direction;
  Int_t maxStep = 1000;

  if (inPos.Z() >= finalZ) direction = -1.;
  else direction = 1.;
  
  delta_z = direction*(finalZ - inPos.Z());
  
  if (delta_z < h){
      h = delta_z; //make sure the initial step size is less than the distance between two planes
      hh = h*0.5;
      h6 = h/6.0;
  }
  
  h*=direction;
  hh*=direction;
  h6*=direction;

  yIn[0] = inPos.X(); yIn[1] = inPos.Y(); yIn[2] = inPos.Z(); 
  yIn[3] = inMom.X(); yIn[4] = inMom.Y(); yIn[5] = inMom.Z();

  Int_t countStep = 0;  
  while(do_loop && countStep < maxStep){
      countStep++;
      //the classcial 4th order Runge-Kutta method happends here
      Int_t i;
      RightHandSide(yIn, charge, inMom.Mag(), dydx);
      for(i=0;i<nvar;i++){
          yt[i] = yIn[i] + hh*dydx[i] ;             // 1st Step K1=h*dydx
      }

      RightHandSide(yt,charge,inMom.Mag(), dydxt);

      for(i=0;i<nvar;i++){
          yt[i] = yIn[i] + hh*dydxt[i] ;
      }

      RightHandSide(yt,charge,inMom.Mag(), dydxm);

      for(i=0;i<nvar;i++){
           yt[i]   = yIn[i] + h*dydxm[i] ;
           dydxm[i] += dydxt[i] ;                    // now dydxm=(K2+K3)/h
      }

      RightHandSide(yt,charge,inMom.Mag(), dydxt) ;
      // Final RK4 output
      for(i=0;i<nvar;i++)    {
         yOut[i] = yIn[i]+h6*(dydx[i]+dydxt[i]+2.0*dydxm[i]); //+K1/6+K4/6+(K2+K3)/3
      }
      //decide if more propagation is needed and if the step size needs to be changed
      delta_z_save = delta_z;
      
      delta_z = direction*(finalZ- yOut[2]);
      if (delta_z > direction*h){
           //the particle is still far away from the target plane
           //continue the loop without chaning delta_s
           do_loop = kTRUE;
           for (i=0; i<nvar; i++){
               yIn[i] = yOut[i];
           }
      }
       else if (delta_z < direction*h && delta_z > distance2plane){
           //the particle is approaching the target plane
           //change the step size to delta_z so that the propagation will never pass the target plane
           do_loop = kTRUE;
           h = direction*delta_z;
           hh = h*0.5;
           h6 = h/6.0;
           for (i=0; i<nvar; i++){
               yIn[i] = yOut[i];
           }
      }
      else if (delta_z < distance2plane && delta_z > 0.){
           //stop the propagation
           do_loop = kFALSE;
      }else{
           cout<<"SoLKalFieldStepper::PropagationClassicalRK4, propagation pass the target plane"<<endl;
           do_loop = kFALSE;
      }

    }
    //prepare for final output
    delta_z = finalZ - yOut[2];
    
    fiPos.SetXYZ(yOut[0] + delta_z*(yOut[3]/yOut[5]), 
                 yOut[1] + delta_z*(yOut[4]/yOut[5]), 
                 yOut[2] + delta_z);
    fiMom.SetXYZ(yOut[3], yOut[4], yOut[5]);
}
//________________________________________________________________________________________________
void SoLKalFieldStepper::RightHandSide(const Double_t y[], const Double_t charge, 
                                       const Double_t mom_mag, Double_t dydx[])
{
  Double_t momentum_mag_square = y[3]*y[3] + y[4]*y[4] + y[5]*y[5];
  Double_t inv_momentum_magnitude = 1.0 / std::sqrt( momentum_mag_square );
  Double_t cof = (charge*TMath::C()/1.e10)/sqrt(momentum_mag_square);
  TVector3 B = fFieldMap->GetBField(y[0], y[1], y[2]);
  
  dydx[0] = y[3]*inv_momentum_magnitude;       //  (d/ds)x = Vx/V
  dydx[1] = y[4]*inv_momentum_magnitude;       //  (d/ds)y = Vy/V
  dydx[2] = y[5]*inv_momentum_magnitude;       //  (d/ds)z = Vz/V

  dydx[3] = cof*(y[4]*B.Z() - y[5]*B.Y()) ;   // Ax = a*(Vy*Bz - Vz*By)
  dydx[4] = cof*(y[5]*B.X() - y[3]*B.Z()) ;   // Ay = a*(Vz*Bx - Vx*Bz)
  dydx[5] = cof*(y[3]*B.Y() - y[4]*B.X()) ;   // Az = a*(Vx*By - Vy*Bx)
}
//__________________________________________________________________________________________________
Double_t SoLKalFieldStepper::Distance2Points(const TVector3 &vec1, const TVector3 &vec2)
{
  // Calculates the distance between two points.
    TVector3 distV = vec1 - vec2;
    return distV.Mag();
}
//__________________________________________________________________________________________________
void SoLKalFieldStepper::InitTrack(Double_t &mass, Double_t &charge, Bool_t &isElectron, Bool_t dir)
{
  fMass = mass;
  fCharge = charge;
  fIsElectron = isElectron;
  fIsBackward = dir;
}
//__________________________________________________________________________________________________
void SoLKalFieldStepper::Transport(const SoLKalTrackSite  &from, // site from
                                   SoLKalTrackSite   &to,   // sit to
                                   SoLKalMatrix      &sv,   // state vector
                                   SoLKalMatrix      &F,    // propagator matrix
                                   SoLKalMatrix      &Q)   // process noise matrix
{
  assert(from.GetCurState().GetType() == SoLKalTrackSite::kFiltered);
  Double_t z = to.GetZ();
  Transport(from.GetCurState(), z, sv, F, Q);
}
//__________________________________________________________________________________________________
void SoLKalFieldStepper::Transport(const SoLKalTrackSite  &from, // site from
                                   Double_t         &finalZ, // z position of the destination
                                   SoLKalMatrix       &sv,   // state vector
                                   SoLKalMatrix       &F,    // propagator matrix
                                   SoLKalMatrix       &Q)   // process noise matrix
{
  assert(from.GetCurState().GetType() == SoLKalTrackSite::kFiltered);
  Transport(from.GetCurState(), finalZ, sv, F, Q);
}
//__________________________________________________________________________________________________
void SoLKalFieldStepper::Transport(const SoLKalTrackState  &sv_from, // site from
                                   Double_t         &finalZ, // z position of the destination
                                   SoLKalMatrix       &sv,   // state vector
                                   SoLKalMatrix       &F,    // propagator matrix
                                   SoLKalMatrix       &Q)   // process noise matrix
{
  trackPosAtZ    = sv_from.GetZ0();
  
  trackLength    = 0.;
  stepLength     = 0.;
  jstep          = 0;
  energyLoss     = 0.;
  Double_t beta;
  Bool_t bCalcJac = kTRUE;
  
  if (trackPosAtZ >= finalZ) fIsBackward = kTRUE;
  else fIsBackward = kFALSE;
  
  SoLKalMatrix sv_to(kSdim, 1);
  SoLKalMatrix sv_PreStep(kSdim, 1);
  SoLKalMatrix DF(kSdim, kSdim);                  // propagator matrix segment
  for (Int_t i=0; i<kSdim; i++){
      sv_to(i,0) = sv_from(i,0);
  }
  F.UnitMatrix(); // initialize F to unity
  Q.Zero();       // initialize Q to zero
  TVector3 posPreStep;
  TVector3 posAt; 
  
  posAt.SetXYZ(sv_to(0,0), sv_to(1,0), trackPosAtZ);
  
  Double_t step = initialStepSize; 
  Double_t nextStep = step;
  Double_t stepFac = 1.;
  Bool_t doNextStep = kTRUE;
  // Reduce step size if too close to the target plane.
  TVector3 dirAt;
  SoLKalTrackState::CalcDir(dirAt, sv_to);
  TVector3 pointIntersect;
  if (!FindTargetPlaneIntersection(pointIntersect, finalZ, dirAt, posAt)){
	   cout<<"track is parallel to the GEM tracker"<<endl;
	}

	Double_t d = Distance2Points(posAt, pointIntersect);
  if(d < step) {
	 step     = d;
	 nextStep = step;
  }
  
  // Steps are taken in z-direction. Compensate for track inclination.
  Double_t stepz = step * dirAt.Z();
  //-----------------------initialization end--------------------------//
  
  //-----------------------begin Runge-Kutta stepping------------------//
  
   while (d >= maxDist && doNextStep && jstep < maxNumSteps && step){
	   // Store some properties before stepping.
	   
	   posPreStep.SetXYZ(sv_to(kIdxX0, 0), sv_to(kIdxY0, 0), trackPosAtZ);
	   sv_PreStep = sv_to;
	   // Calculate step size in z-direction.
	   TVector3 dirAt;
	   SoLKalTrackState::CalcDir(dirAt, sv_to);
	   stepz = step * dirAt.z();
	   
	   
	   stepFac = RKPropagation(sv_to, DF, stepz, bCalcJac, fIsBackward); // do one step
	   
	   posAt.SetXYZ(sv_to(kIdxX0, 0), sv_to(kIdxY0, 0), trackPosAtZ);//update position vector after RK propagation
       
	   Double_t sd =  trackPosAtZ - finalZ;
	   Double_t prec = 1.e-4;
	   if(((sd > 0. + prec) && fIsBackward == kFALSE) || ((sd < 0. - prec) && fIsBackward == kTRUE)) {
	     // Track went past target plane during propagtion. Repeating last step with reduced step size.
	     step       *= stepSizeDec;
	     nextStep   *= stepSizeDec;
	     sv_to  = sv_PreStep;
	     trackPosAtZ = posPreStep.Z();
	     posAt.SetXYZ(posPreStep.X(), posPreStep.Y(), posPreStep.Z());
	     trackLength -= stepLength;
	     continue;
	   }
       
	   // After each Runge-Kutta step i the covariance matrix C_i would have to be
	   // transported with the propagator matrix F_i:
	   // C_1 = F_1 * C_0   * F_1^T
	   // C_2 = F_2 * C_1   * F_2^T
	   // ...
	   // C_n = F_n * C_n-1 * F_n^T
	   // This can be transformed to:
	   // C_n = F_n * ... * F_2 * F_1 * C_0 * F_1^T * F_2^T * ... * F_n^T
	   //     =             F         * C_0 *         F^T 
       
	   F = DF * F; // update propagator
       

	   //calculating process noice and energy loss 
	   //SoLKalMatrix DFt  = SoLKalMatrix(TMatrixD::kTransposed, DF);
	   //SoLKalMatrix Qms(kSdim, kSdim);
	   //Qms.Zero();
	   
	   if ((IsMSOn() || IsDEDXOn()) && stepLength > minLengthCalcQ )
	     {
	       // Track inclination = sqrt(1 + tx^2 + ty^2).
	       Double_t trackIncl = TMath::Sqrt(1. + sv_PreStep(kIdxTX, 0)*sv_PreStep(kIdxTX, 0) 
						+ sv_PreStep(kIdxTY, 0)*sv_PreStep(kIdxTY, 0));
	       // Total step length = step length in z direction * track inclination.
	       
	       stepLength = stepz * trackIncl;
	       
	       beta = 1. / TMath::Sqrt(1. + (fMass*fMass*1e-6) * sv_PreStep(kIdxQP, 0)*sv_PreStep(kIdxQP, 0));
	       
	       if ( IsMSOn() )     /*Double_t cms2 = */CalcMultScat(Q, sv_PreStep, stepLength, beta);
	       if ( IsDEDXOn() )   energyLoss += CalcEnergyLoss(Q, sv_to, stepLength, sv_PreStep(kIdxQP, 0), beta);
	     }
       //end calculating process noice and energy loss
       
       //Q = DF * (Q + Qms) * DFt;
       

	   // Decide if more steps must be done and calculate new step size.
	   // -----------
	   SoLKalTrackState::CalcDir(dirAt, sv_to);
	   if(!FindTargetPlaneIntersection(pointIntersect,finalZ, dirAt, posAt)) {
	     cout<<"No intersection with target plane found."<<endl;
	   }
	   d = Distance2Points(posAt, pointIntersect);
       
	   if (step < nextStep) {  // last step was rest step to chamber
	     if (stepFac < 1.) nextStep = step * stepFac;
	   } else {
	     nextStep *= stepFac;
	   }
       
	   if (d > nextStep || d < maxDist) step = nextStep;
	   else step = d;
     
	   if(fIsBackward == kFALSE) {
	     if(posAt.z() < pointIntersect.z()) { doNextStep = kTRUE; }
	     else { doNextStep = kFALSE; }
	   } else {
	     if(posAt.z() > pointIntersect.z()) { doNextStep = kTRUE; }
	     else { doNextStep = kFALSE; }
	   }
     
	 }
       //-----------------------end Runge-Kutta stepping--------------------//
  
   // To make sure the track position is on the target layer propagate to the target plane
   // using a straight line.
   sv_PreStep = sv_to;
   posPreStep.SetXYZ(sv_PreStep(kIdxX0, 0), sv_PreStep(kIdxY0, 0), trackPosAtZ);
   
   if(!PropagateStraightLine(sv_to, DF, trackPosAtZ, finalZ, fIsBackward)) {
	   cout<<"TKalDetCradle::Transport: final propagation to target plane failed"<<endl;
   }
   
       //calculate the energy loss and multiple scattering for this straight line propagation
   
   if ((IsMSOn() || IsDEDXOn()) && stepLength > minLengthCalcQ ){ 
	   
	   beta = 1. / TMath::Sqrt(1. + (fMass*fMass*1e-6) * sv_PreStep(kIdxQP, 0)*sv_PreStep(kIdxQP, 0));
	   
	   if ( IsMSOn() )     /*Double_t cms2 = */CalcMultScat(Q, sv_PreStep, stepLength, beta);
	   if ( IsDEDXOn() )   energyLoss += CalcEnergyLoss(Q, sv_to, stepLength, sv_PreStep(kIdxQP, 0), beta);
	   
	 }
   //make a correction for the GEM detector for now
   //SoLKalMatrix DFt  = SoLKalMatrix(TMatrixD::kTransposed, DF);
	 //SoLKalMatrix Qms(kSdim, kSdim);
	 //Qms.Zero();
   
   beta = 1. / TMath::Sqrt(1. + (fMass*fMass*1e-6) * sv_PreStep(kIdxQP, 0)*sv_PreStep(kIdxQP, 0));
   
   Double_t theta = atan( sqrt(pow(sv_to(kIdxTX, 0), 2) + pow(sv_to(kIdxTY, 0), 2)) );
	 Double_t distInGEM = 1.5525e-2 / cos(theta);
	 
	 if ( IsMSOn() )     CalcMultScat(Q, sv_PreStep, distInGEM, beta, kGEM);
	 if ( IsDEDXOn() )   energyLoss += CalcEnergyLoss(Q, sv_to, distInGEM, sv_PreStep(kIdxQP, 0), beta, kGEM);
	 
	 //Q = DF * (Q + Qms) * DFt;
   
   F = DF * F; // final propagator matrix
   sv = sv_to;
   
}
//___________________________________________________________________________________________________
Double_t SoLKalFieldStepper::RKPropagation(SoLKalMatrix &stateVec, SoLKalMatrix &fPropStep, Double_t stepSize, 
                                           Bool_t bCalcJac, Bool_t dir)
{
    // One step of track tracing from track state.
    //
    // Input and output:
    // stateVec:  Starting track parameters.
    // fPropStep: Change in propagator matrix for this step.
    //
    // Input:
    // totStep:   Step size.
    // bCalcJac:  Update the Jacobian (propagator matrix).
    
    const Int_t numPars = 4; // x, y, tx, ty
    const Double_t kappa = TMath::C() / 1.e10; // in GeV/(c * kG * m)
    
    //  Changes of state vector a by stepping in z direction:
    //  da/dz, a = (x, y, tx, ty, qp), z = stepping dir
    //  x'  = tx
    //  y'  = ty
    // tx'  = k * (tx * ty   * B[0] - (1 + tx2) * B[1] + ty * B[2])
    // ty'  = k * ((1 + ty2) * B[0] - txty      * B[1] - tx * B[2])
    //  k   = c * (q/p) * sqrt(1 + tx2 + ty2)

    // f    = f(z, a(z))
    // k1   = f(zn        , an)
    // k2   = f(zn + 1/2*h, an + 1/2*k1)
    // k3   = f(zn + 1/2*h, an + 1/2*k2)
    // k4   = f(zn + h    , an + k3)
    // fn+1 = fn + 1/6*k1 + 1/3*k2 + 1/3*k3 + 1/6*k4 + O(h^5)

    // Constants for RK stepping.
    static Double_t a[numPars] = { 0.0    , 0.5    , 0.5    , 1.0     };
    static Double_t b[numPars] = { 0.0    , 0.5    , 0.5    , 1.0     };
    static Double_t c[numPars] = { 1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0 };
    // for rk stepping
    //Int_t step4;
    const Int_t rksteps = 4; // The 4 points used in Runge-Kutta stepping: start, 2 mid points and end
    const Int_t rkStart = 0;
    const Int_t rkMid1  = 1;
    const Int_t rkMid2  = 2;
    const Int_t rkEnd   = 3;

    // Change of track parameters in each step.
    // First index: step point (0 = start; 1,2 = mid, 3 = end), second index: state parameter.
    Double_t k[rksteps][numPars];

    // for propagator matrix
    Double_t F_tx[numPars], F_ty[numPars], F_tx_tx[numPars], F_ty_tx[numPars], F_tx_ty[numPars], F_ty_ty[numPars];
    Double_t F2_tx[numPars], F2_ty[numPars]; // @2tx/@z^2 and @2ty/@z^2

    //----------------------------------------------------------------
    Double_t est     = 0.; // error estimation
    Double_t stepFac = 1.;

    TVector3 B;                  // B-field
    Double_t h = stepSize;       // step size
    if(fIsBackward == kTRUE) {
        h *= -1;                 // stepping in negative z-direction
    }
    Double_t half    = h * 0.5;  // half step interval for fourth order RK
    
    Double_t qp_in   = stateVec(kIdxQP, 0);
    TVector3 posFrom = TVector3(stateVec(kIdxX0, 0), stateVec(kIdxY0, 0), trackPosAtZ);
    Double_t z_in    = posFrom.z();
    trackPosAtZ      = z_in;
    TVector3 posAt   = posFrom;

    Double_t hC      = h * kappa;

    // Input track state vector and state vector during stepping.
    Double_t sv_in[numPars], sv_step[numPars];
    sv_in[kIdxX0]       = stateVec(kIdxX0, 0);
    sv_in[kIdxY0]       = stateVec(kIdxY0, 0);
    sv_in[kIdxTX]   = stateVec(kIdxTX, 0);
    sv_in[kIdxTY] = stateVec(kIdxTY, 0);

    fPropStep.UnitMatrix();

    //------------------------------------------------------------------------
    //   Runge-Kutta step
    //
    Int_t istep;
    Int_t ipar;
    
    do {
        half  = h * 0.5;
        for (istep = 0; istep < rksteps; ++istep) {  // k1,k2,k3,k4  (k1=start, k2,k3 = half, k4=end of interval)
            for(ipar = 0; ipar < numPars; ++ipar) {  // 4 track parameters
                if(istep == 0) {
                    sv_step[ipar] = sv_in[ipar];     // in first step copy input track parameters (x,y,tx,ty)
                } else {
                    sv_step[ipar] = sv_in[ipar] + b[istep] * k[istep-1][ipar];     // do step
                }
            }
            trackPosAtZ = z_in + a[istep] * h; // move z along with track
            posAt.SetXYZ(sv_step[kIdxX0], sv_step[kIdxY0], trackPosAtZ ); // update z value for current position

            
            //get the magnatic field 
	          B = fFieldMap->GetBField(posAt.X(), posAt.Y(), posAt.Z());

            Double_t tx        = sv_step[kIdxTX];
            Double_t ty        = sv_step[kIdxTY];
            Double_t tx2       = tx * tx;
            Double_t ty2       = ty * ty;
            Double_t txty      = tx * ty;
            Double_t tx2ty21   = 1.0 + tx2 + ty2;

            Double_t I_tx2ty21 = 1.0 / tx2ty21 * qp_in;
            Double_t tx2ty2    = sqrt(tx2ty21);
                     tx2ty2   *= hC;
            Double_t tx2ty2qp  = tx2ty2 * qp_in;

            // for state propagation
            F_tx[istep] = ( txty        *B.x() - ( 1.0 + tx2 )*B.y() + ty*B.z()) * tx2ty2;    // h * @tx/@z / (qp) = h * tx' / (qp)
            F_ty[istep] = (( 1.0 + ty2 )*B.x() - txty         *B.y() - tx*B.z()) * tx2ty2;    // h * @ty/@z / (qp) = h * ty' / (qp)

            //------------------------------------------------------------------------
            // for transport matrix
            F_tx_tx[istep] = F_tx[istep]*tx*I_tx2ty21 + ( ty*B.x()-2.0*tx*B.y() ) * tx2ty2qp; // h * @tx'/@tx
            F_tx_ty[istep] = F_tx[istep]*ty*I_tx2ty21 + ( tx*B.x()+B.z()        ) * tx2ty2qp; // h * @tx'/@ty
            F_ty_tx[istep] = F_ty[istep]*tx*I_tx2ty21 + (-ty*B.y()-B.z()        ) * tx2ty2qp; // h * @ty'/@tx
            F_ty_ty[istep] = F_ty[istep]*ty*I_tx2ty21 + ( 2.0*ty*B.x()-tx*B.y() ) * tx2ty2qp; // h * @ty'/@ty

            // Change of track parameters in each step.
            k[istep][kIdxX0]       = tx * h;              // dx
            k[istep][kIdxY0]       = ty * h;              // dy
            k[istep][kIdxTX]   = F_tx[istep] * qp_in; // dtx
            k[istep][kIdxTY] = F_ty[istep] * qp_in; // dty

            // h * @tx'/@z = h * (@tx/@dz)/@dz
            F2_tx[istep]   = qp_in *
                (tx*k[istep][kIdxTX]/TMath::Sqrt(tx2ty21)*(txty*B.X() - (1. + tx2)*B.Y() + ty*B.Z())
                  + TMath::Sqrt(tx2ty21) *
                 (k[istep][kIdxTX]*ty*B.X() + tx*k[istep][kIdxTY]*B.X()
                  - 2.*tx*k[istep][kIdxTX]*B.Y() + k[istep][kIdxTY]*B.Z())
                );

            // h * @ty'/@z = h * (@ty/@dz)/@dz
            F2_ty[istep]   = qp_in *
                (ty*k[istep][kIdxTY]/TMath::Sqrt(tx2ty21)*(( 1.0 + ty2 )*B.X() - txty*B.Y() - tx*B.Z())
                  + TMath::Sqrt(tx2ty21) *
                 (2.*ty*k[istep][kIdxTY]*B.X() - k[istep][kIdxTX]*ty*B.Y() - tx*k[istep][kIdxTY]*B.Y()
                 - k[istep][kIdxTX]*B.Z()));

        }  // end of Runge-Kutta steps
        //------------------------------------------------------------------------

        //------------------------------------------------------------------------
        // error estimation ala Geant
        est = 0.;

        est += fabs(k[rkStart][kIdxX0]       + k[rkEnd][kIdxX0]       - k[rkMid1][kIdxX0]       - k[rkMid2][kIdxX0])       * half;
        est += fabs(k[rkStart][kIdxY0]       + k[rkEnd][kIdxY0]       - k[rkMid1][kIdxY0]       - k[rkMid2][kIdxY0])       * half;
        est += fabs(k[rkStart][kIdxTX]       + k[rkEnd][kIdxTX]       - k[rkMid1][kIdxTX]       - k[rkMid2][kIdxTX])       * half;
        est += fabs(k[rkStart][kIdxTY]       + k[rkEnd][kIdxTY]       - k[rkMid1][kIdxTY]       - k[rkMid2][kIdxTY])       * half;
        //------------------------------------------------------------------------
        if (fabs(est) < minPrecision || fabs(h) <= minStepSize || stepFac <= minStepSize) {
            // we found a step size with good precision
            jstep ++;
            break;
        } else {
            // precision not good enough. make smaller step
            stepFac *= stepSizeDec;
            h       *= stepSizeDec;
            hC       = h * kappa;

#if rkDebug > 2
            cout<<"Precision not good enough. Reducing step size to "<<h<<endl;
#endif
        }

    } while (jstep < maxNumSteps);
    
    if (est < maxPrecision && fabs(h) < maxStepSize) {
        stepFac *= stepSizeInc;
    }
    
    //------------------------------------------------------------------------
    // set output track parameters
    for(ipar = 0; ipar < numPars; ++ ipar) {
        // yn+1        = yn         + 1/6*k1                    + 1/3*k2                  + 1/3*k3                  + 1/6*k4
      stateVec(ipar, 0) = sv_in[ipar]+c[rkStart]*k[rkStart][ipar]+c[rkMid1]*k[rkMid1][ipar]+c[rkMid2]*k[rkMid2][ipar]+c[rkEnd]*k[rkEnd][ipar];
    }


    //------------------------------------------------------------------------
    
    stepLength = fabs(Distance2Points(posFrom, posAt));
    trackLength += stepLength;  // calculate track length

    if (est < maxPrecision && fabs(h) < maxStepSize) {
        stepFac *= stepSizeInc;
    }

    if(!bCalcJac) {
        return stepFac;
    }
    
    //------------------------------------------------------------------------
    //
    //     Derivatives

    //     x  y tx ty qp
    //  x  0  5 10 15 20    J[ 0 ...] = da/dx       = 0,   1 = 1
    //  y  1  6 11 16 21    J[ 5 ...] = da/dy       = 0,   6 = 1
    // tx  2  7 12 17 22    J[10 ...] = da/dtx            12 = 1, 14 = 0
    // ty  3  8 13 18 23    J[15 ...] = da/dty            18 = 1, 19 = 0
    // qp  4  9 14 19 24    J[20 ...] = da/dqp            24 = 1
    //


    //------------------------------------------------------------------------
    //
    //     Derivatives    dx/dqp
    //

    // Update of Jacobian matrix for each step:
    // F_i = @a(zi)/@a(z0) =  @f(a,zi)/@a(zi) * (I + F_i-1 * (zi - z0)/dz) * dz
    //
    // z0      : starting z position
    // zi      : z position after step i = (z0, z0 + 0.5*dz, z0 + 0.5*dz, z0 + dz)
    // dz      : step size
    // (zi - z0) / dz = (0, 0.5, 0.5, 1)
    //
    // F_i = @a(zi)/@a(z0) : Jacobian matrix of RK step i
    // f(a,zi) = @a(zi)/@zi: derivation of track state vector a at step i by z
    // E                   : unit matrix
    // a(zi)               : track state vector at step i
    // a(z0)               : initial track state vector

    // x0: contribution of the Unit matrix to Jacobian
    // x:  contribution of Jacobian from previous step
    Double_t x0[numPars], x[numPars];
    
    // Elements of the Jacobian matrix F_i of one step i for the derivation of one specific track parameter.
    Double_t k1[rksteps][numPars];

    x0[kIdxX0] = 0.0; x0[kIdxY0] = 0.0; x0[kIdxTX] = 0.0; x0[kIdxTY] = 0.0;

    //   Runge-Kutta step for derivatives dx/dqp

    for(istep = 0; istep < rksteps; ++ istep) {
        for(ipar = 0; ipar < numPars; ++ ipar) {
            if(istep == 0) {
                x[ipar] = x0[ipar];
            } else {
                x[ipar] = x0[ipar] + b[istep] * k1[istep-1][ipar];
            }
        }
        //     F_i(x, qp)       = @x'/@tx * F_(i-1)(tx,qp) * dz
        //                      = F_(i-1)(tx,qp) * dz
        k1[istep][kIdxX0]       = x[kIdxTX]   * h;
        //     F_i(y, qp)       = @y'/@ty * F_(i-1)(ty,qp) * dz
        //                      = F_(i-1)(ty,qp) * dz
        k1[istep][kIdxY0]       = x[kIdxTY] * h;
        //     F_i(tx, qp)      = @tx'/@qp    + @tx'/@tx * F_(i-1)(tx,qp) * dz + @tx'/@ty * F_(i-1)(ty,qp) * dz
        k1[istep][kIdxTX]   = F_tx[istep] + F_tx_tx[istep] * x[kIdxTX] + F_tx_ty[istep] * x[kIdxTY];
        //     F_i(ty, qp)      = @ty'/@qp    + @ty'/@tx * F_(i-1)(tx,qp) * dz + @ty'/@ty * F_(i-1)(ty,qp) * dz
        k1[istep][kIdxTY] = F_ty[istep] + F_ty_tx[istep] * x[kIdxTX] + F_ty_ty[istep] * x[kIdxTY];
    }  // end of Runge-Kutta steps for derivatives dx/dqp


    for (ipar = 0; ipar < numPars; ++ipar ) {
        fPropStep(ipar, kIdxQP) = x0[ipar] + c[rkStart]*k1[rkStart][ipar] + c[rkMid1]*k1[rkMid1][ipar]
                                           + c[rkMid2] *k1[rkMid2][ipar]  + c[rkEnd] *k1[rkEnd][ipar];
    }
    fPropStep(kIdxQP,kIdxQP) = 1.;
    //------------------------------------------------------------------------
    
    //------------------------------------------------------------------------
    //     Derivatives    dx/tx
    //

    x0[kIdxX0] = 0.0; x0[kIdxY0] = 0.0; x0[kIdxTX] = 1.0; x0[kIdxTY] = 0.0;

    //
    //   Runge-Kutta step for derivatives dx/dtx
    //

    for (istep = 0; istep < 4; ++ istep) {
        for(ipar = 0; ipar < numPars; ++ipar) {
            if(istep == 0) {
                x[ipar] = x0[ipar];
            } else if ( ipar != kIdxTX ){
                x[ipar] = x0[ipar] + b[istep] * k1[istep-1][ipar];
            }
        }

        k1[istep][kIdxX0]       = x[kIdxTX]   * h;
        k1[istep][kIdxY0]       = x[kIdxTY] * h;
        //k1[istep][kIdxTX] = F_tx_tx[istep] * x[kIdxTX] + F_tx_ty[istep] * x[kIdxTY]; // not needed
        k1[istep][kIdxTY] = F_ty_tx[istep] * x[kIdxTX] + F_ty_ty[istep] * x[kIdxTY];
    }  // end of Runge-Kutta steps for derivatives dx/dtx

    for(ipar = 0; ipar < numPars; ++ipar ) {
        if(ipar != kIdxTX) {
            fPropStep(ipar, kIdxTX) = x0[ipar] + c[rkStart]*k1[rkStart][ipar] + c[rkMid1]*k1[rkMid1][ipar]
                                                   + c[rkMid2] *k1[rkMid2][ipar]  + c[rkEnd] *k1[rkEnd][ipar];
        }
    }
    //      end of derivatives dx/dtx
    fPropStep(kIdxTX, kIdxTY) = 1.;
    fPropStep(kIdxQP,     kIdxTX)   = 0.;
    //------------------------------------------------------------------------
    
    //------------------------------------------------------------------------
    //     Derivatives    dx/ty
    //

    x0[kIdxX0] = 0.0; x0[kIdxY0] = 0.0; x0[kIdxTX] = 0.0; x0[kIdxTY] = 1.0;

    //
    //   Runge-Kutta step for derivatives dx/dty
    //

    for (istep = 0; istep < 4; ++ istep) {
        for(ipar = 0; ipar < numPars; ++ipar) {
            if(istep == 0) {
                x[ipar] = x0[ipar];           // ty fixed
            } else if(ipar != kIdxTY) {
                x[ipar] = x0[ipar] + b[istep] * k1[istep-1][ipar];
                //x[ipar] = x0[ipar] + b[istep] * k1[istep*4-4+ipar];
            }
        }

        k1[istep][kIdxX0]         = x[kIdxTX]   * h;
        k1[istep][kIdxY0]         = x[kIdxTY] * h;
        k1[istep][kIdxTX]     = F_tx_tx[istep] * x[kIdxTX] + F_tx_ty[istep] * x[kIdxTY];
        //k1[istep][kIdxTY] = F_ty_tx[istep] * x[kIdxTX] + F_ty_ty[istep] * x[kIdxTY]; // not needed
    }  // end of Runge-Kutta steps for derivatives dx/dty

    for(ipar = 0; ipar < 3; ++ipar ) {
        fPropStep(ipar, kIdxTY) = x0[ipar] + c[rkStart]*k1[rkStart][ipar] + c[rkMid1]*k1[rkMid1][ipar]
                                                 + c[rkMid2] *k1[rkMid2][ipar]  + c[rkEnd] *k1[rkEnd][ipar];
    }
    //      end of derivatives dx/dty
    fPropStep(kIdxTY,kIdxTY) = 1.;
    fPropStep(kIdxQP      ,kIdxTY) = 0.;
    //------------------------------------------------------------------------
    
    //------------------------------------------------------------------------
    //
    //    derivatives dx/dx and dx/dy

    for(ipar = 0; ipar < numPars + 1; ipar++) {
        fPropStep(ipar, kIdxX0) = 0.;
        fPropStep(ipar, kIdxY0) = 0.;
    }
    fPropStep(kIdxX0, kIdxX0) = 1.;
    fPropStep(kIdxY0, kIdxY0) = 1.;


    // Propagator has entry for z-coordinate as well.

    if(fPropStep.GetNrows() == 6) {
        Double_t x0z[numPars+1], xz[numPars+1];

        Double_t k1z[rksteps][numPars+1];

        Int_t idxZ0 = 4;
        x0z[kIdxX0] = 0.; x0z[kIdxY0] = 0.; x0z[kIdxTX] = 0.; x0z[kIdxTY] = 0.; x0z[idxZ0] = 1.;

        for (istep = 0; istep < 4; ++ istep) {
            for(ipar = 0; ipar < numPars+1; ++ipar) {
                if(istep == 0) {
                    xz[ipar] = x0z[ipar];
                } else
                    xz[ipar] = x0z[ipar] + b[istep] * k1z[istep-1][ipar];
            }

            // @x'/@z = @tx/@z
            k1z[istep][kIdxX0]       = xz[kIdxTX]   * h + F_tx[istep] * qp_in * xz[idxZ0];
            k1z[istep][kIdxY0]       = xz[kIdxTY] * h + F_ty[istep] * qp_in * xz[idxZ0];
            k1z[istep][kIdxTX]   = F_tx_ty[istep] * xz[kIdxTX] + F_tx_ty[istep] * xz[kIdxTY]
                + F2_tx[istep] * xz[idxZ0];
            k1z[istep][kIdxTY] = F_ty_tx[istep] * xz[kIdxTX] + F_ty_ty[istep] * xz[kIdxTY]
                + F2_ty[istep] * xz[idxZ0];
            k1z[istep][idxZ0]        = 0.;
        }
        for (ipar = 0; ipar < numPars; ++ipar ) {
            fPropStep(ipar, kIdxZ0) = x0z[ipar] + c[rkStart]*k1z[rkStart][ipar] + c[rkMid1]*k1z[rkMid1][ipar]
                                                + c[rkMid2] *k1z[rkMid2][ipar]  + c[rkEnd] *k1z[rkEnd][ipar];
        }

        fPropStep(kIdxZ0, kIdxZ0) = 1.;
    }

    return stepFac;
}
//_____________________________________________________________________________________________________________
Bool_t SoLKalFieldStepper::FindTargetPlaneIntersection(TVector3 &intersection, 
                                                       Double_t target_z, TVector3 &dir, TVector3 &pos)
{
  Double_t delta_z = target_z - pos.Z();
  Double_t cos_theta = TMath::Abs(dir.Z());
  if (cos_theta == 0.0){
      return kFALSE;
  }else{
      Double_t t = delta_z/cos_theta;
      intersection = pos + (t*dir);
      return kTRUE;
  }
}
//______________________________________________________________________________________________________________
void SoLKalFieldStepper::PropagateStraightLine(SoLKalMatrix &stateVec, SoLKalMatrix &fPropChange, 
                             Double_t &zPos, Double_t dz)
{
  // Propagate the track state along a straight line in its current direction.
    // (x',y',z') = (x,y,z) + dz * (tx,ty,1)
    //
    // Output:
    // fPropChange: Change in propagator matrix.
    //
    // Input & output:
    // stateVec:    Current track state vector (x,y,tx,ty,qp).
    //
    // Input:
    // zPos:        Current z position of track.
    // dz:          Step length in z coordinate.

  Double_t tx = stateVec(kIdxTX, 0);
  Double_t ty = stateVec(kIdxTY, 0);

    // Update state vector.
  stateVec(kIdxX0, 0) = stateVec(kIdxX0, 0) + dz * tx;
  stateVec(kIdxY0, 0) = stateVec(kIdxY0, 0) + dz * ty;
  zPos             = zPos             + dz;
  trackPosAtZ      = zPos;

  // Update propagator matrix.
  fPropChange.UnitMatrix();

  fPropChange(kIdxX0, kIdxTX)   = dz;
  fPropChange(kIdxY0, kIdxTY) = dz;

  stepLength   = TMath::Abs(dz) * TMath::Sqrt(1. + tx*tx + ty*ty);
  trackLength += stepLength;
}
//______________________________________________________________________________________________________________
Bool_t SoLKalFieldStepper::PropagateStraightLine(SoLKalMatrix &stateVec, SoLKalMatrix &fPropChange, 
                             Double_t &zPos, const Double_t target_z, Bool_t propDir)
{
    // From the position and direction stored in the track state vector, propagate the track
    // to a target plane using a straight line. The track state and reference layer are updated
    // and the propagator matrix is calculated. The function returns the length of the straight line.
    // The class variable trackPosAtZ must contain the current z-value of the track.
    //
    // Output:
    // fPropChange: Change in propagator matrix
    //
    // Input and output:
    // stateVec:    Current Track state vector.
    // zPos:        Current z-position of track.
    //
    // Input:
    // target_z:    z coordinate of the target plane
    // propDir:     Propagation direction.

    stepLength = 0.;
    TVector3 pos(stateVec(kIdxX0, 0), stateVec(kIdxY0, 0), zPos);
    TVector3 dir;
    SoLKalTrackState::CalcDir(dir, stateVec);
    TVector3 pointIntersect;
    FindTargetPlaneIntersection(pointIntersect, target_z, dir, pos);
    Double_t dz = (pointIntersect.Z() - pos.Z());
    

    if((dz > 0. && propDir == kFALSE) || (dz < 0. && propDir == kTRUE)) {
        PropagateStraightLine(stateVec, fPropChange, zPos, dz);
        stepLength = (pos - pointIntersect).Mag();
    } else {
        fPropChange.UnitMatrix();
        if(TMath::Abs(dz) > 0.001) {
          //if(bPrintWarn) {
          //    Warning("propagateStraightLine()", Form("Track already past target plane by dz = %f.", TMath::Abs(dz)));
          //  }
            return kFALSE;
        }
    }
    return kTRUE;

}
//______________________________________________________________________________________________________________
void SoLKalFieldStepper::InitDetMaterial()
{
  //property array: effective atomic weight A, effective atomic number Z, 
  //radiation length(m), mean excitation energy (Mev), density (g/cm^3)
  fDetMatProperties[kAir][kAtomicNum] = 14.6046;
  fDetMatProperties[kAir][kProtonNum] = 7.3;
  fDetMatProperties[kAir][kExcitEnergy] = 85.7e-6;  //MeV
  fDetMatProperties[kAir][kDensity] = 1.29e-3; //g/cm^3
  fDetMatProperties[kAir][kRadLength] = 300; //m
  
  fDetMatProperties[kGEM][kAtomicNum] = 21.8;
  fDetMatProperties[kGEM][kProtonNum] = 10.603;
  fDetMatProperties[kGEM][kExcitEnergy] = 106.6e-6;  //MeV
  fDetMatProperties[kGEM][kDensity] = 0.1117433; //g/cm^3
  fDetMatProperties[kGEM][kRadLength] = 3.022; //m
  //fDetMatProperties[kGEM][kRadLength] = 2.5;
}
//_______________________________________________________________________________________________________________
Double_t SoLKalFieldStepper::CalcMultScat(SoLKalMatrix &Q, SoLKalMatrix &sv_to, Double_t length, 
                                          Double_t beta, SoLMatType type )
{
  // Add multiple scattering to the process noise covariance.
  //
  // Input and output:
  // fProc:     Process noise covariance. Contribution of multiple scattering is added to this matrix.
  //
  // Input:
  // stateVec:  Track state vector at start of an RK step.
  // length:    Track length in cm.
  // radLength: Radiation length of passed material in cm.
  // beta:      v/c of particle.
  // pid:       Geant particle ID.

  Double_t tx = sv_to(kIdxTX, 0);
  Double_t ty = sv_to(kIdxTY, 0);
  Double_t t  = 1. + tx*tx + ty*ty;
  // 1/beta^2
  Double_t beta2Inv = 1. / (beta*beta);

  // 1/momentum^2
  Double_t mom2Inv  = TMath::Power(sv_to(kIdxQP, 0), 2);

  // Squared scatter angle cms2.
  // cms = 13.6 MeV / (beta * c * p) * sqrt(l/X0) * (1 + 0.038 * ln(l/X0))
  // with l/X0 = length of particle track in units of radiation length.
  
  Double_t lx0 = length / fDetMatProperties[type][kRadLength];
  Double_t cms2 = (0.0136 * 0.0136 * beta2Inv * mom2Inv * lx0 * TMath::Power((1 + .038 * TMath::Log(lx0)),2));
  
  // Update process noise.
  Q(kIdxTX, kIdxTX) += (1 + tx*tx) * t * cms2;
  Q(kIdxTY, kIdxTY) += (1 + ty*ty) * t * cms2;
  Q(kIdxTX, kIdxTY) += tx * ty     * t * cms2;
  Q(kIdxTY, kIdxTX)    = Q(kIdxTX, kIdxTY); // matrix is symmetric

  return cms2;

}
//_____________________________________________________________________________________________________________
Double_t SoLKalFieldStepper::CalcEnergyLoss(SoLKalMatrix &Q, SoLKalMatrix &sv_to, Double_t length, 
                                            Double_t qp, Double_t beta, SoLMatType type)
{
  Double_t ZoverA = fDetMatProperties[type][kProtonNum]/fDetMatProperties[type][kAtomicNum];
  Double_t p    = fCharge / qp;
  Double_t ElossRad = 0.;
  Double_t ElossIon = 0.;
  
   if(fIsElectron) {
    // Radiation loss for electrons/positrons.
      ElossRad = CalcRadLoss(Q, length, qp, fDetMatProperties[type][kRadLength]);
      if(fIsBackward) {
      ElossRad *= -1.;
    }

    // Critical energy for gases:
    // E_c = 710 MeV / (Z + 0.92)
    // From "Passage of particles through matter", Particle Data Group, 2009
    
    
    ElossIon = length * CalcDEDXIonLepton(qp, ZoverA, fDetMatProperties[type][kDensity], 
                                          fDetMatProperties[type][kExcitEnergy]);
    
  } else { // Energy loss for heavy particles.
    // Energy loss due to ionization.
    ElossIon = length * CalcDEDXBetheBloch(beta, ZoverA, fDetMatProperties[type][kDensity],
                                           fDetMatProperties[type][kExcitEnergy]);
    
  }
  if(fIsBackward == kTRUE) {
    ElossIon *= -1.;
  }
  Double_t Eloss = ElossIon; // delta(E)
  
  if (fIsElectron){
    // For electrons: E ~ p
    // p' = p * (1 + delta(E)/p)
    // Track state parameter change is:
    // q/p' = q/p / (1 + delta(E)/p)
    
    sv_to(kIdxQP, 0) = qp / (1. +  Eloss / p);
  }
  else{
    // E' = E + delta(E)
    // E'^2 = E^2 + 2*E*delta(E) + delta(E)^2
    // p'^2 + m^2 = p^2 + m^2 + 2*E*delta(E) + delta(E)^2
    // p'^2 = p^2 + 2*E*delta(E) + delta(E)^2
    Double_t p2 = p * p;
	  Double_t E = TMath::Sqrt(p2 + fMass*fMass*1e-6);
    Double_t pnew = TMath::Sqrt(p2 + 2.*E*Eloss + Eloss*Eloss);
    if(pnew > 0.) {
      sv_to(kIdxQP, 0) = fCharge / pnew;
    }
  }
  return Eloss;
}
//________________________________________________________________________________________________________________
Double_t SoLKalFieldStepper::CalcRadLoss(SoLKalMatrix &Q, Double_t length, Double_t qp, Double_t radLength)
{
  // Calculates energy loss due to bremsstrahlung for electrons with Bethe-Heitler formula.
  // Sign will be negative if propagating in forward direction, i.e. particle looses energy.
  //
  // Input and Output:
  // fProc:     process noise matrix (will be modified)
  //
  // Input:
  // length:    track length in mm
  // mass:      particle mass in MeV/c^2
  // qp:        charge/momentum before energy loss in MeV/c
  // radLength: radiation length of material in mm.

  /*if(mass != HPhysicsConstants::mass("e-")) {
    if(bPrintWarn) {
      Warning("calcRadLoss()", Form("Particle with mass is %f MeV not an electron/positron.", mass));
    }
    return 0.;
    }*/

  Double_t t = length / radLength; // t := l/xr. Track length in units of radiation length.

  // Bethe and Heitler formula:
  // - dE/dx = E / xr
  // Integrating yields:
  // E' = E * exp(-t)
  // E'/E = exp(-t)
  // with t := l/xr
  // l = track length
  // xr = radiation length

  // delta(E) = E' - E
  //          = E * (E'/E - 1)
  //          = (E'/E - 1) / qp since E ~ p for electrons
  Double_t ElossRad = (TMath::Exp(-t) - 1.) / TMath::Abs(qp);

  // The variance of E'/E as done in:
  // D. Stampfer, M. Regler and R. Fruehwirth,
  // Track fitting with energy loss.
  // Comp. Phys. Commun. 79 (1994) 157-164
  // http://www-linux.gsi.de/~ikisel/reco/Methods/Fruehwirth-FittingWithEnergyLoss-CPC79-1994.pdf
  Double_t varElossRadFac = TMath::Exp(- t * TMath::Log(3.) / TMath::Log(2.)) - TMath::Exp(-2. * t);

  // Update process noise.
  //Q(kIdxQP, kIdxQP) += qp * qp * varElossRadFac;

  return ElossRad;

}
//______________________________________________________________________________________________________________
Double_t SoLKalFieldStepper::CalcDEDXIonLepton(Double_t qp, Double_t ZoverA, Double_t density, Double_t I)
{
  // Calculates energy loss dE/dx in MeV/mm due to ionization for relativistic electrons/positrons.
  //
  // For formula used, see:
  // Track fitting with energy loss
  // Stampfer, Regler and Fruehwirth
  // Computer Physics Communication 79 (1994), 157-164
  //
  // qp:         particle charge divided by momentum
  // ZoverA:     atomic number / atomic mass of passed material
  // density:    density of material in g/mm^3
  // I : mean excitation energy in GeV


  Double_t de       = 5.0989 * 1.e-23;                 // 4*pi*re*me*c^2 in MeV * mm^2 with re = classical electron radius
  Double_t avogadro = TMath::Na();                     // Avogadro constant in 1/mol
  Double_t ne       = ZoverA * avogadro * density/1000.;     // electron density
  Double_t me       = 0.5109989181;   // electron mass in MeV/c^2
  // gamma = E/m
  Double_t E        = 1. / TMath::Abs(qp);             // Energy for relativistic electrons.
  Double_t gamma    = E*1e3 / me;                          // Relativistic gamma-factor.

  // Formula is slightly different for electrons and positrons.
  // Factors for positrons:
  Double_t gammaFac = 3.;
  Double_t corr     = 1.95;

  Double_t dedx = 0.5 * de * ne * (2 * TMath::Log(2*me/I) + gammaFac * TMath::Log(gamma) - corr);
  
  return (- dedx);//convert from MeV/mm to GeV/m
}
//________________________________________________________________________________________________________________
Double_t SoLKalFieldStepper::CalcDEDXBetheBloch(Double_t beta, Double_t ZoverA, Double_t density, Double_t I)
{
  // Returns the ionization energy loss for thin materials with Bethe-Bloch formula.
  // - dE/dx = K/A * Z * z^2 / beta^2 * (0.5 * ln(2*me*beta^2*gamma^2*T_max/I^2) - beta^2)
  // T_max = (2*me*beta^2*gamma^2) / (1 + 2*gamma*me/mass + me^2/mass^2)
  //
  // beta:       v/c of particle
  // mass:       mass of the particle in MeV/c^2 traversing through the material
  // ZoverA:     atomic number / atomic mass of passed material
  // density:    density of passed material in g/cm^3
  // exciteEner: mean excitation energy of passed material in MeV
  // z:          atomic number of particle

  Double_t me         = 0.510998918; // electron mass in MeV/c^2
  Double_t K          = 0.307075;               // in MeV*cm^2/g for A = 1 g/mol

  Double_t Krho       = K * density; // MeV/cm
  Double_t beta2      = beta*beta;
  Double_t gamma      = 1./TMath::Sqrt(1 - beta2);
  Double_t betagamma  = beta * gamma;
  Double_t betagamma2 = betagamma * betagamma;

  if(betagamma < 0.1 || betagamma > 1000.) {
    /*if(bPrintErr && bPrintElossErr) { // Print only once.
      Error("calcDXDXBetheBloch()",
            Form("beta*gamma of particle is %f. Bethe-Bloch formula is only good for values between 0.1 and 1000.", betagamma));
            }
    bElossErr      = kTRUE;
    bPrintElossErr = kFALSE;*/
    //assert(betagamma<1000.);
    cout<<"Bethe-Bloch formula is only good for values between 0.1 and 1000."<<" beta= "<<beta<<" gamma="<<gamma<<endl;
    return 0.;
  }
  // Maximum kinetic energy that can be imparted on a free electron in a single collision.
  Double_t tmax = (2. * me * betagamma2) / (1 + 2*gamma*me/(fMass) + (me*me)/(fMass*fMass));

  return (-((Krho * fCharge*fCharge * ZoverA) / beta2)
          * (0.5 * TMath::Log((2. * me * betagamma2 * tmax)/(I*I)) - beta2))/10.;//convert from MeV/cm to GeV/m

}





















