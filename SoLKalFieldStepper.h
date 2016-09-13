#ifndef ROOT_SOL_KAL_FIELD_STEPPER
#define ROOT_SOL_KAL_FIELD_STEPPER
//c++
//ROOT
#include "TVector3.h"
//SoLIDTracking
#include "SoLIDUtility.h"
#include "SoLIDFieldMap.h"
#include "SoLKalMatrix.h"
class SoLKalTrackSystem;
class SoLKalTrackSite;
class SoLKalTrackState;

class SoLKalFieldStepper
{
  public:
  ~SoLKalFieldStepper();
  static SoLKalFieldStepper * GetInstance() {
    if (fSoLKalFieldStepper == NULL) fSoLKalFieldStepper = new SoLKalFieldStepper();
    return fSoLKalFieldStepper;
  }
  void TurnOnMS()       { fIsMSOn = kTRUE;  }
  void TurnOffMS()      { fIsMSOn = kFALSE; }
  void TurnOnDEDX ()    { fIsDEDXOn = kTRUE;  }
  void TurnOffDEDX ()   { fIsDEDXOn = kFALSE; }
  
  inline Bool_t IsMSOn()   const { return fIsMSOn;     }
  inline Bool_t IsDEDXOn() const { return fIsDEDXOn;   }
  inline Double_t GetEnergyLoss()  const { return energyLoss; }
  inline Double_t GetTrackLength() const { return trackLength; }
  inline Double_t GetTrackPosAtZ() const { return trackPosAtZ; }
  
  void Transport(const SoLKalTrackSite  &from, // site from
                       SoLKalTrackSite   &to,   // sit to
                       SoLKalMatrix      &sv,   // state vector
                       SoLKalMatrix      &F,    // propagator matrix
                       SoLKalMatrix      &Q);   // process noise matrix
                        
  void Transport(const SoLKalTrackSite  &from, // site from
                       Double_t         &finalZ, // z position of the destination
                       SoLKalMatrix       &sv,   // state vector
                       SoLKalMatrix       &F,    // propagator matrix
                       SoLKalMatrix       &Q);   // process noise matrix
  void Transport(const SoLKalTrackState &sv_from,
                       Double_t         &finalZ, // z position of the destination
                       SoLKalMatrix       &sv,   // state vector
                       SoLKalMatrix       &F,    // propagator matrix
                       SoLKalMatrix       &Q);   // process noise matrix
                       
  Double_t RKPropagation(SoLKalMatrix &stateVec, SoLKalMatrix &fPropStep, Double_t stepSize, 
                         Bool_t bCalcJac, Bool_t dir);
  void PropagationClassicalRK4(TVector3 &inMom, TVector3 &inPos, Double_t &finalZ, Double_t &charge, 
                               Double_t &stepSize, TVector3 &fiMom, TVector3 &fiPos);
  void RightHandSide(const Double_t y[], const Double_t charge, 
                     const Double_t mom_mag, Double_t dydx[]);
  Double_t Distance2Points(const TVector3 &vec1, const TVector3 &vec2);
  Bool_t FindTargetPlaneIntersection(TVector3 &intersection, 
                                     Double_t target_z, TVector3 &dir, TVector3 &pos);
  void PropagateStraightLine(SoLKalMatrix &stateVec, SoLKalMatrix &fPropChange, 
                             Double_t &zPos, Double_t dz);
  Bool_t PropagateStraightLine(SoLKalMatrix &stateVec, SoLKalMatrix &fPropChange, 
                               Double_t &zPos, const Double_t target_z, Bool_t propDir);
  void InitTrack(Double_t &mass, Double_t &charge, Bool_t &isElectron, Bool_t dir);
  void UseDefaultStep();
  void UseFineStep();
  
  //Process Noise calculation
  Double_t CalcMultScat(SoLKalMatrix &Q, SoLKalMatrix &sv_to, Double_t length, 
                        Double_t beta, SoLMatType type = kAir);
  Double_t CalcEnergyLoss(SoLKalMatrix &Q, SoLKalMatrix &sv_to, Double_t length, 
                          Double_t qp, Double_t beta, SoLMatType type = kAir);
  Double_t CalcRadLoss(SoLKalMatrix &Q, Double_t length, Double_t qp, Double_t radLength);
  Double_t CalcDEDXIonLepton(Double_t qp, Double_t ZoverA, Double_t density, Double_t I);
  Double_t CalcDEDXBetheBloch(Double_t beta, Double_t ZoverA, Double_t density, Double_t I);
  protected:
  SoLKalFieldStepper();
  
  void InitDetMaterial();
  SoLIDFieldMap* fFieldMap;
  
  static SoLKalFieldStepper* fSoLKalFieldStepper;
  
  Bool_t    fIsMSOn;         //! switch for multiple scattering
  Bool_t    fIsDEDXOn;       //! switch for energy loss
  Bool_t    fIsBackward;     //kTRUE: backward propagation; kFALSE: forward propagation
  Double_t  trackPosAtZ;  //z position of the track
  Double_t  minPrecision; // minimum precision for the Runge-Kutta method
  Double_t  maxPrecision; // maximum precision for the Runge-Kutta method
  Double_t  minStepSize;  //minimum step size for the Runge-Kutta method
  Double_t  maxStepSize;  //maximum step size for the Runge-Kutta method
  Double_t  maxNumSteps;  //maximum number of steps for the Runge-Kutta method
  Double_t  stepSizeDec;  //step size reduction factor for the Runge-Kutta method
  Double_t  stepSizeInc;  //step size increasment factor for the Runge-Kutta method
  Double_t  stepLength;   //distance between the Runge-Kutta start and end points
  Double_t  trackLength;  //total track length
  Double_t  energyLoss;   //total energy loss
  Double_t  maxDist;      //maximum distance for straight line approximation
  Int_t     jstep;        //Runge-Kutta step number
  Double_t  initialStepSize; //initial step size for the Runge-Kutta method
  Double_t  minLengthCalcQ; //when the step size of propagation is larger than this value, process noise will be calculated
  
  Double_t  fMass;
  Double_t  fCharge;
  Bool_t    fIsElectron;
  
  Double_t  fDetMatProperties[2][5];
  
};

#endif
