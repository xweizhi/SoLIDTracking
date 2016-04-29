#ifndef ROOT_SoLID_ECAL
#define ROOT_SoLID_ECAL
//c++
#include <map>
#include <vector>
//ROOT
#include "TRandom.h"
//Hall A Analyzer
#include "THaDetMap.h"
#include "THaSubDetector.h"

using namespace std;

enum ECalDataType{ kLAECPos = 0, kLAECEdp, kFAECPos, kFAECEdp};
enum ECalMapKey { kECalXPos = 0, kECalYPos, kECalEdp };
class SoLIDECal : public THaSubDetector
{
  public:
  SoLIDECal();
  SoLIDECal( const char* name, const char* description = "",
           THaDetectorBase* parent = 0 );
  ~SoLIDECal();
  
  virtual void      Clear( Option_t* opt="" );
  virtual Int_t     Decode( const THaEvData& );
  virtual EStatus   Init( const TDatime& date );
  virtual void      Print( Option_t* opt="" ) const;
  virtual void      PrintDataBase() const;
  virtual Int_t     Begin( THaRunBase* r=0 );
  virtual Int_t     End( THaRunBase* r=0 );
  
  Bool_t  IsFAECTriggered() const { return fIsFAECTriggered; }
  Bool_t  IsLAECTriggered() const { return fIsLAECTriggered; }
  UInt_t  GetNLAECHits()    const { return fNLAECHits; }
  UInt_t  GetNFAECHits()    const { return fNFAECHits; }
  map< Int_t, vector<Float_t> > * GetLAECHits()  { return &fLAECHitMap; }
  map< Int_t, vector<Float_t> > * GetFAECHits()  { return &fFAECHitMap; }                
  
  void    SmearPosition(Float_t *x, Float_t *y);
  void    SmearEnergy(Float_t *energy);
  protected:
  virtual Int_t     ReadDatabase( const TDatime& date );
  virtual Int_t     DefineVariables( EMode mode = kDefine );
  
  Bool_t                          fIsLAECTriggered;
  Bool_t                          fIsFAECTriggered;
  UInt_t                          fNLAECHits;
  UInt_t                          fNFAECHits;
  Double_t                        fLAECZ;
  Double_t                        fFAECZ;
  Double_t                        fLAECEdpCut;
  Double_t                        fFAECEdpCut;
  Double_t                        fPosReso;
  Double_t                        fEReso;
  map< Int_t, vector<Float_t> >   fLAECHitMap;
  map< Int_t, vector<Float_t> >   fFAECHitMap;
  
  ClassDef(SoLIDECal,0)
};

#endif