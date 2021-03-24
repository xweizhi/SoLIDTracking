#ifndef ROOT_SoLID_Utility
#define ROOT_SoLID_Utility

#include "THashTable.h"
#include "TObject.h"
#include <iostream>
#include <algorithm>
#include <cassert>
using namespace std;
#define kSdim 5
#define kMdim 2
#define kGiga 1.e9
enum ECType{kLAEC = 0, kFAEC};
enum SoLKalIndex{ kIdxX0 = 0, kIdxY0, kIdxTX, kIdxTY, kIdxQP, kIdxZ0};
enum SoLMatType{kAir = 0, kGEM};
enum SeedType{kMidBack = 0, kFrontMid, kFrontBack, kTriplet};
enum SoLMatPropertyIdx{ kAtomicNum = 0, kProtonNum, kExcitEnergy, kDensity, kRadLength }; 
static const double kElectronMass = 0.51099891;
static const double kPimMass = 139.57018;
static const double kProtonMass = 938.2720813;


class CSpair : public TObject {
public:
  CSpair( UShort_t crate, UShort_t slot ) : fCrate(crate), fSlot(slot) {}
  virtual ULong_t Hash() const {
    UInt_t cs = (static_cast<UInt_t>(fCrate)<<16) + static_cast<UInt_t>(fSlot);
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,16,0)
    return TString::Hash( &cs, sizeof(cs) );
#else
    return TMath::Hash( &cs, sizeof(cs) );
#endif
  }
  virtual Bool_t IsEqual( const TObject* obj ) const {
    const CSpair* m;
    if( !obj or !(m = dynamic_cast<const CSpair*>(obj)) ) return kFALSE;
    return ( fCrate == m->fCrate and fSlot == m->fSlot );
  }
  UShort_t  fCrate;
  UShort_t  fSlot;
};

class DAQmodule : public CSpair {
public:
  DAQmodule( UShort_t crate, UShort_t slot, UInt_t model, UInt_t nchan )
    : CSpair(crate, slot), fModel(model), fNchan(nchan),
      fHasResolution(false) {}
  DAQmodule( UShort_t crate, UShort_t slot, UInt_t model, UInt_t nchan,
             Double_t res )
    : CSpair(crate, slot), fModel(model), fNchan(nchan),
      fResolution(res*1e-12) /* NB: resolution in ps! */,
      fHasResolution(true) {}
  virtual ~DAQmodule() {}
  virtual void Copy( TObject& obj ) const {
    TObject::Copy(obj);
    assert( dynamic_cast<DAQmodule*>(&obj) );
    DAQmodule* m = static_cast<DAQmodule*>(&obj);
    m->fCrate = fCrate; m->fSlot = fSlot; m->fModel = fModel;
    m->fNchan = fNchan; m->fResolution = fResolution;
    m->fHasResolution = fHasResolution;
  }
  virtual void Print( Option_t* ) const {
    cout << "DAQmodule: "
         << " crate = " << fCrate
         << " slot = "  << fSlot
         << " model = " << fModel
         << " nchan = " << fNchan;
    if( fHasResolution )
      cout << " res = "   << fResolution*1e12 << " ps";
    cout << endl;
  }
  UInt_t    fModel;
  UInt_t    fNchan;
  Double_t  fResolution;
  bool      fHasResolution;
};



struct DeleteObject {
  template< typename T >
  void operator() ( const T* ptr ) const { delete ptr; }
};

//___________________________________________________________________________
template< typename Container >
inline void DeleteContainer( Container& c )
{
  // Delete all elements of given container of pointers
  for_each( c.begin(), c.end(), DeleteObject() );
  c.clear();
}
//___________________________________________________________________________
template< typename ContainerOfContainers >
inline void DeleteContainerOfContainers( ContainerOfContainers& cc )
{
  // Delete all elements of given container of containers of pointers
  for_each( cc.begin(), cc.end(),
          DeleteContainer<typename ContainerOfContainers::value_type> );
  cc.clear();
}
//___________________________________________________________________________
inline Int_t NumberOfSetBits( UInt_t v )
{
  // Count number of bits set in 32-bit integer. From
  // http://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetParallel
  v = v - ((v >> 1) & 0x55555555);
  v = (v & 0x33333333) + ((v >> 2) & 0x33333333);
  return (((v + (v >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
}

//_____________________________________________________________________________
class SoLIDCaloHit : public TObject{
  public:
  SoLIDCaloHit() {}
  SoLIDCaloHit(Double_t xpos, Double_t ypos, Int_t plane, Double_t edp, Int_t subid) :
  fXPos(xpos), fYPos(ypos), fECID(plane), fEdp(edp), fSubDetID(subid) {}
  ~SoLIDCaloHit() {}
  Double_t fXPos;
  Double_t fYPos;
  Int_t    fECID;
  Double_t fEdp;
  Int_t    fSubDetID;
  
  ClassDef(SoLIDCaloHit, 1)
};
ClassImp(SoLIDCaloHit)
#endif















