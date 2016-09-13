#ifndef ROOT_SOL_KAL_MATRIX
#define ROOT_SOL_KAL_MATRIX
//ROOT
#include "TMatrixD.h"
#include "TVector3.h"

class SoLKalMatrix : public TMatrixD {
public:
   SoLKalMatrix(Int_t rowdim = 1, Int_t coldim = 1);

   SoLKalMatrix(const SoLKalMatrix &orig);
   SoLKalMatrix(const TMatrixD &orig);

   SoLKalMatrix(TMatrixD::EMatrixCreatorsOp1 op,
              const SoLKalMatrix &prototype);
   SoLKalMatrix(TMatrixD::EMatrixCreatorsOp1 op,
              const TMatrixD &prototype);

   SoLKalMatrix(const SoLKalMatrix &a,
              TMatrixD::EMatrixCreatorsOp2 op,
              const SoLKalMatrix &b) ;
   SoLKalMatrix(const TMatrixD &a,
              TMatrixD::EMatrixCreatorsOp2 op,
              const TMatrixD &b) ;

   SoLKalMatrix(const TVector3 &v);

   virtual ~SoLKalMatrix() {}

   virtual void      DebugPrint(Option_t *opt = "", Int_t nc = 5) const;

   static SoLKalMatrix ToKalMat  (const TVector3 &vec);
   static TVector3   ToThreeVec(const TMatrixD &mat);

private:

   ClassDef(SoLKalMatrix,1)      // Base class for Kalman matrix

};

#endif
