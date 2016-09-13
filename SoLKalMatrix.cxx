//c++
#include <iostream>
#include <iomanip>
//ROOT
#include "TString.h"
//SoLIDTracking
#include "SoLKalMatrix.h"

ClassImp(SoLKalMatrix)
//___________________________________________________________
SoLKalMatrix::SoLKalMatrix(Int_t rowdim, Int_t coldim)
              :TMatrixD(rowdim, coldim)
{
}
//___________________________________________________________
SoLKalMatrix::SoLKalMatrix(const SoLKalMatrix &orig)
              :TMatrixD(orig)
{
}
//___________________________________________________________
SoLKalMatrix::SoLKalMatrix(const TMatrixD &orig)
              :TMatrixD(orig)
{
}
//___________________________________________________________
SoLKalMatrix::SoLKalMatrix(TMatrixD::EMatrixCreatorsOp1 op,
                             const SoLKalMatrix &prototype)
              :TMatrixD(op, prototype)
{
}
//___________________________________________________________
SoLKalMatrix::SoLKalMatrix(TMatrixD::EMatrixCreatorsOp1 op,
                             const TMatrixD &prototype)
              :TMatrixD(op, prototype)
{
}
//___________________________________________________________
SoLKalMatrix::SoLKalMatrix(const SoLKalMatrix &a,
                             TMatrixD::EMatrixCreatorsOp2 op,
                             const SoLKalMatrix &b)
              :TMatrixD(a, op, b)
{
}
//___________________________________________________________
SoLKalMatrix::SoLKalMatrix(const TMatrixD &a,
                             TMatrixD::EMatrixCreatorsOp2 op,
                             const TMatrixD &b)
              :TMatrixD(a, op, b)
{
}
//___________________________________________________________
SoLKalMatrix::SoLKalMatrix(const TVector3 &v)
              :TMatrixD(SoLKalMatrix::ToKalMat(v))
{
}
//___________________________________________________________
void SoLKalMatrix::DebugPrint(Option_t *opt, Int_t ncolsps) const
{
   using namespace std;

   Int_t ncols   = GetNcols();
   Int_t nrows   = GetNrows();
   Int_t nsheets = (ncols-1)/ncolsps + 1;
   Int_t lastcw  = ncols - ncolsps*(nsheets-1);
   Int_t ns      = 0;

   TString title(opt);
   TString off(title.Length());
   for (Int_t c=0; c<title.Length(); c++) off += " ";

   for (Int_t sc = 1; sc <= ncols; sc += ncolsps) {
      ns++;
      if (ns == 1) cerr << title << "+-";
      if (ns == nsheets) {
         for (Int_t ib = 1; ib <= 11*lastcw; ib++) cerr << " ";
         cerr << "-+" << endl;
      } else {
         cerr << endl;
      }
      for (Int_t i = 1; i <= nrows; i++) {
         if (ns == 1) cerr << off << "| ";
         for (Int_t j = sc; j < sc+ncolsps && j <= ncols; j++) {
            cerr // << setiosflags(ios::scientific) 
                 << setw(11) << setprecision(4)
                 << (*this)(i-1,j-1);
            if (j == ncols) cerr << " |";
         }
         cerr << endl;
      }
      if (ns == 1) cerr << off << "+-";
      if (ns == nsheets) {
         for (Int_t ib = 1; ib <= 11*lastcw; ib++) cerr << " ";
         cerr << "-+" << endl;
      } else {
         cerr << endl;
      }
   }
}
//___________________________________________________________
SoLKalMatrix SoLKalMatrix::ToKalMat(const TVector3 &a)
{
   SoLKalMatrix md(3,1);
   md(0,0) = a.x();
   md(1,0) = a.y();
   md(2,0) = a.z();
   return md;
}
//___________________________________________________________
TVector3 SoLKalMatrix::ToThreeVec(const TMatrixD &a)
{
   TVector3 v;
   v.SetXYZ(a(0,0), a(1,0),a(2,0));
   return v;
}
//___________________________________________________________










