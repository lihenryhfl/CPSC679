#include "SNH.h"
#include <iostream>
using namespace std;

SNH::SNH(const REAL mu, const REAL lambda) :
  _mu(mu), _lambda(lambda), _alpha(1 - mu / lambda)
{
  _name = std::string("SNH");
  _hessJ = makeHessJ();
}

MATRIX4 SNH::makeHessJ() const
{
  // computes the cofactor matrix of F
  MATRIX4 hessJ;
  hessJ.setZero();
  hessJ(3, 0) = 1.0;
  hessJ(2, 1) = -1.0;
  hessJ(1, 2) = -1.0;
  hessJ(0, 3) = 1.0;
  return hessJ;
}

MATRIX2 SNH::DJDF(const MATRIX2& F) const
{
  // computes the cofactor matrix of F
  MATRIX2 pJpF;
  pJpF(0, 0) = F(1, 1);
  pJpF(0, 1) = -F(1, 0);
  pJpF(1, 0) = -F(0, 1);
  pJpF(1, 1) = F(0, 0);
  return pJpF;
}

VECTOR4 SNH::flatten(const MATRIX2& A) const
{
  VECTOR4 column;

  unsigned int index = 0;
  for (unsigned int j = 0; j < A.cols(); j++)
    for (unsigned int i = 0; i < A.rows(); i++, index++)
      column[index] = A(i,j);

  return column;
}

///////////////////////////////////////////////////////////////////////
// derivative of PK1 w.r.t. F
///////////////////////////////////////////////////////////////////////
MATRIX SNH::hessian(const MATRIX& F) const
{
  const VECTOR4 pJpF = flatten(DJDF(F));
  REAL scale = _lambda * (F.determinant() - _alpha);
  return _mu * MATRIX4::Identity() + _lambda * pJpF * pJpF.transpose() + scale * _hessJ;
}

///////////////////////////////////////////////////////////////////////
// get the strain energy
///////////////////////////////////////////////////////////////////////
REAL SNH::psi(const MATRIX2& F) const
{
  const REAL Ic = F.squaredNorm();
  const REAL Jminus1 = F.determinant() - _alpha;
  const REAL spring_term = _mu * (Ic - 3.0);
  const REAL volume_term = _lambda * Jminus1 * Jminus1;
  return 0.5 * (spring_term + volume_term);
}

///////////////////////////////////////////////////////////////////////
// P = first Piola-Kirchoff stress tensor
// P = F * S
///////////////////////////////////////////////////////////////////////
MATRIX SNH::PK1(const MATRIX2& F) const
{
  const MATRIX2 pJpF = DJDF(F);
  const REAL Jminus1 = F.determinant() - _alpha;
  MATRIX2 _PK1 = _mu * F + _lambda * Jminus1 * pJpF;
  return _PK1;
}

