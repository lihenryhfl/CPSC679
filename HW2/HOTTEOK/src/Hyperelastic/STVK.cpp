#include "STVK.h"
#include <iostream>

STVK::STVK(const REAL mu, const REAL lambda) :
  _mu(mu), _lambda(lambda)
{
  _name = std::string("StVK");
}

///////////////////////////////////////////////////////////////////////
// P = first Piola-Kirchoff stress tensor
// P = F * S
///////////////////////////////////////////////////////////////////////
MATRIX STVK::PK1(const MATRIX2& F) const
{
  return F * PK2(F);
  //return PK2(F) * F;
}

///////////////////////////////////////////////////////////////////////
// S = second Piola-Kirchoff stress tensor
///////////////////////////////////////////////////////////////////////
MATRIX STVK::PK2(const MATRIX2& F) const
{
  MATRIX2 E = 0.5 * (F.transpose() * F - MATRIX2::Identity());
  MATRIX2 S = _lambda * E.trace() * MATRIX2::Identity() + 2.0 * _mu * E;
  return S;
}

///////////////////////////////////////////////////////////////////////
// compute the derivative of S, the second Piola-Kirchhoff
// with respect to E = 1/2 (F^T * F - I)
///////////////////////////////////////////////////////////////////////
MATRIX STVK::DSDE() const
{
  MATRIX result(4,4);
  result.setZero();

  result(0,0) = _lambda + 2.0 * _mu;
  result(3,0) = _lambda;

  result(1,1) = 2.0 * _mu;

  result(2,2) = 2.0 * _mu;

  result(0,3) = _lambda;
  result(3,3) = _lambda + 2.0 * _mu;

  return result;
}

///////////////////////////////////////////////////////////////////////
// take the gradient of FS (i.e. PK1) w.r.t F assuming that S is frozen
///////////////////////////////////////////////////////////////////////
MATRIX STVK::DFSDF(const MATRIX& S) const
{
  MATRIX result(4,4);
  result.setZero();

  result(0,0) = S(0,0);
  result(1,1) = S(0,0);
  result(2,2) = S(1,1);
  result(3,3) = S(1,1);

  result(2,0) = S(0,1);
  result(3,1) = S(0,1);

  result(0,2) = S(1,0);
  result(1,3) = S(1,0);

  return result;
}

///////////////////////////////////////////////////////////////////////
// derivative of PK1 w.r.t. F
///////////////////////////////////////////////////////////////////////
MATRIX STVK::hessian(const MATRIX& F) const
{
  MATRIX S = PK2(F);

  MATRIX IS = DFSDF(S);

  MATRIX diagF = blockDiag(F, 2);

  // no getting around computing this one -- changes with F
  MATRIX DEDF = computeDEDF(F);

  return IS + diagF * DSDE() * DEDF;
}

///////////////////////////////////////////////////////////////////////
// compute the derivative of E = 1/2 (F^T * F - I)
// with respect to F
///////////////////////////////////////////////////////////////////////
MATRIX STVK::computeDEDF(const MATRIX2& F) const
{
  const REAL& f00 = F(0,0);
  const REAL& f10 = F(1,0);
  const REAL& f01 = F(0,1);
  const REAL& f11 = F(1,1);
  MATRIX result(4,4);

  result(0,0) = 2.0 * f00;
  result(1,0) = f01;
  result(2,0) = f01;
  result(3,0) = 0.0;

  result(0,1) = 2.0 * f10;
  result(1,1) = f11;
  result(2,1) = f11;
  result(3,1) = 0.0;

  result(0,2) = 0.0;
  result(1,2) = f00;
  result(2,2) = f00;
  result(3,2) = 2.0 * f01;

  result(0,3) = 0.0;
  result(1,3) = f10;
  result(2,3) = f10;
  result(3,3) = 2.0 * f11;

  // E = 0.5 (F^T * F - I), so a 0.5 factor is needed in front
  result *= 0.5;

  return result;
}

///////////////////////////////////////////////////////////////////////
// repeat the given matrix "repeat" times along the block diagonal
///////////////////////////////////////////////////////////////////////
MATRIX STVK::blockDiag(const MATRIX& A, const int repeats) const
{
  const int rows = A.rows();
  const int cols = A.cols();
  MATRIX result(rows * repeats, cols * repeats);
  result.setZero();

  for (int y = 0; y < cols; y++)
    for (int x = 0; x < rows; x++)
    {
      for (int i = 0; i < repeats; i++)
        result(i * rows + x, i * cols + y) = A(x,y);
    }

  return result;
}

///////////////////////////////////////////////////////////////////////
// get the strain energy
///////////////////////////////////////////////////////////////////////
REAL STVK::psi(const MATRIX2& F) const
{
  MATRIX2 E = 0.5 * (F.transpose() * F - MATRIX2::Identity());
  return _mu * E.squaredNorm() + _lambda * 0.5 * pow(E.trace(), 2.0);
}
