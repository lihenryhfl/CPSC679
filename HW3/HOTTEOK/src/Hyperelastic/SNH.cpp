#include "SNH.h"
#include <iostream>
using namespace std;

SNH::SNH(const REAL mu, const REAL lambda, const REAL eps) :
  _eps(eps), _mu(mu), _lambda(lambda), _alpha(1 - mu / lambda)
{
  _name = std::string("SNH");
  _hessJ = makeHessJ();
  _R << 0., 1., -1., 0.;
}

REAL SNH::getEps() const
{
  return _eps;
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
  MATRIX2 pJpF;
  pJpF(0, 0) = F(1, 1);
  pJpF(0, 1) = -F(1, 0);
  pJpF(1, 0) = -F(0, 1);
  pJpF(1, 1) = F(0, 0);
  return pJpF;
}

MATRIX SNH::dtdx(const VECTOR6& x) const
{
  MATRIX dtdx = MATRIX(2, 6);
  dtdx.setZero();
  dtdx(0, 0) = 1;
  dtdx(1, 1) = 1;
  dtdx(0, 2) = -1;
  dtdx(1, 3) = -1;
  return dtdx;
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

MATRIX SNH::reshape(const MATRIX& A, int n0, int n1) const
{
  int num_entries = n0 * n1;

  if (num_entries != A.rows() * A.cols())
    throw std::invalid_argument("SNH::reshape: Num entries in A should be equal to n0 * n1.");

  MATRIX newA = MATRIX(n0, n1);

  int nx, ny, x, y;
  for (int i = 0; i < num_entries; i++) {
    nx = i / n1;
    ny = i % n1;
    x = i / A.cols();
    y = i % A.cols();
    newA(nx, ny) = A(x, y);
  }

  return newA;
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

///////////////////////////////////////////////////////////////////////
// collision energy
///////////////////////////////////////////////////////////////////////
REAL SNH::cpsi(const VECTOR6& x) const
{
  VECTOR2 x0, x1, x2, n, t;
  x0 << x[0], x[1];
  x1 << x[2], x[3];
  x2 << x[4], x[5];
  n = normal(x);
  t = x0 - x1;
  return _mu * pow((n.transpose() * t - _eps), 2);
}

///////////////////////////////////////////////////////////////////////
// first Piola-Kirchoff tensor of collision energy
///////////////////////////////////////////////////////////////////////
MATRIX SNH::cPK1(const VECTOR6& x) const
{
  MATRIX g = MATRIX(1, 6);
  VECTOR2 x0, x1, x2, n, t;
  x0 << x[0], x[1];
  x1 << x[2], x[3];
  x2 << x[4], x[5];
  n = normal(x);
  t = x0 - x1;
  g = n.transpose() * dtdx(x) + t.transpose() * normalJacobian(x);
  return 2 * _mu * (n.transpose() * t - _eps) * g;
}

///////////////////////////////////////////////////////////////////////
// hessian of collision energy
///////////////////////////////////////////////////////////////////////
MATRIX SNH::cHessian(const VECTOR6& x) const
{
  MATRIX H = MATRIX(6, 6);
  MATRIX g, JnJt;
  MATRIX nHessian = normalHessian(x);
  VECTOR2 x0, x1, x2, n, t;
  x0 << x[0], x[1];
  x1 << x[2], x[3];
  x2 << x[4], x[5];
  n = normal(x);
  t = x0 - x1;
  g = n.transpose() * dtdx(x) + t.transpose() * normalJacobian(x);
  JnJt = dtdx(x).transpose() * normalJacobian(x);

  // compute t * nHessian tensor product
  nHessian = t.transpose() * reshape(nHessian, 2, 36);
  nHessian = reshape(nHessian, 6, 6);

  H = nHessian + JnJt + JnJt.transpose();
  H = H * (t.transpose() * n - _eps) + g.transpose() * g;
  return 2 * _mu * H;
}

///////////////////////////////////////////////////////////////////////
// hessian of the normal
///////////////////////////////////////////////////////////////////////
MATRIX SNH::normalHessian(const VECTOR6& x) const
{
  MATRIX H = MATRIX(12, 6);
  H.setZero();

  VECTOR2 x1, x2, a;
  x1 << x[2], x[3];
  x2 << x[4], x[5];
  a = x2 - x1;
  REAL norm = a.norm();

  // temp matrix used to fill in entries of A
  MATRIX B = MATRIX(2, 2), C = MATRIX(4, 1);

  MATRIX A1 = MATRIX(4, 2);
  B << 2 * a[0], a[1], a[1], 0;
  B = _R * B;
  C = reshape(B, 4, 1);
  A1.block(0, 0, 4, 1) = C;
  B << 0, a[0], a[0], 2 * a[1];
  B = _R * B;
  C = reshape(B, 4, 1);
  A1.block(0, 1, 4, 1) = C;
  A1 = A1 * pow(norm, -3);


  MATRIX A2 = MATRIX(4, 2);
  B = a * a.transpose();
  B = _R * B;
  C = reshape(B, 4, 1) * a[0];
  A2.block(0, 0, 4, 1) = C;
  C = reshape(B, 4, 1) * a[1];
  A2.block(0, 1, 4, 1) = C;
  A2 = A2 * -3 * pow(norm, -5);

  MATRIX A3 = MATRIX(4, 2);
  B = MATRIX2::Identity();
  B = _R * B;
  C = reshape(B, 4, 1) * a[0];
  A3.block(0, 0, 4, 1) = C;
  C = reshape(B, 4, 1) * a[1];
  A3.block(0, 1, 4, 1) = C;
  A3 = A3 * pow(norm, -3);

  MATRIX A = A1 + A2 + A3;
  A = reshape(A, 2, 4);

  MATRIX A_1 = A.block(0, 0, 1, 4), A_2 = A.block(1, 0, 1, 4);
  A_1 = reshape(A_1, 2, 2);
  A_2 = reshape(A_2, 2, 2);

  H.block(2, 2, 2, 2) = -A_1;
  H.block(8, 2, 2, 2) = -A_2;
  H.block(2, 4, 2, 2) = A_1;
  H.block(8, 4, 2, 2) = A_2;

  H.block(4, 2, 2, 2) = A_1;
  H.block(10, 2, 2, 2) = A_2;
  H.block(4, 4, 2, 2) = -A_1;
  H.block(10, 4, 2, 2) = -A_2;

  return H;
}

///////////////////////////////////////////////////////////////////////
// jacobian of the normal
///////////////////////////////////////////////////////////////////////
MATRIX SNH::normalJacobian(const VECTOR6& x) const
{
  MATRIX dndx = MATRIX(2, 6);
  dndx.setZero();

  VECTOR2 x1, x2, a;
  x1 << x[2], x[3];
  x2 << x[4], x[5];
  a = x2 - x1;

  REAL norm = a.norm();
  MATRIX2 C = a * a.transpose() / pow(norm, 3) - MATRIX2::Identity() / norm;
  MATRIX2 RC = _R * C;

  dndx(0, 2) = RC(0, 0);
  dndx(0, 3) = RC(0, 1);
  dndx(1, 2) = RC(1, 0);
  dndx(1, 3) = RC(1, 1);

  dndx(0, 4) = -RC(0, 0);
  dndx(0, 5) = -RC(0, 1);
  dndx(1, 4) = -RC(1, 0);
  dndx(1, 5) = -RC(1, 1);

  return dndx;
}

///////////////////////////////////////////////////////////////////////
// normal
///////////////////////////////////////////////////////////////////////
VECTOR2 SNH::normal(const VECTOR6& x) const
{
  VECTOR2 x1, x2, n;
  x1 << x[2], x[3];
  x2 << x[4], x[5];
  n = _R * (x2 - x1);
  n = n / n.norm();
  return n;
}

