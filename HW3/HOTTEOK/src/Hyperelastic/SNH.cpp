#include "SNH.h"
#include <iostream>
using namespace std;

SNH::SNH(const REAL mu, const REAL lambda) :
  _mu(mu), _lambda(lambda), _alpha(1 - mu / lambda)
{
  _name = std::string("SNH");
  _hessJ = makeHessJ();
  _R = makeR();
}

MATRIX2 SNH::makeR() const
{
  // computes a 2x2 rotation matrix for 90 deg rotation
  MATRIX2 R;
  R.setZero();
  R(0, 1) = -1.0;
  R(1, 0) = 1.0;
  return R;
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

// copy the values of B into a slice of A, i.e.
// A[xstart:xend, ystart:yend] = B
void SNH::assign(MATRIX& A, int xstart, int xend, int ystart, int yend, MATRIX& B, int bxstart, int bystart) const
{
  int nr = xend - xstart, nc = yend - ystart;
  int nrb = B.rows() - bxstart, ncb = B.cols() - bystart;

  if ((nr > nrb) || (nc > ncb) || (nr > A.rows()) || (nc > A.cols())) {
    cout << (nr > nrb) << (nc > ncb) << (nr > A.rows()) << (nc > A.cols()) << endl;
    cout << "nr: " << nr << " nrb: " << nrb << endl;
    cout << "nc: " << nc << " ncb: " << ncb << endl;
    cout << "A.rows(): " << A.rows() << " A.cols(): " << A.cols() << endl;
    throw std::invalid_argument("SNH::assign: A slice does not map B shape!");
  }

  for (int j = ystart; j < yend; j++)
    for (int i = xstart; i < xend; i++) {
      A(i,j) = B(i-xstart+bxstart,j-ystart+bystart);
    }
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
  assign(A1, 0, 4, 0, 1, C);
  B << 0, a[0], a[0], 2 * a[1];
  B = _R * B;
  C = reshape(B, 4, 1);
  assign(A1, 0, 4, 1, 2, C);
  A1 = A1 * pow(norm, -3);


  MATRIX A2 = MATRIX(4, 2);
  B = a * a.transpose();
  B = _R * B;
  C = reshape(B, 4, 1) * a[0];
  assign(A2, 0, 4, 0, 1, C);
  C = reshape(B, 4, 1) * a[1];
  assign(A2, 0, 4, 1, 2, C);
  A2 = A2 * -3 * pow(norm, -5);

  MATRIX A3 = MATRIX(4, 2);
  B = MATRIX2::Identity();
  B = _R * B;
  C = reshape(B, 4, 1) * a[0];
  assign(A3, 0, 4, 0, 1, C);
  C = reshape(B, 4, 1) * a[1];
  assign(A3, 0, 4, 1, 2, C);
  A3 = A3 * pow(norm, -3);

  MATRIX A = A1 + A2 + A3;
  A = reshape(A, 2, 4);

  MATRIX A_1 = MATRIX(1, 4), A_2 = MATRIX(1, 4);
  assign(A_1, 0, 1, 0, 4, A, 0, 0);
  assign(A_2, 0, 1, 0, 4, A, 1, 0);
  A_1 = reshape(A_1, 2, 2);
  A_2 = reshape(A_2, 2, 2);
  MATRIX nA_1 = -A_1;
  MATRIX nA_2 = -A_2;

  assign(H, 2, 4, 2, 4, nA_1);
  assign(H, 8, 10, 2, 4, nA_2);
  assign(H, 2, 4, 4, 6, A_1);
  assign(H, 8, 10, 4, 6, A_2);

  assign(H, 4, 6, 2, 4, A_1);
  assign(H, 10, 12, 2, 4, A_2);
  assign(H, 4, 6, 4, 6, nA_1);
  assign(H, 10, 12, 4, 6, nA_2);

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

