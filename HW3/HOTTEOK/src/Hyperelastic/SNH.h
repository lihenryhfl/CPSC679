#ifndef SNH_H
#define SNH_H

#include "MATERIAL.h"

class SNH : public MATERIAL
{
public:
  // settings from the 2008 cubature paper
  SNH(const REAL mu = 5000, const REAL lambda = 1000, const REAL eps = 0.02);
  ~SNH() {};

  // P = first Piola-Kirchoff stress tensor
  // P = F * S
  virtual MATRIX PK1(const MATRIX2& F) const override;

  VECTOR4 flatten(const MATRIX2& A) const;

  // compute partial J / partial F
  MATRIX2 DJDF(const MATRIX2& F) const;

  MATRIX4 makeHessJ() const;

  // derivative of PK1 w.r.t. F
  virtual MATRIX hessian(const MATRIX& F) const override;

  // get the strain energy
  virtual REAL psi(const MATRIX2& F) const override;

  MATRIX reshape(const MATRIX& A, int n0 = 1, int n1 = 1) const;

  // collision energies
  MATRIX dtdx(const VECTOR6& x) const;
  REAL cpsi(const VECTOR6& x) const override;
  MATRIX cPK1(const VECTOR6& x) const override;
  MATRIX cHessian(const VECTOR6& x) const override;

  // normals
  MATRIX normalHessian(const VECTOR6& x) const;
  MATRIX normalJacobian(const VECTOR6& x) const;
  VECTOR2 normal(const VECTOR6& x) const;

  // unprotected object variables
  REAL getEps() const override;

protected:
  REAL _eps;
  REAL _mu;
  REAL _lambda;
  REAL _alpha;
  MATRIX4 _hessJ;
  MATRIX2 _R;
};

#endif
