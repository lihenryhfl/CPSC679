#ifndef SNH_H
#define SNH_H

#include "MATERIAL.h"

class SNH : public MATERIAL
{
public:
  // settings from the 2008 cubature paper
  SNH(const REAL mu = 5000, const REAL lambda = 1000);
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

protected:
  REAL _mu;
  REAL _lambda;
  REAL _alpha;
  MATRIX4 _hessJ;
};

#endif
