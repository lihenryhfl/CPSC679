#ifndef STVK_H
#define STVK_H

#include "MATERIAL.h"

class STVK : public MATERIAL
{
public:
  // settings from the 2008 cubature paper
  STVK(const REAL mu = 5000, const REAL lambda = 1000);
  ~STVK() {};

  // P = first Piola-Kirchoff stress tensor
  // P = F * S
  virtual MATRIX PK1(const MATRIX2& F) const override;

  // S = second Piola-Kirchoff stress tensor
  MATRIX PK2(const MATRIX2& F) const;
  
  // derivative of PK2 w.r.t E = 0.5 * (F^T * F - I)
  MATRIX DSDE() const;

  // derivative of PK1 w.r.t. F
  virtual MATRIX hessian(const MATRIX& F) const override;

  // get the strain energy
  virtual REAL psi(const MATRIX2& F) const override;

protected:
  // take the gradient of FS (i.e. PK1) w.r.t F assuming that S is frozen
  // \frac{\partial F S}{\partial F}
  MATRIX DFSDF(const MATRIX& S) const;

  // repeat the given matrix "repeats" times along the block diagonal
  MATRIX blockDiag(const MATRIX& A, const int repeats) const;

  // compute the derivative of E = 1/2 (F^T * F - I)
  // with respect to F
  MATRIX computeDEDF(const MATRIX2& F) const;

  REAL _mu;
  REAL _lambda;
};

#endif
