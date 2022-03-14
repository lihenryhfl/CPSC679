#ifndef MATERIAL_H
#define MATERIAL_H

#include "SETTINGS.h"
#include <iostream>

class MATERIAL
{
public:
  MATERIAL() : _name("None") {};
  virtual ~MATERIAL() {};

  // P = first Piola-Kirchoff stress tensor
  // P = F * S
  virtual MATRIX PK1(const MATRIX2& F) const = 0;

  // derivative of PK1 w.r.t. F
  virtual MATRIX hessian(const MATRIX& F) const = 0;

  // get the strain energy
  virtual REAL psi(const MATRIX2& F) const = 0;

  // collision energies
  virtual REAL getEps() const = 0;
  virtual REAL cpsi(const VECTOR6& x) const = 0;
  virtual MATRIX cPK1(const VECTOR6& x) const = 0;
  virtual MATRIX cHessian(const VECTOR6& x) const = 0;

  const std::string& name() const { return _name; };

  // convert Young's modulus (E) and Poisson's ratio (nu) to Lam\'{e} parameters
  static REAL computeMu(const REAL E, const REAL nu)     { return E / (2.0 * (1.0 + nu)); };
  static REAL computeLambda(const REAL E, const REAL nu) { return (E * nu) / ((1.0 + nu) * (1.0 - 2.0 * nu)); };

protected:
  // the name of the material, for display purposes
  std::string _name;
};

#endif
