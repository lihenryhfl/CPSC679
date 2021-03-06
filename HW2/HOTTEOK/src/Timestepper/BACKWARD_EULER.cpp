#include "BACKWARD_EULER.h"
#include "Hyperelastic/STVK.h"
#include <float.h>
#include <iostream>

using namespace std;

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
BACKWARD_EULER::BACKWARD_EULER(TRIANGLE_MESH& triangleMesh, MATERIAL& hyperelastic) :
  TIMESTEPPER(triangleMesh, hyperelastic)
{
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// clang sure cares that this has to be here
///////////////////////////////////////////////////////////////////////////////////////////////////////
BACKWARD_EULER::~BACKWARD_EULER()
{
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
bool BACKWARD_EULER::solve(const bool verbose)
{
  cout << "========================================================= " << endl;
  cout << "  Backward Euler step: " << _totalSteps << endl;
  cout << "========================================================= " << endl;

  _position = _triangleMesh.getPosition();

  _triangleMesh.computeFs();
  VECTOR fInternal = _triangleMesh.computeMaterialForces(&_hyperelastic);
  cout << " R norm: " << fInternal.norm() << " _dt: " << _dt << endl;


  //VECTOR update = (_dt * _dt) * _Minv * (fInternal + _externalForces + _C * _velocity) + _dt * _velocity;

  //MATRIX dfdx = _triangleMesh.computeStiffnessMatrix(&_hyperelastic);
  //MATRIX A = _M - (_dt * _dt) * dfdx;
  //VECTOR b = (_dt * _dt) * (fInternal + _externalForces) + _dt * _M * _velocity;
  //VECTOR update = A.colPivHouseholderQr().solve(b);

  // build filter matrix (S matrix from lecture)
  const int DOFs = _triangleMesh.DOFs();
  MATRIX filter(DOFs, DOFs);
  filter.setIdentity();
  VECTOR z(DOFs);
  z.setZero();
  MATRIX eye(DOFs, DOFs);
  eye.setIdentity();

  // apply pinned vertex constraints
  for (unsigned int x = 0; x < _pinnedVertices.size(); x++)
  {
    const int pin2 = _pinnedVertices[x] * 2;
    filter(pin2, pin2) = 0;
    filter(pin2 + 1, pin2 + 1) = 0;
    z(pin2) = 0;
    z(pin2 + 1) = 0;
  }


  VECTOR cp(DOFs);
  cp.setZero();
  VECTOR2 vertexVelocity;
  MATRIX2 outerProd;
  REAL normalVelocity;
  // apply collided vertex constraints
  for (unsigned int x = 0; x < _collidedVertices.size(); x++)
  {
    const int pin2 = _collidedVertices[x] * 2;
    vertexVelocity << _velocity[pin2], _velocity[pin2 + 1];
    normalVelocity = _collidedNormals[x].transpose() * vertexVelocity;
    if (normalVelocity < 0)
    {
      //cout << "RUNNING COLLISION CODE" << endl;
      outerProd = (_collidedNormals[x] * _collidedNormals[x].transpose());
      filter(pin2, pin2) -= outerProd(0, 0);
      filter(pin2 + 1, pin2) -= outerProd(1, 0);
      filter(pin2, pin2 + 1) -= outerProd(0, 1);
      filter(pin2 + 1, pin2 + 1) -= outerProd(1, 1);
      z(pin2) = _collidedDeltas[x][0];
      z(pin2 + 1) = _collidedDeltas[x][1];
    }
  }

  MATRIX dfdx = _triangleMesh.computeStiffnessMatrix(&_hyperelastic);
  MATRIX A = _M - (_dt * _dt) * dfdx;
  MATRIX A_ = filter * A * filter.transpose() + (eye - filter);
  VECTOR b = (_dt * _dt) * (fInternal + _externalForces) + _dt * _M * _velocity;
  VECTOR b_ = filter * (b - A * z);
  VECTOR y = A_.colPivHouseholderQr().solve(b_);

  //VECTOR zeros = filter * y;
  //cout << " Sy norm should be zero! " << zeros.norm() << endl;

  VECTOR update = y + z;

  update += _position;

  VECTOR diff = (update - _position);
  _triangleMesh.setPosition(update);

  // update velocity
  _velocity = diff / _dt;

  _totalSteps++;

  return true;
}

