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


  cout << " Applying " << _collidedVertices.size() << " vertex collisions " << endl;
  VECTOR cp(DOFs);
  cp.setZero();
  VECTOR2 vertexVelocity, normal, closestPoint, delta;
  MATRIX2 outerProd;
  REAL normalVelocity;
  // apply collided vertex constraints
  for (unsigned int x = 0; x < _collidedVertices.size(); x++)
  {
    const int pin2 = _collidedVertices[x] * 2;
    vertexVelocity << _velocity[pin2], _velocity[pin2 + 1];
    normalVelocity = _collidedNormals[x].transpose() * vertexVelocity;
    //normal = VECTOR2(_collidedNormals[x * 2], _collidedNormals[x * 2 + 1]);
    //normalVelocity = normal.transpose() * vertexVelocity;
    //if (_signedDistances[x] < -0.05 ||  normalVelocity < 0)
    if (normalVelocity < 0)
    //if (true)
    {
      cout << "RUNNING COLLISION CODE" << endl;
      outerProd = (_collidedNormals[x] * _collidedNormals[x].transpose());
      filter(pin2, pin2) -= outerProd(0, 0);
      filter(pin2 + 1, pin2) -= outerProd(1, 0);
      filter(pin2, pin2 + 1) -= outerProd(0, 1);
      filter(pin2 + 1, pin2 + 1) -= outerProd(1, 1);
      z(pin2) = _collidedDeltas[x][0];
      z(pin2 + 1) = _collidedDeltas[x][1];
      //cout << " integrator delta norm " << _collidedDeltas[x].norm() << endl;
      //cout << " integrator vertices + delta" << _position(pin2) + z(pin2) << " " << _position(pin2 + 1) + z(pin2 + 1) << endl;

      cp(pin2) = _collidedClosestPoints[x][0];
      cp(pin2 + 1) = _collidedClosestPoints[x][1];
      //delta = VECTOR2(_collidedDeltas[x * 2], _collidedDeltas[x * 2 + 1]);
      //closestPoint = VECTOR2(_collidedClosestPoints[x * 2], _collidedClosestPoints[x * 2 + 1]);
      //outerProd = (normal * normal.transpose());
      //filter(pin2, pin2) -= outerProd(0, 0);
      //filter(pin2 + 1, pin2) -= outerProd(1, 0);
      //filter(pin2, pin2 + 1) -= outerProd(0, 1);
      //filter(pin2 + 1, pin2 + 1) -= outerProd(1, 1);
      //z(pin2) = delta[0];
      //z(pin2 + 1) = delta[1];
      //cout << " integrator delta norm " << delta.norm() << endl;
      //cout << " integrator delta " << delta[0] << " " << delta[1] << endl;
      //cout << " integrator vertices + delta " << _position(pin2) + z(pin2) << " " << _position(pin2 + 1) + z(pin2 + 1) << endl;

      //cp(pin2) = closestPoint[0];
      //cp(pin2 + 1) = closestPoint[1];
    }
  }

  //cout << filter << endl;
  //cout << z << endl;
  //cout << _position << endl;
  //cout << "------------------------" << endl;
  //cout << _position + z << endl;
  //cout << "------------------------" << endl;
  //cout << cp << endl;
  //cout << _position << endl;

  MATRIX dfdx = _triangleMesh.computeStiffnessMatrix(&_hyperelastic);
  MATRIX A = _M - (_dt * _dt) * dfdx;
  MATRIX A_ = filter * A * filter.transpose() + (eye - filter);
  VECTOR b = (_dt * _dt) * (fInternal + _externalForces) + _dt * _M * _velocity;
  VECTOR b_ = filter * (b - A * z);
  VECTOR y = A_.colPivHouseholderQr().solve(b_);
  //VECTOR y = A.inverse() * b_;

  VECTOR zeros = filter * y;
  cout << " Sy norm should be zero! " << zeros.norm() << endl;
  //for (unsigned int x = 0; x < _collidedVertices.size(); x++)
  //{
    //const int pin2 = _collidedVertices[x] * 2;
    //if (true)
    //{
      //b_(pin2) = _collidedDeltas[x][0];
      //b_(pin2 + 1) = _collidedDeltas[x][1];
    //}
  //}

  VECTOR update = y + z;

  update += _position;

  VECTOR diff = (update - _position);
  _triangleMesh.setPosition(update);

  // update velocity
  _velocity = diff / _dt;

  _totalSteps++;

  return true;
}

