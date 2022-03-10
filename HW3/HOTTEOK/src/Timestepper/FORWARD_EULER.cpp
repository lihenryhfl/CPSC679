#include "FORWARD_EULER.h"
#include "Hyperelastic/STVK.h"
#include <float.h>
#include <iostream>

using namespace std;

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
FORWARD_EULER::FORWARD_EULER(TRIANGLE_MESH& triangleMesh, MATERIAL& hyperelastic) :
  TIMESTEPPER(triangleMesh, hyperelastic)
{
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// clang sure cares that this has to be here
///////////////////////////////////////////////////////////////////////////////////////////////////////
FORWARD_EULER::~FORWARD_EULER()
{
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
bool FORWARD_EULER::solve(const bool verbose)
{
  cout << "========================================================= " << endl;
  cout << "  Forward Euler step: " << _totalSteps << endl;
  cout << "========================================================= " << endl;

  _position = _triangleMesh.getPosition();

  _triangleMesh.computeFs();
  VECTOR fInternal = _triangleMesh.computeMaterialForces(&_hyperelastic);
  cout << " R norm: " << fInternal.norm() << " _dt: " << _dt << endl;

  VECTOR update = (_dt * _dt) * _Minv * (fInternal + _externalForces + _C * _velocity) + _dt * _velocity;

  //const int DOFs = _triangleMesh.DOFs();
  //MATRIX filter(DOFs, DOFs);
  //filter.setIdentity();

  //for (unsigned int x = 0; x < _pinnedVertices.size(); x++)
  //{
    //const int pin2 = _pinnedVertices[x] * 2;
    //filter(pin2, pin2) = 0;
    //filter(pin2 + 1, pin2 + 1) = 0;
  //}

  //update = filter * update;
  update += _position;

  VECTOR diff = (update - _position);
  _triangleMesh.setPosition(update);

  // update velocity
  _velocity = diff / _dt;

  _totalSteps++;

  return true;
}
