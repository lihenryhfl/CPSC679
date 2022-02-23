#include "TIMESTEPPER.h"
#include "Hyperelastic/STVK.h"
#include <float.h>
#include <iostream>

using namespace std;

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
TIMESTEPPER::TIMESTEPPER(TRIANGLE_MESH& triangleMesh, MATERIAL& hyperelastic) :
  _triangleMesh(triangleMesh),
  _hyperelastic(hyperelastic)
{
  _DOFs = _triangleMesh.DOFs();
  _forces.resize(_DOFs);
  _externalForces.resize(_DOFs);

  _forces.setZero();
  _externalForces.setZero();

  _position.resize(_DOFs);
  _velocity.resize(_DOFs);

  _position.setZero();
  _velocity.setZero();

  _name = string("TIMESTEPPER");

  //_dt = 1.0 / 77.0;  // waits a while before going kablooie
  //_dt = 1.0 / 100.0;  // waits a while before going kablooie
  //_dt = 1.0 / 60.0; // goes kablooie very quickly
  //_dt = 1.0 / 300.0;  // stable, but slow
  //_dt = 1.0 / 3000.0;  // stable, but slow
  //_dt = 1.0 / 30000.0;  // stable, but slow
  _dt = 1.0 / 30.0;
  //_dt = 1.0 / 90.0;
  //_rayleighAlpha = 0.1;
  //_rayleighBeta  = 0.1;
  //_rayleighAlpha = 0.001;
  //_rayleighBeta  = 0.001;
  _rayleighAlpha = 0.0001;
  _rayleighBeta  = 0.0001;
  _rayleighAlpha = 0.0;
  _rayleighBeta  = 0.0;

  // build the mass matrix once and for all
  buildMassMatrix();

  // cache frozen Rayleigh in case we want it later
  buildRayleighDampingMatrix();

  _totalSteps = 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// clang sure cares that this has to be here
///////////////////////////////////////////////////////////////////////////////////////////////////////
TIMESTEPPER::~TIMESTEPPER()
{
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// all vertices that start out inside this box are pinned
///////////////////////////////////////////////////////////////////////////////////////////////////////
void TIMESTEPPER::applyPinConstraints(const SQUARE& square)
{
  _pinnedVertices.clear();

  const vector<VECTOR2>& vertices = _triangleMesh.vertices();
  for (int x = 0; x < vertices.size(); x++)
  {
    if (square.inside(vertices[x]))
      _pinnedVertices.push_back(x);
  }

  cout << " Pinned " << _pinnedVertices.size() << " vertices " << endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// all vertices that start out inside this box are pinned
///////////////////////////////////////////////////////////////////////////////////////////////////////
void TIMESTEPPER::findCollidedVertices(const vector<SQUARE>& squares)
{
  _collidedVertices.clear();
  _collidedNormals.clear();
  _collidedVertices.clear();

  VECTOR2 closestPoint, normal;
  const vector<VECTOR2>& vertices = _triangleMesh.vertices();
  for (int x = 0; x < vertices.size(); x++)
  {
    for (int s = 0; s < squares.size(); s++)
    {
      if (squares[s].inside(vertices[x]))
      {
        squares[s].getClosestPoint(vertices[x], closestPoint, normal);
        _collidedVertices.push_back(x);
        _collidedDeltas.push_back(closestPoint - vertices[x]);
        _collidedNormals.push_back(normal);
      }
    }
  }

  cout << " Found " << _collidedVertices.size() << " vertex collisions " << endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// build the mass matrix based on the one-ring volumes
///////////////////////////////////////////////////////////////////////////////////////////////////////
void TIMESTEPPER::buildMassMatrix()
{
  const int DOFs = _triangleMesh.DOFs();
  _M.resize(DOFs, DOFs);
  _Minv.resize(DOFs, DOFs);

  for (int x = 0; x < _triangleMesh.vertices().size(); x++)
  {
    const REAL area = _triangleMesh.oneRingArea(x);
    const int twoX = 2 * x;
    _M(twoX, twoX)         = area;
    _M(twoX + 1, twoX + 1) = area;

    _Minv(twoX, twoX)         = 1.0 / area;
    _Minv(twoX + 1, twoX + 1) = 1.0 / area;
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// build the damping matrix based on the rest pose stiffness
///////////////////////////////////////////////////////////////////////////////////////////////////////
void TIMESTEPPER::buildRayleighDampingMatrix()
{
  // back up current state
  _temp = _triangleMesh.getDisplacement();

  // set to zero displacement
  VECTOR zero(_DOFs);
  zero.setZero();
  _triangleMesh.setDisplacement(zero);

  // get stiffness matrix at that state
  _triangleMesh.computeFs();
  //_triangleMesh.computeSVDs();
  MATRIX K = _triangleMesh.computeStiffnessMatrix(&_hyperelastic);

  // restore old state
  _triangleMesh.setDisplacement(_temp);

  // build out the Rayleigh damping
  _C = _rayleighAlpha * _M;
  _C += _rayleighBeta * K;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// add a gravity body force to the simulation
///////////////////////////////////////////////////////////////////////////////////////////////////////
void TIMESTEPPER::addGravity(const VECTOR2& bodyForce)
{
  for (int x = 0; x < _triangleMesh.vertices().size(); x++)
  {
    const REAL area = _triangleMesh.oneRingArea(x);
    const VECTOR2 scaledForce = area * bodyForce;
    _externalForces[2 * x]     += scaledForce[0];
    _externalForces[2 * x + 1] += scaledForce[1];
  }
}
