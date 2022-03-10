#ifndef TIMESTEPPER_H
#define TIMESTEPPER_H

#include "Geometry/TRIANGLE_MESH.h"
#include "Geometry/SQUARE.h"
#include "Hyperelastic/MATERIAL.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
class TIMESTEPPER
{
public:
  TIMESTEPPER(TRIANGLE_MESH& triangleMesh, MATERIAL& hyperelastic);
  virtual ~TIMESTEPPER();

  VECTOR& externalForces()             { return _externalForces; };
  const VECTOR& externalForces() const { return _externalForces; };
  const REAL dt() const                { return _dt; };
  REAL& dt()                           { return _dt; };

  const TRIANGLE_MESH& triangleMesh() const       { return _triangleMesh; };
  const VECTOR position() const         { return _position; };
  const VECTOR velocity() const         { return _velocity; };
  VECTOR& position()                    { return _position; };
  VECTOR& velocity()                    { return _velocity; };

  // all vertices that start out inside this box are pinned
  void applyPinConstraints(const SQUARE& square);

  // checks for collisions given squares
  void findCollidedVertices(const std::vector<SQUARE>& squares);

  // take a timestep
  virtual bool solve(const bool verbose) = 0;

  // add a gravity body force to the simulation
  void addGravity(const VECTOR2& bodyForce);
  void clearExternalForces()            { _externalForces.setZero(); };

protected:
  // build the mass matrix based on the one-ring volumes
  void buildMassMatrix();

  // build the damping matrix based on the rest pose stiffness
  void buildRayleighDampingMatrix();

  TRIANGLE_MESH& _triangleMesh;
  MATERIAL& _hyperelastic;

  int _DOFs;
  VECTOR _forces;
  VECTOR _externalForces;

  // variables to solve for
  VECTOR _position;
  VECTOR _velocity;
  VECTOR _temp;

  // timestep
  REAL _dt;
  REAL _rayleighAlpha;
  REAL _rayleighBeta;

  // integrator name
  std::string _name;

  // mass matrix
  MATRIX _M;

  // cached inverse mass matrix
  MATRIX _Minv;

  // Rayleigh damping matrix
  MATRIX _C;

  // how many steps has the integrator taken?
  int _totalSteps;

  // which vertices have been pinned by a box?
  std::vector<int> _pinnedVertices;

  // which vertices are colliding with boxes?
  std::vector<int> _collidedVertices;

  // and how far are they from the boxes?
  std::vector<REAL> _signedDistances;

  // which way do we go, to project to the closest exterior point on the square?
  std::vector<VECTOR2> _collidedDeltas;

  // what are their normals?
  std::vector<VECTOR2> _collidedNormals;

  // what are their closestPoints?
  std::vector<VECTOR2> _collidedClosestPoints;

  //// which way do we go, to project to the closest exterior point on the square?
  //std::vector<REAL> _collidedDeltas;

  //// what are their normals?
  //std::vector<REAL> _collidedNormals;

  //// what are their closestPoints?
  //std::vector<REAL> _collidedClosestPoints;
};

#endif
