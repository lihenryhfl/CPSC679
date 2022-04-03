///////////////////////////////////////////////////////////////////////////////
// This is an extensively reworked version of the "solver.c" file
// from Jos Stam's original "Stable Fluids" code:
//
// http://www.dgp.toronto.edu/people/stam/reality/Research/zip/CDROM_GDC03.zip
//
///////////////////////////////////////////////////////////////////////////////

#ifndef FLUID_2D_H
#define FLUID_2D_H

#ifdef USING_OSX
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif
#include "FIELD_2D.h"

class FLUID_2D {

public:
  FLUID_2D(int xRes, int yRes, float dt);
  virtual ~FLUID_2D() {};

  // timestep the fluid simulation
  void step();

  // user input density
  void addDensity(int x, int y, float amount) { _densityOld(x,y) += amount; };

  // user input force
  void addForce(int x, int y, float xForce, float yForce);

  // add a smoke packet to the center
  void addSource();

  // add a buoyancy force from the density
  void addBuoyancy();

  // add a vorticity force
  void addVorticity();

  // clear old fields
  void clearOlds();

  // clear all fields
  void clear();

  // draw the density field to GL
  void drawDensity();

  // draw the velocity field to GL
  void drawVelocity();

  void wipeBoundaries() {
    fillBoundary(_xVelocity);
    fillBoundary(_yVelocity);
    fillBoundary(_density);
  }

  void solvePressure(FIELD_2D& field, FIELD_2D& b, int iterations=10);

  void stepTurbulence(float dtOrg, FIELD_2D& xVelocity, FIELD_2D& yVelocity);

protected:
  // simulation constants
  int _xRes;
  int _yRes;
  int _iterations;
  float _dt;
  float _vorticityEps;

  // simulation variables
  FIELD_2D _density;
  FIELD_2D _densityOld;
  FIELD_2D _xVelocity;
  FIELD_2D _xVelocityOld;
  FIELD_2D _yVelocity;
  FIELD_2D _yVelocityOld;
  FIELD_2D _zVorticity;
  FIELD_2D _vorticity;
  FIELD_2D _residual;
  FIELD_2D _direction;
  FIELD_2D _q;

  // fill boundary with fillWith
  void fillBoundary(FIELD_2D& field, bool fillX=true, bool fillY=true, float fillWith=0.0);

  // swap the pointers for fields 'left' and 'right'
  void swapFields(FIELD_2D& left, FIELD_2D& right);

  // add the contents of 'source' to 'field'
  void addSource(FIELD_2D& field, FIELD_2D& source);

  // solve linear system with Gauss-Seidel iteration
  void gaussSeidel(FIELD_2D& pressure, FIELD_2D& divergence, int iterations=10);

  // advect field 'old' into 'current' using velocity field
  // 'xVelocity' and 'yVelocity' and periodic boundary conditions
  virtual void advect(FIELD_2D& current, FIELD_2D& old, FIELD_2D& xVelocity, FIELD_2D& yVelocity) = 0;

  // perform projection using periodic boundary conditions
  virtual void project() = 0;

  // step density field using periodic boundaries
  virtual void stepDensity() = 0;

  // step velocity field using periodic boundaries
  virtual void stepVelocity() = 0;
};

#endif
