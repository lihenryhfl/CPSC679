///////////////////////////////////////////////////////////////////////////////
// This is an extensively reworked version of the "solver.c" file
// from Jos Stam's original "Stable Fluids" code:
//
// http://www.dgp.toronto.edu/people/stam/reality/Research/zip/CDROM_GDC03.zip
//
///////////////////////////////////////////////////////////////////////////////

#ifndef FLUID_2D_BOUNDED_H
#define FLUID_2D_BOUNDED_H

#ifdef USING_OSX
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif
#include "FIELD_2D.h"
#include "FLUID_2D.h"

class FLUID_2D_BOUNDED : public FLUID_2D {

public:
  FLUID_2D_BOUNDED(int xRes, int yRes, float dt, bool WT);
  virtual ~FLUID_2D_BOUNDED() {};

private:
  // perform projection
  virtual void project();

  // step density field
  virtual void stepDensity();

  // step velocity field
  virtual void stepVelocity();
};

#endif
