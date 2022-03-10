#ifndef BACKWARD_EULER_H
#define BACKWARD_EULER_H

#include "Geometry/TRIANGLE_MESH.h"
#include "Geometry/SQUARE.h"
#include "Hyperelastic/MATERIAL.h"
#include "Timestepper/TIMESTEPPER.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
class BACKWARD_EULER : public TIMESTEPPER
{
public:
  BACKWARD_EULER(TRIANGLE_MESH& triangleMesh, MATERIAL& hyperelastic);
  virtual ~BACKWARD_EULER();

  // take a timestep
  virtual bool solve(const bool verbose) override;
};

#endif
