#ifndef FORWARD_EULER_H
#define FORWARD_EULER_H

#include "Geometry/TRIANGLE_MESH.h"
#include "Geometry/SQUARE.h"
#include "Hyperelastic/MATERIAL.h"
#include "Timestepper/TIMESTEPPER.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
class FORWARD_EULER : public TIMESTEPPER
{
public:
  FORWARD_EULER(TRIANGLE_MESH& triangleMesh, MATERIAL& hyperelastic);
  virtual ~FORWARD_EULER(); 

  // take a timestep
  virtual bool solve(const bool verbose) override;
};

#endif
