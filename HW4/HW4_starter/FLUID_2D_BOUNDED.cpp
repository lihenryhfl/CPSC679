///////////////////////////////////////////////////////////////////////////////
// This is an extensively reworked version of the "solver.c" file
// from Jos Stam's original "Stable Fluids" code:
//
// http://www.dgp.toronto.edu/people/stam/reality/Research/zip/CDROM_GDC03.zip
//
///////////////////////////////////////////////////////////////////////////////

#include "FLUID_2D_BOUNDED.h"
#include <cmath>
#define DIRICHLET true

///////////////////////////////////////////////////////////////////////
// Constructor / Destructor
///////////////////////////////////////////////////////////////////////
FLUID_2D_BOUNDED::FLUID_2D_BOUNDED(int xRes, int yRes, float dt) :
  FLUID_2D(xRes, yRes, dt)
{
}

///////////////////////////////////////////////////////////////////////
// perform projection using periodic boundary conditions
///////////////////////////////////////////////////////////////////////
void FLUID_2D_BOUNDED::project()
{
  // applying neumann condition
  copyBoundary(_xVelocity, true, false);
  copyBoundary(_yVelocity, false, true);

  int N = _xRes - 2;
  FIELD_2D& pressure = _xVelocityOld;
  FIELD_2D& divergence = _yVelocityOld;

  for (int y = 1; y < _yRes - 1; y++)
    for (int x = 1; x < _xRes - 1; x++)
    {
      divergence(x,y) = -0.5f * (_xVelocity(x + 1,y)  - _xVelocity(x - 1,y) +
                                 _yVelocity(x, y + 1) - _yVelocity(x, y - 1)) / N;
      pressure(x,y) = 0;
    }
  copyBoundary(pressure, true, true, 1);
  //gaussSeidel(pressure, divergence, 10);
  solvePressure(pressure, divergence, 10);
  copyBoundary(pressure, true, true, 1);

  for (int y = 1; y < _yRes - 1; y++)
    for (int x = 1; x < _xRes - 1; x++)
    {
      _xVelocity(x, y) -= 0.5f * N * (pressure(x + 1, y) - pressure(x - 1, y));
      _yVelocity(x, y) -= 0.5f * N * (pressure(x, y + 1) - pressure(x, y - 1));
    }

}

///////////////////////////////////////////////////////////////////////
// step density field using periodic boundaries
///////////////////////////////////////////////////////////////////////
void FLUID_2D_BOUNDED::stepDensity()
{
  // this adds mouse densities, they are stored in old
  addSource(_density, _densityOld);
  swapFields(_density, _densityOld);
  advect(_dt0, _density, _densityOld, _xVelocity, _yVelocity);
  fillBoundary(_density);
}

///////////////////////////////////////////////////////////////////////
// step velocity field using periodic boundaries
///////////////////////////////////////////////////////////////////////
void FLUID_2D_BOUNDED::stepVelocity()
{
  // this adds mouse forces, they are stored in old
  addSource(_xVelocity, _xVelocityOld);
  addSource(_yVelocity, _yVelocityOld);

  project();

  swapFields(_xVelocityOld, _xVelocity);
  swapFields(_yVelocityOld, _yVelocity);

  advect(_dt0, _xVelocity, _xVelocityOld, _xVelocityOld, _yVelocityOld);
  advect(_dt0, _yVelocity, _yVelocityOld, _xVelocityOld, _yVelocityOld);
}
