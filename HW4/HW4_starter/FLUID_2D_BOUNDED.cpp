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
// the following functions help enforce bounded boundary conditions
///////////////////////////////////////////////////////////////////////

void FLUID_2D_BOUNDED::copyBoundary(FIELD_2D& field, bool copyX, bool copyY, int stepsAway)
{
  fillBoundary(field);
  int N = _xRes - 2;

  // copy edges
  if (copyX) {
    for (int i = 1; i <= N; i++)
    {
      field(0,     i) = field(stepsAway, i);
      field(N + 1, i) = field(N + 1 - stepsAway, i);
    }
  }
  if (copyY) {
    for (int i = 1; i <= N; i++)
    {
      field(i,     0) = field(i, stepsAway);
      field(i, N + 1) = field(i, N + 1 - stepsAway);
    }
  }

  // corners are averages of neighbors
  field(0, 0    ) = 0.5 * (field(0, 1) + field(1, 0));
  field(0, N + 1) = 0.5 * (field(0, N) + field(1, N + 1));
  field(N + 1, 0  ) = 0.5 * (field(N, 0) + field(N + 1, 1));
  field(N + 1, N + 1) = 0.5 * (field(N, N + 1) + field(N + 1, N));
}

///////////////////////////////////////////////////////////////////////
// advect field 'old' into 'current' using velocity field
// 'xVelocity' and 'yVelocity' and periodic boundary conditions
///////////////////////////////////////////////////////////////////////
void FLUID_2D_BOUNDED::advect(FIELD_2D& current, FIELD_2D& old, FIELD_2D& xVelocity, FIELD_2D& yVelocity)
{
  copyBoundary(xVelocity, true, false, 1);
  copyBoundary(yVelocity, false, true, 1);
  int N = _xRes - 2;
  float dt0 = _dt * N; // convert velocity units to grid units (N grid cells = 1 model unit)

  for (int y = 2; y < _yRes - 2; y++)
    for (int x = 2; x < _xRes - 2; x++)
    {
      // trace backwards through the velocity field
      float velX = -dt0 * xVelocity(x, y);
      float velY = -dt0 * yVelocity(x, y);
      float tempX = x + velX;
      float tempY = y + velY;

      // if relevant, keep track of times
      float t2LW, t2RW, t2TW, t2BW, clippedTime;
      if (abs(velX) > 0.0 || abs(velY) > 0.0) {
        // calculate time to left/right/top/bottom walls, i.e. t2(L/R/T/B)W
        t2LW = -(x - 1) / velX;
        t2RW = (N - x) / velX;
        t2TW = (N - y) / velY;
        t2BW = -(y - 1) / velY;

        // set time to "infinity" if it is negative (i.e. we are going the opposite
        // direction and will never reach said wall)
        t2LW = (t2LW >= 0.0) ? t2LW : 1e8;
        t2RW = (t2RW >= 0.0) ? t2RW : 1e8;
        t2TW = (t2TW >= 0.0) ? t2TW : 1e8;
        t2BW = (t2BW >= 0.0) ? t2BW : 1e8;

        // if there is a collision...
        if (t2LW <= 1.0 || t2RW <= 1.0 || t2TW <= 1.0 || t2BW <= 1.0) {
          //std::cout << "there may be a collision..." << std::endl;
          // then check collision against all four walls
          if (t2LW >= 0 && t2LW <= t2RW && t2LW <= t2TW && t2LW <= t2BW) {
            // if we are closest to the left wall...
            clippedTime = t2LW;
            //std::cout << "left " << t2LW * velX << std::endl;
          } else if (t2RW >= 0 && t2RW <= t2LW && t2RW <= t2TW && t2RW <= t2BW) {
            // or right wall...
            clippedTime = t2RW;
            //std::cout << "right " << t2RW * velX << std::endl;
          } else if (t2TW >= 0 && t2TW <= t2RW && t2TW <= t2LW && t2TW <= t2BW) {
            // or top wall...
            clippedTime = t2TW;
            //std::cout << "top " << t2TW * velY << std::endl;
          } else if (t2BW >= 0 && t2BW <= t2RW && t2BW <= t2TW && t2BW <= t2LW) {
            // or bottom wall...
            clippedTime = t2BW;
            //std::cout << "bottom " << t2BW * velY << std::endl;
          } else
            assert (false);

          tempX = x + clippedTime * velX;
          tempY = y + clippedTime * velY;
        }
      }

      if ((tempX < 1 && tempX > N) || (tempY < 1 && tempY > N)) {
        std::cout << "t2LW: " << t2LW << ", t2RW: " << t2RW << ", t2TW: " << t2TW << ", t2BW: " << t2BW << std::endl;
        std::cout << "N: " << N << ", velX: " << velX << ", velY: " << velY << std::endl;
        std::cout << "N: " << N << ", tempX: " << tempX << ", tempY: " << tempY << std::endl;
        assert (false);
      }

      // retrieve the coordinates of the grid cells to interpolate
      int x0 = (int) tempX;
      int y0 = (int) tempY;
      int x1 = x0 + 1;
      int y1 = y0 + 1;

      if (x0 >= _xRes || y0 >= _yRes || x1 >= _xRes || y1 >= _yRes) {
        std::cout << "t2LW: " << t2LW << ", t2RW: " << t2RW << ", t2TW: " << t2TW << ", t2BW: " << t2BW << std::endl;
        std::cout << "x: " << x << ", y: " << y << ", velX: " << velX << ", velY: " << velY << std::endl;
        std::cout << "N: " << N << ", tempX: " << tempX << ", tempY: " << tempY << std::endl;
        std::cout << "_xRes: " << _xRes << ", _yRes: " << _yRes <<
          ", x0: " << x0 << ", y0: " << y0 <<
          ", x1: " << x1 << ", y1: " << y1
          << std::endl;
        assert (false);
      }

      // compute the interpolation weights
		  float s1 = tempX - x0;
      float s0 = 1 - s1;
      float t1 = tempY - y0;
      float t0 = 1 - t1;

      // compute the final interpolation
      current(x,y) = s0 * (t0 * old(x0, y0) + t1 * old(x0, y1)) +
                     s1 * (t0 * old(x1, y0) + t1 * old(x1, y1));
    }
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
  gaussSeidel(pressure, divergence, 10);
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
  advect(_density, _densityOld, _xVelocity, _yVelocity);
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

  advect(_xVelocity, _xVelocityOld, _xVelocityOld, _yVelocityOld);
  advect(_yVelocity, _yVelocityOld, _xVelocityOld, _yVelocityOld);
}

///////////////////////////////////////////////////////////////////////
// solve the system using conjugate gradients
///////////////////////////////////////////////////////////////////////
void FLUID_2D_BOUNDED::conjugateGradient(FIELD_2D& pressure, FIELD_2D& divergence)
{
  // IMPLEMENT ME
}

///////////////////////////////////////////////////////////////////////
// solve the system using conjugate gradients
///////////////////////////////////////////////////////////////////////
void FLUID_2D_BOUNDED::preconditionedConjugateGradient(FIELD_2D& pressure, FIELD_2D& divergence)
{
  // IMPLEMENT ME
}
