///////////////////////////////////////////////////////////////////////////////
// This is an extensively reworked version of the "solver.c" file
// from Jos Stam's original "Stable Fluids" code:
//
// http://www.dgp.toronto.edu/people/stam/reality/Research/zip/CDROM_GDC03.zip
//
///////////////////////////////////////////////////////////////////////////////

#include "FLUID_2D.h"
#include <cmath>
#define SOLVER_ACCURACY 1e-08

///////////////////////////////////////////////////////////////////////
// Constructor / Destructor
///////////////////////////////////////////////////////////////////////
FLUID_2D::FLUID_2D(int xRes, int yRes, float dt, bool WT) :
  _xRes(xRes), _yRes(yRes), _dt(dt),
  _density(xRes, yRes),
  _densityOld(xRes, yRes),
  _xVelocity(xRes, yRes),
  _xVelocityOld(xRes, yRes),
  _yVelocity(xRes, yRes),
  _yVelocityOld(xRes, yRes),
  _zVorticity(xRes, yRes),
  _vorticity(xRes, yRes),
  _residual(xRes, yRes),
  _direction(xRes, yRes),
  _q(xRes, yRes)
{
  _vorticityEps = 2000.0f;
  //_vorticityEps = 1280.0f;
  _N = _xRes - 2;
  _dt0 = _dt * _N;
  _WT = WT;
  if (_WT)
    _wTurbulence = new WTURBULENCE(_xRes, _yRes, 2);
}

///////////////////////////////////////////////////////////////////////
// swap the pointers for fields 'left' and 'right'
///////////////////////////////////////////////////////////////////////
void FLUID_2D::swapFields(FIELD_2D& left, FIELD_2D& right)
{
  float*& leftData = left.data();
  float*& rightData = right.data();

  float* temp = leftData;
  leftData = rightData;
  rightData = temp;
}

///////////////////////////////////////////////////////////////////////
// add the contents of 'source' to 'field'
///////////////////////////////////////////////////////////////////////
void FLUID_2D::addSource(FIELD_2D& field, FIELD_2D& source)
{
  for (int y = 0; y < _xRes; y++)
    for (int x = 0; x < _yRes; x++)
      field(x,y) += _dt * source(x,y);
}

///////////////////////////////////////////////////////////////////////
// user input force
///////////////////////////////////////////////////////////////////////
void FLUID_2D::addForce(int x, int y, float xForce, float yForce)
{
  _xVelocityOld(x,y) = xForce;
  _yVelocityOld(x,y) = yForce;
}

///////////////////////////////////////////////////////////////////////
// clear all fields
///////////////////////////////////////////////////////////////////////
void FLUID_2D::clear()
{
  _xVelocity.clear();
  _yVelocity.clear();
  _density.clear();

  _xVelocityOld.clear();
  _yVelocityOld.clear();
  _densityOld.clear();
}

///////////////////////////////////////////////////////////////////////
// clear old fields
///////////////////////////////////////////////////////////////////////
void FLUID_2D::clearOlds()
{
  _xVelocityOld.clear();
  _yVelocityOld.clear();
  _densityOld.clear();
}

///////////////////////////////////////////////////////////////////////
// draw the density field to GL
///////////////////////////////////////////////////////////////////////
void FLUID_2D::drawDensity()
{
  glBegin(GL_QUADS);
  int xRes, yRes, offset;
  FIELD_2D* density;
  if (_WT) {
    xRes = _wTurbulence->getXResBig();
    yRes = _wTurbulence->getYResBig();
    offset = 2;
    density = _wTurbulence->getDensityBig();
  } else {
    xRes = _xRes;
    yRes = _yRes;
    offset = 1;
    density = &_density;
  }
  float h = 1.0f / (xRes - 2 * offset);

    for (int i = offset; i < xRes - offset; i++)
    {
      float x = (i - 0.5f) * h;
      for (int j = offset; j < yRes - offset; j++)
      {
        float y = (j - 0.5f) * h;

        float tmp = (*density)(i, j);

        glColor3f(tmp, tmp, tmp);
        glVertex2f(x, y);
        glVertex2f(x + h, y);
        glVertex2f(x + h, y + h);
        glVertex2f(x, y + h);
      }
    }
  glEnd();
}

///////////////////////////////////////////////////////////////////////
// draw the velocity field to GL
///////////////////////////////////////////////////////////////////////
void FLUID_2D::drawVelocity()
{
	float h = 1.0f / (_xRes - 2);

	glColor3f(1.0f, 1.0f, 1.0f);
	glLineWidth(1.0f);

	glBegin(GL_LINES);

		for(int i = 1; i < _xRes - 1; i++)
    {
			float x = (i - 0.5f) * h;
			for(int j = 1; j < _yRes - 1; j++)
      {
				float y = (j - 0.5f) * h;
        float xDiff = _xVelocity(i,j);
        float yDiff = _yVelocity(i,j);

				glVertex2f(x, y);
				glVertex2f(x + xDiff, y + yDiff);
			}
		}

	glEnd();
}

///////////////////////////////////////////////////////////////////////
// timestep the fluid simulation with periodic boundary conditions
///////////////////////////////////////////////////////////////////////
void FLUID_2D::step()
{
  stepVelocity();
  stepDensity();
  if (_WT)
    _wTurbulence->stepTurbulenceReadable(_dt0, _xVelocity, _yVelocity);
}

///////////////////////////////////////////////////////////////////////
// add a smoke packet to the center
///////////////////////////////////////////////////////////////////////
void FLUID_2D::_addSource(int xRes, int yRes, FIELD_2D* density)
{
  float dx = 1.0 / xRes;
  float dy = 1.0 / yRes;
  for (int y = 0; y < yRes; y++)
    for (int x = 0; x < xRes; x++)
    {
      float xReal = x * dx;
      float yReal = y * dy;

      if (xReal > 0.475 && xReal < 0.525 &&
          yReal > 0.1 && yReal < 0.15)
      {
        (*density)(x, y) = 1.0;
      }
    }
}

void FLUID_2D::addSource()
{
  _addSource(_xRes, _yRes, &_density);

  if (_WT) {
    int xRes = _wTurbulence->getXResBig();
    int yRes = _wTurbulence->getYResBig();
    FIELD_2D* density = _wTurbulence->getDensityBig();
    _addSource(xRes, yRes, density);
  }
}


///////////////////////////////////////////////////////////////////////
// add some smoke here
///////////////////////////////////////////////////////////////////////
void FLUID_2D::addDensity(int x, int y, float amount)
{
  _densityOld(x, y) += amount;

  if (_WT) {
    int xRes = _wTurbulence->getXResBig();
    int yRes = _wTurbulence->getYResBig();
    FIELD_2D* density = _wTurbulence->getDensityBig();
    (*density)(2 * x, 2 * y) += amount;
    (*density)(2 * x + 1, 2 * y) += amount;
    (*density)(2 * x, 2 * y + 1) += amount;
    (*density)(2 * x + 1, 2 * y + 1) += amount;
  }
}

///////////////////////////////////////////////////////////////////////
// add a buoyancy force from the density
///////////////////////////////////////////////////////////////////////
void FLUID_2D::addBuoyancy()
{
  for (int y = 0; y < _yRes; y++)
    for (int x = 0; x < _xRes; x++)
      _yVelocity(x,y) += 0.001 * _density(x,y);
      //_yVelocity(x,y) += 0.001 * _density(x,y);
}

///////////////////////////////////////////////////////////////////////
// add a vorticity force
///////////////////////////////////////////////////////////////////////
void FLUID_2D::addVorticity()
{
  int N = _xRes - 2;
  // calculate vorticity
  for (int y = 1; y < _yRes - 1; y++) {
    for (int x = 1; x < _xRes - 1; x++) {
      float duydx = _yVelocity(x + 1, y) - _yVelocity(x - 1, y);
      float duxdy = _xVelocity(x, y + 1) - _xVelocity(x, y - 1);
      _zVorticity(x, y) = 0.5f * (duydx - duxdy) / N;
      _vorticity(x, y) = abs(_zVorticity(x, y));
    }
  }

  float dV[2];
  for (int y = 1; y < _yRes - 1; y++) {
    for (int x = 1; x < _xRes - 1; x++) {
      dV[0] = 0.5f * (_vorticity(x + 1, y) - _vorticity(x - 1, y)) / N;
      dV[1] = 0.5f * (_vorticity(x, y + 1) - _vorticity(x, y - 1)) / N;
      float magnitude = sqrt(dV[0] * dV[0] + dV[1] * dV[1]);

      if (magnitude > 0.0f) {
        dV[0] /= magnitude;
        dV[1] /= magnitude;
        _xVelocity(x, y) += (dV[1] * _zVorticity(x, y)) * _vorticityEps / N;
        _yVelocity(x, y) -= (dV[0] * _zVorticity(x, y)) * _vorticityEps / N;
      }
    }
  }
}

void FLUID_2D::fillBoundary(FIELD_2D& field, bool fillX, bool fillY, float fillWith)
{
  int xRes = field.xRes();
  int N = xRes - 2;

  // zero X boundaries (so the leftmost and rightmost columns)
  if (fillX) {
    for (int i = 1; i <= N; i++)
    {
      field(0,     i) = fillWith;
      field(N + 1, i) = fillWith;
    }
  }

  // zero Y boundaries (so the topmost and bottommost rows)
  if (fillY) {
    for (int i = 1; i <= N; i++)
    {
      field(i,     0) = fillWith;
      field(i, N + 1) = fillWith;
    }
  }

  // zero corners
  field(0, 0    ) = fillWith;
  field(0, N + 1) = fillWith;
  field(N + 1, 0  ) = fillWith;
  field(N + 1, N + 1) = fillWith;
}

//////////////////////////////////////////////////////////////////////
// solve the poisson equation with CG
//////////////////////////////////////////////////////////////////////
void FLUID_2D::solvePressure(FIELD_2D& field, FIELD_2D& b, int iterations)
{
  int x, y;

  // i = 0
  int i = 0;

  // r = b - Ax
  for (y = 1; y < _yRes - 1; y++)
    for (x = 1; x < _xRes - 1; x++)
    {
      _residual(x, y) = b(x, y) + (-4.0f * field(x, y) +
        field(x + 1, y) + field(x - 1, y) +
        field(x, y + 1) + field(x, y - 1));
    }

  // d = r
  for (y = 1; y < _yRes - 1; y++)
    for (x = 1; x < _xRes - 1; x++)
      _direction(x, y) = _residual(x, y);

  // deltaNew = transpose(r) * r
  float deltaNew = 0.0f;
  for (y = 1; y < _yRes - 1; y++)
    for (x = 1; x < _xRes - 1; x++)
      deltaNew += _residual(x, y) * _residual(x, y);

  // delta0 = deltaNew
  float delta0 = deltaNew;

  // While deltaNew > (eps^2) * delta0
  const float eps  = SOLVER_ACCURACY;
  float maxR = 2.0f * eps;
  while ((i < iterations) && (maxR > eps))
  {
    // q = Ad
    for (y = 1; y < _yRes - 1; y++)
      for (x = 1; x < _xRes - 1; x++)
      {
        _q(x, y) = -(-4.0f * _direction(x, y) +
          _direction(x + 1, y) + _direction(x - 1, y) +
          _direction(x, y + 1) + _direction(x, y - 1));
      }

    // alpha = deltaNew / (transpose(d) * q)
    float alpha = 0.0f;
    for (y = 1; y < _yRes - 1; y++)
      for (x = 1; x < _xRes - 1; x++)
        alpha += _direction(x, y) * _q(x, y);
    if (fabs(alpha) > 0.0f)
      alpha = deltaNew / alpha;

    // x = x + alpha * d
    for (y = 1; y < _yRes - 1; y++)
      for (x = 1; x < _xRes - 1; x++)
        field(x, y) += alpha * _direction(x, y);

    // r = r - alpha * q
    maxR = 0.0f;
    for (y = 1; y < _yRes - 1; y++)
      for (x = 1; x < _xRes - 1; x++)
      {
        _residual(x, y) -= alpha * _q(x, y);
        maxR = (_residual(x, y) > maxR) ? _residual(x, y) : maxR;
      }

    // deltaOld = deltaNew
    float deltaOld = deltaNew;

    // deltaNew = transpose(r) * r
    deltaNew = 0.0f;
    for (y = 1; y < _yRes - 1; y++)
      for (x = 1; x < _xRes - 1; x++)
        deltaNew += _residual(x, y) * _residual(x, y);

    // beta = deltaNew / deltaOld
    float beta = deltaNew / deltaOld;

    // d = r + beta * d
    for (y = 1; y < _yRes - 1; y++)
      for (x = 1; x < _xRes - 1; x++)
        _direction(x, y) = _residual(x, y) + beta * _direction(x, y);

    // i = i + 1
    i++;
  }
  //cout << i << " iterations converged to " << maxR << ", deltaNew " << deltaNew << endl;
  //cout << i << " iterations, deltaNew " << deltaNew << endl;
  cout << "deltaNew " << deltaNew << endl;
}


///////////////////////////////////////////////////////////////////////
// solve linear system with Gauss-Seidel iteration
///////////////////////////////////////////////////////////////////////
void FLUID_2D::gaussSeidel(FIELD_2D& pressure, FIELD_2D& divergence, int iterations)
{
	for (int k = 0; k < 10; k++)
  {
    for (int y = 1; y < _yRes - 1; y++)
      for (int x = 1; x < _xRes - 1; x++)
			  pressure(x,y) = (divergence(x,y) + pressure(x-1,y) + pressure(x+1,y) + pressure(x,y-1) + pressure(x,y+1)) * 0.25;
      // i.e.: p = d - Dp, where p is pressure, d is divergence, and D is the divergence operator
	}
}

void FLUID_2D::copyBoundary(FIELD_2D& field, bool copyX, bool copyY, int stepsAway)
{
  fillBoundary(field);
  int xRes = field.xRes();
  int N = xRes - 2;

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
void FLUID_2D::advect(float dt0, FIELD_2D& current, FIELD_2D& old, FIELD_2D& xVelocity, FIELD_2D& yVelocity)
{
  //copyBoundary(xVelocity, true, false, 1);
  //copyBoundary(yVelocity, false, true, 1);
  copyBoundary(xVelocity, false, true, 1);
  copyBoundary(yVelocity, true, false, 1);

  int xRes = current.xRes();
  int yRes = current.yRes();

  assert (xRes == old.xRes() && xRes == xVelocity.xRes() && xRes == yVelocity.xRes());
  assert (yRes == old.yRes() && yRes == xVelocity.yRes() && yRes == yVelocity.yRes());

  int N = xRes - 2;

  for (int y = 1; y < yRes - 1; y++)
    for (int x = 1; x < xRes - 1; x++)
    {
      // trace backwards through the velocity field
      float velX = -dt0 * xVelocity(x, y);
      float velY = -dt0 * yVelocity(x, y);
      float tempX = x + velX;
      float tempY = y + velY;

      if (tempX < 0.5) tempX = 0.5;
      if (tempX > xRes - 1.5) tempX = xRes - 1.5;
      if (tempY < 0.5) tempY = 0.5;
      if (tempY > yRes - 1.5) tempY = yRes - 1.5;

      if ((tempX < 1 && tempX > N) || (tempY < 1 && tempY > N)) {
        std::cout << "N: " << N << ", velX: " << velX << ", velY: " << velY << std::endl;
        std::cout << "N: " << N << ", tempX: " << tempX << ", tempY: " << tempY << std::endl;
        assert (false);
      }

      // retrieve the coordinates of the grid cells to interpolate
      int x0 = (int) tempX;
      int y0 = (int) tempY;
      int x1 = x0 + 1;
      int y1 = y0 + 1;

      if (x0 >= xRes || y0 >= yRes || x1 >= xRes || y1 >= yRes) {
        std::cout << "x: " << x << ", y: " << y << ", velX: " << velX << ", velY: " << velY << std::endl;
        std::cout << "N: " << N << ", tempX: " << tempX << ", tempY: " << tempY << std::endl;
        std::cout << "_xRes: " << xRes << ", _yRes: " << yRes <<
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

//void FLUID_2D::advect(float dt0, FIELD_2D& current, FIELD_2D& old, FIELD_2D& xVelocity, FIELD_2D& yVelocity)
//{
  //copyBoundary(xVelocity, true, false, 1);
  //copyBoundary(yVelocity, false, true, 1);

  //int xRes = current.xRes();
  //int yRes = current.yRes();

  //assert (xRes == old.xRes() && xRes == xVelocity.xRes() && xRes == yVelocity.xRes());
  //assert (yRes == old.yRes() && yRes == xVelocity.yRes() && yRes == yVelocity.yRes());

  //int N = xRes - 2;

  //for (int y = 2; y < yRes - 2; y++)
    //for (int x = 2; x < xRes - 2; x++)
    //{
      //// trace backwards through the velocity field
      //float velX = -dt0 * xVelocity(x, y);
      //float velY = -dt0 * yVelocity(x, y);
      //float tempX = x + velX;
      //float tempY = y + velY;

      //// if relevant, keep track of times
      //float t2LW, t2RW, t2TW, t2BW, clippedTime;
      //if (abs(velX) > 0.0 || abs(velY) > 0.0) {
        //// calculate time to left/right/top/bottom walls, i.e. t2(L/R/T/B)W
        //t2LW = -(x - 1) / velX;
        //t2RW = (N - x) / velX;
        //t2TW = (N - y) / velY;
        //t2BW = -(y - 1) / velY;

        //// set time to "infinity" if it is negative (i.e. we are going the opposite
        //// direction and will never reach said wall)
        //t2LW = (t2LW >= 0.0) ? t2LW : 1e8;
        //t2RW = (t2RW >= 0.0) ? t2RW : 1e8;
        //t2TW = (t2TW >= 0.0) ? t2TW : 1e8;
        //t2BW = (t2BW >= 0.0) ? t2BW : 1e8;

        //// if there is a collision...
        //if (t2LW <= 1.0 || t2RW <= 1.0 || t2TW <= 1.0 || t2BW <= 1.0) {
          ////std::cout << "there may be a collision..." << std::endl;
          //// then check collision against all four walls
          //if (t2LW >= 0 && t2LW <= t2RW && t2LW <= t2TW && t2LW <= t2BW) {
            //// if we are closest to the left wall...
            //clippedTime = t2LW;
            ////std::cout << "left " << t2LW * velX << std::endl;
          //} else if (t2RW >= 0 && t2RW <= t2LW && t2RW <= t2TW && t2RW <= t2BW) {
            //// or right wall...
            //clippedTime = t2RW;
            ////std::cout << "right " << t2RW * velX << std::endl;
          //} else if (t2TW >= 0 && t2TW <= t2RW && t2TW <= t2LW && t2TW <= t2BW) {
            //// or top wall...
            //clippedTime = t2TW;
            ////std::cout << "top " << t2TW * velY << std::endl;
          //} else if (t2BW >= 0 && t2BW <= t2RW && t2BW <= t2TW && t2BW <= t2LW) {
            //// or bottom wall...
            //clippedTime = t2BW;
            ////std::cout << "bottom " << t2BW * velY << std::endl;
          //} else
            //assert (false);

          //tempX = x + clippedTime * velX;
          //tempY = y + clippedTime * velY;
        //}
      //}

      //if ((tempX < 1 && tempX > N) || (tempY < 1 && tempY > N)) {
        //std::cout << "t2LW: " << t2LW << ", t2RW: " << t2RW << ", t2TW: " << t2TW << ", t2BW: " << t2BW << std::endl;
        //std::cout << "N: " << N << ", velX: " << velX << ", velY: " << velY << std::endl;
        //std::cout << "N: " << N << ", tempX: " << tempX << ", tempY: " << tempY << std::endl;
        //assert (false);
      //}

      //// retrieve the coordinates of the grid cells to interpolate
      //int x0 = (int) tempX;
      //int y0 = (int) tempY;
      //int x1 = x0 + 1;
      //int y1 = y0 + 1;

      //if (x0 >= xRes || y0 >= yRes || x1 >= xRes || y1 >= yRes) {
        //std::cout << "t2LW: " << t2LW << ", t2RW: " << t2RW << ", t2TW: " << t2TW << ", t2BW: " << t2BW << std::endl;
        //std::cout << "x: " << x << ", y: " << y << ", velX: " << velX << ", velY: " << velY << std::endl;
        //std::cout << "N: " << N << ", tempX: " << tempX << ", tempY: " << tempY << std::endl;
        //std::cout << "_xRes: " << _xRes << ", _yRes: " << _yRes <<
          //", x0: " << x0 << ", y0: " << y0 <<
          //", x1: " << x1 << ", y1: " << y1
          //<< std::endl;
        //assert (false);
      //}

      //// compute the interpolation weights
      //float s1 = tempX - x0;
      //float s0 = 1 - s1;
      //float t1 = tempY - y0;
      //float t0 = 1 - t1;

      //// compute the final interpolation
      //current(x,y) = s0 * (t0 * old(x0, y0) + t1 * old(x0, y1)) +
                     //s1 * (t0 * old(x1, y0) + t1 * old(x1, y1));
    //}
//}
