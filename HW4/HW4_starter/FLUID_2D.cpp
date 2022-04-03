///////////////////////////////////////////////////////////////////////////////
// This is an extensively reworked version of the "solver.c" file
// from Jos Stam's original "Stable Fluids" code:
//
// http://www.dgp.toronto.edu/people/stam/reality/Research/zip/CDROM_GDC03.zip
//
///////////////////////////////////////////////////////////////////////////////

#include "FLUID_2D.h"
#include <cmath>
#define SOLVER_ACCURACY 1e-06

///////////////////////////////////////////////////////////////////////
// Constructor / Destructor
///////////////////////////////////////////////////////////////////////
FLUID_2D::FLUID_2D(int xRes, int yRes, float dt) :
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
  _vorticityEps = 4000.0f;
  //_vorticityEps = 2.0f;
	//float scaling = 64.0f / _xRes;
	//scaling = (scaling < 1.0f) ? 1.0f : scaling;
	//_vorticityEps /= scaling;
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
	float h = 1.0f / (_xRes - 2);

	glBegin(GL_QUADS);
		for (int i = 1; i < _xRes - 1; i++)
    {
			float x = (i - 0.5f) * h;
			for (int j = 1; j < _yRes - 1; j++)
      {
				float y = (j - 0.5f) * h;

  		  float density = _density(i,j);

				glColor3f(density, density, density);
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
}

///////////////////////////////////////////////////////////////////////
// add a smoke packet to the center
///////////////////////////////////////////////////////////////////////
void FLUID_2D::addSource()
{
  float dx = 1.0 / _xRes;
  float dy = 1.0 / _yRes;
  for (int y = 0; y < _yRes; y++)
    for (int x = 0; x < _xRes; x++)
    {
      float xReal = x * dx;
      float yReal = y * dy;

      if (xReal > 0.475 && xReal < 0.525 &&
          yReal > 0.1 && yReal < 0.15)
      {
        _density(x,y) = 1.0;
      }
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
  int N = _xRes - 2;

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

  //for (y = 1; y < _yRes - 1; y++)
    //for (x = 1; x < _xRes - 1; x++)
      //field(x,y) = (b(x,y) + field(x-1,y) + field(x+1,y) + field(x,y-1) + field(x,y+1)) * 0.25;

  // r = b - Ax
  for (y = 1; y < _yRes - 1; y++)
    for (x = 1; x < _xRes - 1; x++)
    {
      _residual(x, y) = b(x, y) - (4.0f * field(x, y) +
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
        _q(x, y) = (4.0f * _direction(x, y) +
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
