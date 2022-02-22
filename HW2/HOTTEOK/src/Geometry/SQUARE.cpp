#include "SQUARE.h"
#include <iostream>

using namespace std;

///////////////////////////////////////////////////////////////////////
// box positions are defined using R * S * x + t
///////////////////////////////////////////////////////////////////////
SQUARE::SQUARE(const VECTOR2& center, const REAL& scale) 
{
  _scale = MATRIX2::Identity() * scale;
  _rotation = MATRIX2::Identity();
  _translation = center;
  _scaleInverse = _scale.inverse();
}

SQUARE::~SQUARE()
{
}

///////////////////////////////////////////////////////////////////////
// is a point inside the box?
///////////////////////////////////////////////////////////////////////
bool SQUARE::inside(const VECTOR2& point) const 
{
  // transform back to local coordinates
  VECTOR2 transformed = worldVertexToLocal(point);

  if (transformed[0] <= 0.5 && transformed[0] >= -0.5 &&
      transformed[1] <= 0.5 && transformed[1] >= -0.5)
    return true;

  return false;
}

///////////////////////////////////////////////////////////////////////
// get the four vertices of the square
///////////////////////////////////////////////////////////////////////
vector<VECTOR2> SQUARE::vertices() const
{
  vector<VECTOR2> vs(4);
  vs[0] = VECTOR2(-0.5,-0.5);
  vs[1] = VECTOR2(0.5,-0.5);
  vs[2] = VECTOR2(0.5,0.5);
  vs[3] = VECTOR2(-0.5,0.5);

  for (int x = 0; x < 4; x++)
    vs[x] = _rotation * _scale * vs[x] + _translation;

  return vs;
}

///////////////////////////////////////////////////////////////////////
// distance to the box
///////////////////////////////////////////////////////////////////////
REAL SQUARE::distance(const VECTOR2& point) const 
{
  // this isn't implemented yet
  return 0;
}

///////////////////////////////////////////////////////////////////////
// signed distance to the box
///////////////////////////////////////////////////////////////////////
REAL SQUARE::signedDistance(const VECTOR2& point) const
{
  // this isn't implemented yet
  return 0;
}

//////////////////////////////////////////////////////////////////////
// get the closest point on the cube, as well as the normal at the point
//////////////////////////////////////////////////////////////////////
void SQUARE::getClosestPoint(const VECTOR2& query, VECTOR2& closestPoint, VECTOR2& normal) const
{
}
