#include "SQUARE.h"
#include <iostream>
#include <cmath>

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
  VECTOR2 closestPoint, normal;
  return abs(getPointInfo(point, closestPoint, normal));
}

///////////////////////////////////////////////////////////////////////
// signed distance to the box
///////////////////////////////////////////////////////////////////////
REAL SQUARE::signedDistance(const VECTOR2& point) const
{
  VECTOR2 closestPoint, normal;
  if (inside(point))
    return getPointInfo(point, closestPoint, normal);
  else
    return -getPointInfo(point, closestPoint, normal);
}

//////////////////////////////////////////////////////////////////////
// get the closest point on the cube, as well as the normal at the point
//////////////////////////////////////////////////////////////////////
void SQUARE::getClosestPoint(const VECTOR2& query, VECTOR2& closestPoint, VECTOR2& normal) const
{
  getPointInfo(query, closestPoint, normal);
}

REAL SQUARE::getPointInfo(const VECTOR2& query, VECTOR2& closestPoint, VECTOR2& normal) const
{
  VECTOR2 closestPointLocal, normalLocal;
  // transform back to local coordinates
  VECTOR2 transformed = worldVertexToLocal(query);
  REAL distance;

  REAL left = abs(transformed[0] + 0.5);
  REAL right = abs(transformed[0] - 0.5);
  REAL bottom = abs(transformed[1] + 0.5);
  REAL top = abs(transformed[1] - 0.5);

  // find closest side of square, then compute closestPoint and normal wrt that side
  if (left <= right && left <= bottom && left <= top) {
    closestPointLocal << -0.5, transformed[1];
    normalLocal << -1.0, 0.0;
  } else if (right <= left && right <= bottom && right <= top) {
    closestPointLocal << 0.5, transformed[1];
    normalLocal << 1.0, 0.0;
  } else if (bottom <= right && bottom <= left && bottom <= top) {
    closestPointLocal << transformed[0], -0.5;
    normalLocal << 0.0, -1.0;
  } else if (top <= right && top <= left && top <= bottom) {
    closestPointLocal << transformed[0], 0.5;
    normalLocal << 0.0, 1.0;
  } else {
    throw;
  }

  closestPoint = localVertexToWorld(closestPointLocal);
  normal = _rotation * normalLocal;

  distance = (query - closestPoint).norm();

  return distance;
}
