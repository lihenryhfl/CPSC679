#ifndef KINEMATIC_SHAPE_H
#define KINEMATIC_SHAPE_H

#include "SETTINGS.h"

// positions are defined using R * S * x + t,
// starting from a primitive centered at (0,0,0), with radius of 1,
//
// Reminder: to make a new rotation matrix in Eigen, do:
// Eigen::Rotation2D<REAL>(M_PI / 4.0).toRotationMatrix();
class KINEMATIC_SHAPE
{
public:
  KINEMATIC_SHAPE() {};
	~KINEMATIC_SHAPE() {};

  virtual bool inside(const VECTOR2& point) const = 0;
  virtual REAL distance(const VECTOR2& point) const = 0;

  // remember that "inside" is negative with signed distance
  virtual REAL signedDistance(const VECTOR2& point) const = 0;

  const MATRIX2& rotation() const    { return _rotation; };
  const MATRIX2& scale() const       { return _scale; };
  const VECTOR2& translation() const { return _translation; };
  MATRIX2& rotation()    { return _rotation; };
  MATRIX2& scale()       { return _scale; };
  MATRIX2& scaleInverse()       { return _scaleInverse; };
  VECTOR2& translation() { return _translation; };

  // transform a vertex from the local space to world.
  // this is used to track how a constraint has moved. If the object
  // moves, we can then see where a specific point on the object has
  // now moved
  virtual VECTOR2 localVertexToWorld(const VECTOR2& local) const
  {
    return _rotation * _scale * local + _translation;
  };

  // transform a vertex from world space to local.
  // this is used to track constraints. If we want to attach a node
  // to a specific place on the kinematic object, we need a local i
  // coordinate in the object's local space. Then later, when the i
  // object moves, we can get the relative world position by calling
  // localToWorld
  virtual VECTOR2 worldVertexToLocal(const VECTOR2& world) const
  {
    return _scaleInverse * _rotation.transpose() * (world - _translation);
  };

  // transform a normal from the local space to world.
  // The reasoning is identical to localVertexToWorld
  VECTOR2 localNormalToWorld(const VECTOR2& normal) const { return _rotation * normal; };

  // get the closest point on the object, as well as the normal at
  // the point
  virtual void getClosestPoint(const VECTOR2& query,
                               VECTOR2& closestPointLocal,
                               VECTOR2& normalLocal) const = 0;
protected:
  VECTOR2 _center;
  MATRIX2 _scale;
  MATRIX2 _scaleInverse;
  MATRIX2 _rotation;
  VECTOR2 _translation;
};

#endif
