#ifndef SQUARE_H
#define SQUARE_H

#include "KINEMATIC_SHAPE.h"

// box positions are defined using R * S * x + t,
// starting from a square centered at (0,0,0), with sides of length 1,
// so the max corner is at  (0.5, 0.5)
// and the min corner is at (-0.5, -0.5)
//
// Reminder: to make a new rotation matrix in Eigen, do:
// Eigen::Rotation2D<REAL>(M_PI / 4.0).toRotationMatrix();
class SQUARE : public KINEMATIC_SHAPE
{
public:
  SQUARE(const VECTOR2& center, const REAL& scale);
	virtual ~SQUARE();

  // get the four vertices of the square
  std::vector<VECTOR2> vertices() const;

  virtual bool inside(const VECTOR2& point) const override;
  virtual REAL distance(const VECTOR2& point) const override;

  // remember that "inside" is negative with signed distance
  virtual REAL signedDistance(const VECTOR2& point) const override;

  // get the closest point on the cube, as well as the normal at the point
  virtual void getClosestPoint(const VECTOR2& query,
                               VECTOR2& closestPointLocal,
                               VECTOR2& normalLocal) const override;

  REAL getPointInfo(const VECTOR2& query,
                    VECTOR2& closestPoint,
                    VECTOR2& normal) const;
};

#endif
