#ifndef TRIANGLE_MESH_H
#define TRIANGLE_MESH_H

#include "SETTINGS.h"
#include "Hyperelastic/MATERIAL.h"
#include <vector>
#include <map>

// TRIANGLE vertex ordering is CLOCKWISE
class TRIANGLE_MESH
{
public:
  TRIANGLE_MESH(const std::vector<VECTOR2>& restVertices,
                const std::vector<VECTOR3I>& triangles);
  ~TRIANGLE_MESH();

  const std::vector<VECTOR3I>& triangles() { return _triangles; };
  const std::vector<VECTOR2>& vertices()   { return _vertices; };
  const int DOFs() const                   { return _vertices.size() * 2; };

  const REAL triangleArea(const int index) const { return _restAreas[index]; };
  const REAL oneRingArea(const int index) const  { return _oneRingAreas[index]; };

  // get the hyperelastic gradient, negative scaled by area
  VECTOR computeMaterialForces(const MATERIAL* material) const;

  // get the hyperelastic hessian, negative scaled by area
  MATRIX computeStiffnessMatrix(const MATERIAL* material, bool unitTest=false) const;

  // compute deformation gradients
  void computeFs();

  // displacement-based manipulations
  VECTOR getDisplacement() const;
  VECTOR getPosition() const;
  void setDisplacement(const VECTOR& u);
  void setPosition(const VECTOR& x);
  void addDisplacement(const VECTOR& u);

private:
  static VECTOR4 flatten(const MATRIX2& A);
  static MATRIX4x6 pFpx(const MATRIX2& DmInv);
  static MATRIX2 pFpx(const int index, const MATRIX2& DmInv);

  // compute some initial area quantities
  void computeAreas();

  // the geometry
  std::vector<VECTOR2>  _vertices;
  std::vector<VECTOR2>  _restVertices;
  std::vector<VECTOR3I> _triangles;
  std::vector<REAL>     _restAreas;
  std::vector<REAL>     _oneRingAreas;

  // force computation variables
  std::vector<MATRIX2>   _DmInvs;
  std::vector<MATRIX2>   _Fs;
  std::vector<MATRIX4x6> _pFpxs;

  // collision detection / resolution functions
  std::vector<int> _boundaryVertices;
  std::vector<VECTOR2I> _boundaryEdges;
  bool isVertexInTriangle(const VECTOR2& xquery, const VECTOR3I& triangle, REAL eps = 2e-2);
  void computeBoundaryEdgesAndVertices(REAL probeEps = 1e-3);

  // for convenience, the 270 deg rotation matrix
  MATRIX2 _R;

  // did you recompute F after setting the displacements or positions?
  bool _staleFs;
};

#endif
