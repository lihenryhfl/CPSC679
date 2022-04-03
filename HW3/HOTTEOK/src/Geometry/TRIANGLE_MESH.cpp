#include "TRIANGLE_MESH.h"
#include "Hyperelastic/STVK.h"
#include <iostream>

#include <float.h>
#include <random>

using namespace std;

TRIANGLE_MESH::TRIANGLE_MESH(const std::vector<VECTOR2>& restVertices,
                             const std::vector<VECTOR3I>& triangles,
                             REAL eps) :
  _vertices(restVertices),
  _restVertices(restVertices),
  _triangles(triangles),
  _eps(eps)
{
  _DmInvs.resize(_triangles.size());
  _Fs.resize(_triangles.size());
  _restAreas.resize(_triangles.size());
  _oneRingAreas.resize(_vertices.size());
  _pFpxs.resize(_triangles.size());
  _interiorScaling = 1.5;
  _exteriorScaling = 1.5;
  _areaMultiplier = 6.;

  // precompute the Dm inverses
  for (unsigned int x = 0; x < _triangles.size(); x++)
  {
    const VECTOR3I t = _triangles[x];
    MATRIX2 Dm;
    Dm.col(0) = _restVertices[t[1]] - _restVertices[t[0]];
    Dm.col(1) = _restVertices[t[2]] - _restVertices[t[0]];
    _DmInvs[x] = Dm.inverse();
  }

  // precompute the change-of-basis
  for (unsigned int x = 0; x < _triangles.size(); x++)
    _pFpxs[x] = pFpx(_DmInvs[x]);

  // compute triangle and one-ring areas
  computeAreas();

  _R << 0., 1., -1., 0.;
  //_R << 0., -1., 1., 0.;

  // compute boundary edges
  computeBoundaryEdgesAndVertices();

  //for (unsigned int k = 0; k < _boundaryVertices.size(); k++)
    //cout << "k = " << k << ", _boundaryVertex: " << _boundaryVertices[k] << endl;

  //for (unsigned int k = 0; k < _boundaryEdges.size(); k++)
    //cout << "k = " << k << ", _boundaryEdges: " << _boundaryEdges[k] << endl;

  _staleFs = true;
}

TRIANGLE_MESH::~TRIANGLE_MESH()
{
}

///////////////////////////////////////////////////////////////////////
// find the boundary edges, and build requisite data structures
// for collision detection + resolution
///////////////////////////////////////////////////////////////////////
REAL TRIANGLE_MESH::distanceFromEdge(const VECTOR2& x, const VECTOR2& edge_start, const VECTOR2& edge_end) const
{
    VECTOR2 diff = edge_end - edge_start;
    VECTOR2 normal = _R * diff / diff.norm();
    REAL bias = -normal.transpose() * edge_start;
    REAL dist = x.transpose() * normal + bias;
    return dist;
}


bool TRIANGLE_MESH::isVertexInTriangle(const VECTOR2& xquery, const VECTOR3I& triangle, REAL eps) const
{
  // is vertex x_i inside the given triangle?
  // idea: compute three hyperplanes corresponding to the three sides
  // of the triangle. then test the sign of the new point
  //
  if (eps < 0)
    eps = _eps;

  // get points
  VECTOR2 x0, x1, x2;
  x0 = _vertices[triangle[0]];
  x1 = _vertices[triangle[1]];
  x2 = _vertices[triangle[2]];

  // build hyperplanes
  std::vector<VECTOR2> tVerts {x0, x1, x2};
  for (int i = 0; i < 3; i++) {
    VECTOR2 x_start = tVerts[i], x_end = tVerts[(i + 1) % 3];
    REAL dist = distanceFromEdge(xquery, x_start, x_end);
    // sanity check: the third vertex is always on the
    // "negative" side of the hyperplane, right?
    //assert(distanceFromEdge(tVerts[(i + 2) % 3], x_start, x_end) <= 0);

    if (dist > eps)
      return false;
  }
  return true;
}

void TRIANGLE_MESH::computeBoundaryEdgesAndVertices(REAL probeEps)
{
  // for each triangle, check each edge against every other triangle
  // to see if the edge is touching another triangle
  //
  // this check is performed by looking slightly right from
  // the midpoint of the line segment stretching from the two vertices
  // of the edge
  for (unsigned int i = 0; i < _triangles.size(); i++) {
    const VECTOR3I t = _triangles[i];
    VECTOR2 x0, x1, x2;
    x0 = _vertices[t[0]];
    x1 = _vertices[t[1]];
    x2 = _vertices[t[2]];
    std::vector<VECTOR2> tVerts {x0, x1, x2};
    for (int j = 0; j < 3; j++) {
      VECTOR2 x_start = tVerts[j], x_end = tVerts[(j + 1) % 3];
      VECTOR2 diff = x_end - x_start;
      VECTOR2 normal = _R * diff / diff.norm();
      VECTOR2 x_probe = x_start + (diff * 0.5) + (normal * probeEps);
      assert(!isVertexInTriangle(x_probe, t, 0));
      assert(isVertexInTriangle(x_probe, t, probeEps * 1.1));
      bool exterior = true;
      for (unsigned int k = 0; k < _triangles.size(); k++) {
        const VECTOR3I t_ = _triangles[k];
        if (isVertexInTriangle(x_probe, t_, 0)) {
          exterior = false;
          break;
        }
      }
      if (exterior) {
        int idx1 = t[j], idx2 = t[(j + 1) % 3];

        // add to _boundaryEdges
        VECTOR2I edge(idx1, idx2);
        REAL edgeArea = (_vertices[idx2] - _vertices[idx1]).norm();
        std::vector<int> neighbors;
        REAL extEps = _eps * _exteriorScaling;
        REAL intEps = _eps * _interiorScaling;
        for (int k = 0; k < _vertices.size(); k++) {
          if (isVertexInTriangle(_vertices[k], t, extEps)) {
            REAL dist = distanceFromEdge(_vertices[k], _vertices[idx1], _vertices[idx2]);
            if (dist > -intEps) {
              //cout << "current edge vertices: " << idx1 << " and " << idx2 << ", neighbor: " << k <<endl;
              neighbors.push_back(k);
            }
          }
        }
        _boundaryEdges.push_back(edge);
        _boundaryEdgeAreas.push_back(edgeArea);
        _boundaryEdgeNeighbors.push_back(neighbors);
        _boundaryEdgeTriangles.push_back(t);

        // add to _boundaryVertices IF it doesn't already exist
        bool alreadyAdded1 = false, alreadyAdded2 = false;
        for (unsigned int k = 0; k < _boundaryVertices.size(); k++) {
          if (_boundaryVertices[k] == idx1)
            alreadyAdded1 = true;
          if (_boundaryVertices[k] == idx2)
            alreadyAdded2 = true;
        }
        if (!alreadyAdded1) {
          _boundaryVertices.push_back(idx1);
          _boundaryVertexAreas.push_back(0.);
        }
        if (!alreadyAdded2) {
          _boundaryVertices.push_back(idx2);
          _boundaryVertexAreas.push_back(0.);
        }
      }
    }
  }

  // compute vertex areas
  for (unsigned int i = 0; i < _boundaryVertices.size(); i++) {
    for (unsigned int j = 0; j < _boundaryEdges.size(); j++) {
      int idx0 = _boundaryVertices[i], idx1 = _boundaryEdges[j][0], idx2 = _boundaryEdges[j][1];
      if ((idx0 == idx1) || (idx0 == idx2))
        _boundaryVertexAreas[i] += 0.5 * (_vertices[idx1] - _vertices[idx2]).norm();
    }
  }
}

///////////////////////////////////////////////////////////////////////
// compute some initial area quantities
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeAreas()
{
  // precompute the triangle areas
  for (unsigned int x = 0; x < _triangles.size(); x++)
  {
    // get triangle normal in R^3
    const VECTOR3I t = _triangles[x];
    VECTOR3 restPose3[3];
    for (int x = 0; x < 3; x++)
    {
      restPose3[x].setZero();
      for (int y = 0; y < 2; y++)
        restPose3[x][y] = _restVertices[t[x]][y];
    }
    const VECTOR3 edge1 = restPose3[1] - restPose3[0];
    const VECTOR3 edge2 = restPose3[2] - restPose3[0];
    _restAreas[x] = 0.5 * edge1.cross(edge2).norm();
  }

  // distribute the areas to the vertex one-rings
  for (unsigned int x = 0; x < _vertices.size(); x++)
    _oneRingAreas[x] = 0.0;

  for (unsigned int x = 0; x < _triangles.size(); x++)
  {
    // get triangle normal in R^3
    const VECTOR3I t = _triangles[x];
    _oneRingAreas[t[0]] += _restAreas[x] * (1.0 / 3.0);
    _oneRingAreas[t[1]] += _restAreas[x] * (1.0 / 3.0);
    _oneRingAreas[t[2]] += _restAreas[x] * (1.0 / 3.0);
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeFs()
{
  for (unsigned int x = 0; x < _triangles.size(); x++)
  {
    const VECTOR3I t = _triangles[x];
    MATRIX2 Ds;
    Ds.col(0) = _vertices[t[1]] - _vertices[t[0]];
    Ds.col(1) = _vertices[t[2]] - _vertices[t[0]];
    _Fs[x] = Ds * _DmInvs[x];
  }

  _staleFs = false;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR4 TRIANGLE_MESH::flatten(const MATRIX2& A)
{
  VECTOR4 result;
  result[0] = A(0,0);
  result[1] = A(1,0);
  result[2] = A(0,1);
  result[3] = A(1,1);
  return result;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MATRIX4x6 TRIANGLE_MESH::pFpx(const MATRIX2& DmInv)
{
  MATRIX4x6 result;

  for (int x = 0; x < 6; x++)
  {
    MATRIX2 column = pFpx(x, DmInv);
    result.col(x) = flatten(column);
  }
  return result;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MATRIX2 TRIANGLE_MESH::pFpx(const int index, const MATRIX2& DmInv)
{
  MATRIX2 result;
  result.setZero();

  if (index == 0)
  {
    result(0,0) = -1;
    result(0,1) = -1;
    //result << -1, -1,
    //           0, 0;
  }
  if (index == 1)
  {
    result(1,0) = -1;
    result(1,1) = -1;
    //result <<  0,  0,
    //          -1, -1;
  }
  if (index == 2)
  {
    result(0,0) = 1;
    //result <<  1, 0,
    //           0, 0;
  }
  if (index == 3)
  {
    result(1,0) = 1;
    //result <<  0, 0,
    //           1, 0;
  }
  if (index == 4)
  {
    result(0,1) = 1;
    //result <<  0, 1,
    //           0, 0;
  }
  if (index == 5)
  {
    result(1,1) = 1;
    //result <<  0, 0,
    //           0, 1;
  }

  result = result * DmInv;

  return result;
}

///////////////////////////////////////////////////////////////////////
// clamp the eigenvalues of a 6x6 to semi-positive-definite
///////////////////////////////////////////////////////////////////////
MATRIX TRIANGLE_MESH::clampEigenvalues(const MATRIX& A) const
{
  // clamp directly
  Eigen::SelfAdjointEigenSolver<MATRIX> eigensolver(A);
  const MATRIX Q = eigensolver.eigenvectors();
  VECTOR values = eigensolver.eigenvalues();
  for (int x = 0; x < values.rows(); x++)
    values[x] = (values[x] > 0.0) ? values[x] : 0.0;
  MATRIX B = Q * values.asDiagonal() * Q.transpose();

  return B;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR TRIANGLE_MESH::computeMaterialForces(const MATERIAL* material) const
{
  // you recomputed F, right?
  assert(!_staleFs);

  vector<VECTOR6> perElementForces(_triangles.size());
  for (unsigned int index = 0; index < _triangles.size(); index++)
  {
    const MATRIX2& F = _Fs[index];
    material->psi(F);
    const MATRIX2 PK1 = material->PK1(F);
    const VECTOR6 forceDensity = _pFpxs[index].transpose() * flatten(PK1);
    const VECTOR6 force = -_restAreas[index] * forceDensity;
    perElementForces[index] = force;
  }

  // scatter the forces to the global force vector, this can be parallelized
  // better where each vector entry pulls from perElementForce, but let's get
  // the slow preliminary version working first
  const int DOFs = _vertices.size() * 2;
  VECTOR forces(DOFs);
  forces.setZero();

  for (unsigned int index = 0; index < _triangles.size(); index++)
  {
    const VECTOR3I& t = _triangles[index];
    const VECTOR6& triForce = perElementForces[index];
    for (int x = 0; x < 3; x++)
    {
      assert(index < DOFs);
      assert(index + 1 < DOFs);

      unsigned int index = 2 * t[x];
      forces[index]     += triForce[2 * x];
      forces[index + 1] += triForce[2 * x + 1];
    }
  }

  return forces;
}

///////////////////////////////////////////////////////////////////////
// compute the stiffness matrix
///////////////////////////////////////////////////////////////////////
MATRIX TRIANGLE_MESH::computeStiffnessMatrix(const MATERIAL* material, bool unitTest) const
{
  // you recomputed F, right?
  assert(!_staleFs);

  vector<MATRIX6> perElementHessians(_triangles.size());
  for (unsigned int i = 0; i < _triangles.size(); i++)
  {
    const MATRIX2& F      = _Fs[i];
    const MATRIX4x6& pFpx = _pFpxs[i];
    if (unitTest) {
      const MATRIX4 hessian = material->hessian(F);
      perElementHessians[i] = -_restAreas[i] * (pFpx.transpose() * hessian) * pFpx;
    } else {
      const MATRIX4 projectedHessian = clampEigenvalues(material->hessian(F));
      perElementHessians[i] = -_restAreas[i] * (pFpx.transpose() * projectedHessian) * pFpx;
    }
  }

  const int DOFs = _vertices.size() * 2;
  MATRIX K(DOFs, DOFs);
  K.setZero();
  for (unsigned int i = 0; i < _triangles.size(); i++)
  {
    const VECTOR3I& t = _triangles[i];
    const MATRIX6& H = perElementHessians[i];
    for (int y = 0; y < 3; y++)
    {
      int yVertex = t[y];
      for (int x = 0; x < 3; x++)
      {
        int xVertex = t[x];
        K.block(2 * xVertex, 2 * yVertex, 2, 2) += H.block(2 * x, 2 * y, 2, 2);
      }
    }
  }

  return K;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR TRIANGLE_MESH::computeCollisionForces(const MATERIAL* material) const
{
  REAL eps = material->getEps();
  const int DOFs = _vertices.size() * 2;
  VECTOR forces(DOFs);
  forces.setZero();

  // check each boundary vertex against each boundary edge
  for (unsigned int i = 0; i < _boundaryVertices.size(); i++) {
    for (unsigned int j = 0; j < _boundaryEdges.size(); j++) {
      int idxs [3] = {_boundaryVertices[i], _boundaryEdges[j][0], _boundaryEdges[j][1]};

      // check if vertex is an edge neighbor
      bool neighbor = false;
      for (int k = 0; k < _boundaryEdgeNeighbors[j].size(); k++) {
        if (_boundaryEdgeNeighbors[j][k] == idxs[0]) {
          neighbor = true;
          continue;
        }
      }
      if (neighbor)
        continue;

      // if not, continue collision check
      const VECTOR3I& t = _boundaryEdgeTriangles[j];
      VECTOR2 x0 = _vertices[idxs[0]], x1 = _vertices[idxs[1]], x2 = _vertices[idxs[2]];
      if (isVertexInTriangle(x0, t)) {
        REAL dist = distanceFromEdge(x0, x1, x2);
        if (dist > (-_interiorScaling * eps)) {
          VECTOR6 x;
          x << x0, x1, x2;
          REAL area = _boundaryVertexAreas[i] + _boundaryEdgeAreas[j];
          area = area * _areaMultiplier;
          const VECTOR6 forceDensity = material->cPK1(x);
          const VECTOR6 force = -area * forceDensity;

          // fill in corresponding values
          for (int k = 0; k < 3; k++)
            forces.block(2 * idxs[k], 0, 2, 1) += force.block(2 * k, 0, 2, 1);
        }
      }
    }
  }

  return forces;
}

///////////////////////////////////////////////////////////////////////
// compute the stiffness matrix
///////////////////////////////////////////////////////////////////////
MATRIX TRIANGLE_MESH::computeCollisionHessian(const MATERIAL* material, bool unitTest) const
{
  REAL eps = material->getEps();
  const int DOFs = _vertices.size() * 2;
  MATRIX K(DOFs, DOFs);
  K.setZero();

  // check each boundary vertex against each boundary edge
  for (unsigned int i = 0; i < _boundaryVertices.size(); i++) {
    for (unsigned int j = 0; j < _boundaryEdges.size(); j++) {
      int idxs [3] = {_boundaryVertices[i], _boundaryEdges[j][0], _boundaryEdges[j][1]};

      // check if vertex is an edge neighbor
      bool neighbor = false;
      for (int k = 0; k < _boundaryEdgeNeighbors[j].size(); k++) {
        if (_boundaryEdgeNeighbors[j][k] == idxs[0]) {
          neighbor = true;
          continue;
        }
      }
      if (neighbor)
        continue;

      // if not, continue collision check
      const VECTOR3I& t = _boundaryEdgeTriangles[j];
      VECTOR2 x0 = _vertices[idxs[0]], x1 = _vertices[idxs[1]], x2 = _vertices[idxs[2]];
      if (isVertexInTriangle(x0, t)) {
        REAL dist = distanceFromEdge(x0, x1, x2);
        if (dist > (-_interiorScaling * eps)) {
          VECTOR6 x;
          x << x0, x1, x2;
          cout << idxs[0] << " in collision with edge = (" << idxs[1]
            << ", " << idxs[2] << ")" << endl;
          REAL area = _boundaryVertexAreas[i] + _boundaryEdgeAreas[j];
          area = area * _areaMultiplier;
          MATRIX hessian;
          if (unitTest)
            hessian = -area * material->cHessian(x);
          else
            hessian = -area * clampEigenvalues(material->cHessian(x));

          for (int k = 0; k < 3; k++) {
            for (int l = 0; l < 3; l++) {
              K.block(idxs[k] * 2, idxs[l] * 2, 2, 2) += hessian.block(k * 2, l * 2, 2, 2);
            }
          }
        }
      }
    }
  }

  return K;
}

///////////////////////////////////////////////////////////////////////
// displacement-based manipulations
///////////////////////////////////////////////////////////////////////
VECTOR TRIANGLE_MESH::getDisplacement() const
{
  VECTOR u(DOFs());
  for (int x = 0; x < _vertices.size(); x++)
  {
    VECTOR2 diff = _vertices[x] - _restVertices[x];
    u[2 * x]     = diff[0];
    u[2 * x + 1] = diff[1];
  }
  return u;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR TRIANGLE_MESH::getPosition() const
{
  VECTOR u(DOFs());
  for (int x = 0; x < _vertices.size(); x++)
  {
    u[2 * x]     = _vertices[x][0];
    u[2 * x + 1] = _vertices[x][1];
  }
  return u;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::setDisplacement(const VECTOR& u)
{
  assert(u.size() == DOFs());
  for (int x = 0; x < _vertices.size(); x++)
  {
    VECTOR2 diff;
    diff[0] = u[2 * x];
    diff[1] = u[2 * x + 1];
    _vertices[x] = _restVertices[x] + diff;
  }

  _staleFs = true;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::setPosition(const VECTOR& x)
{
  assert(x.size() == DOFs());
  for (int i = 0; i < _vertices.size(); i++)
  {
    VECTOR2 position;
    position[0] = x[2 * i];
    position[1] = x[2 * i + 1];
    _vertices[i] = position;
  }

  _staleFs = true;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::addDisplacement(const VECTOR& u)
{
  assert(u.size() == DOFs());
  for (int x = 0; x < _vertices.size(); x++)
  {
    VECTOR2 diff;
    diff[0] = u[2 * x];
    diff[1] = u[2 * x + 1];
    _vertices[x] += diff;
  }

  _staleFs = true;
}
