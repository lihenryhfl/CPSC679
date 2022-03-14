#include "TRIANGLE_MESH.h"
#include "Hyperelastic/STVK.h"
#include <iostream>

#include <float.h>
#include <random>

using namespace std;

TRIANGLE_MESH::TRIANGLE_MESH(const std::vector<VECTOR2>& restVertices,
                             const std::vector<VECTOR3I>& triangles) :
  _vertices(restVertices),
  _restVertices(restVertices),
  _triangles(triangles)
{
  _DmInvs.resize(_triangles.size());
  _Fs.resize(_triangles.size());
  _restAreas.resize(_triangles.size());
  _oneRingAreas.resize(_vertices.size());
  _pFpxs.resize(_triangles.size());

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
bool TRIANGLE_MESH::isVertexInTriangle(const VECTOR2& xquery, const VECTOR3I& triangle, REAL eps)
{
  // is vertex x_i inside the given triangle?
  // idea: compute three hyperplanes corresponding to the three sides
  // of the triangle. then test the sign of the new point

  // get points
  VECTOR2 x0, x1, x2;
  x0 = _vertices[triangle[0]];
  x1 = _vertices[triangle[1]];
  x2 = _vertices[triangle[2]];

  // build hyperplanes
  std::vector<VECTOR2> tVerts {x0, x1, x2};
  for (int i = 0; i < 3; i++) {
    VECTOR2 x_start = tVerts[i], x_end = tVerts[(i + 1) % 3];
    VECTOR2 diff = x_end - x_start;
    VECTOR2 normal = _R * diff / diff.norm();
    REAL bias = -normal.transpose() * x_start;
    REAL dist = xquery.transpose() * normal + bias;

    // sanity check: the third vertex is always on the
    // "negative" side of the hyperplane, right?
    assert((tVerts[(i + 2) % 3].transpose() * normal + bias) < 0);

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
        int idx1 = t[j], idx2 = t[(j + 1) % 3], idx3 = t[(j + 2) % 3];

        // add to _boundaryEdges
        VECTOR2I edge(idx1, idx2);
        REAL edgeArea = (_vertices[idx2] - _vertices[idx1]).norm();
        _boundaryEdges.push_back(edge);
        _boundaryEdgeAreas.push_back(edgeArea);

        // add to _boundaryVertices IF it doesn't already exist
        bool alreadyAdded1 = false, alreadyAdded2 = false;
        for (unsigned int k = 0; k < _boundaryVertices.size(); k++) {
          if (_boundaryVertices[k] == idx1)
            alreadyAdded1 = true;
          if (_boundaryVertices[k] == idx2)
            alreadyAdded2 = true;
        }
        if (!alreadyAdded1) {
          REAL vertexArea1 = 0.5 * ((_vertices[idx2] - _vertices[idx1]).norm()
            + (_vertices[idx3] - _vertices[idx1]).norm());
          _boundaryVertices.push_back(idx1);
          _boundaryVertexAreas.push_back(vertexArea1);
        }
        if (!alreadyAdded2) {
          REAL vertexArea2 = 0.5 * ((_vertices[idx1] - _vertices[idx2]).norm()
            + (_vertices[idx3] - _vertices[idx2]).norm());
          _boundaryVertices.push_back(idx2);
          _boundaryVertexAreas.push_back(vertexArea2);
        }
      }
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
// clamp the eigenvalues of a 9x9 to semi-positive-definite
///////////////////////////////////////////////////////////////////////
//MATRIX4 clampEigenvalues(const MATRIX4& A)
//{
  //// clamp directly
  //Eigen::SelfAdjointEigenSolver<MATRIX4> eigensolver(A);
  //const MATRIX4 Q = eigensolver.eigenvectors();
  //VECTOR4 values = eigensolver.eigenvalues();
  //for (int x = 0; x < 4; x++)
    //values[x] = (values[x] > 0.0) ? values[x] : 0.0;
  //MATRIX4 B = Q * values.asDiagonal() * Q.transpose();

  //return B;
//}

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
        for (int b = 0; b < 2; b++)
          for (int a = 0; a < 2; a++)
          {
            assert(2 * x + a < 6);
            assert(2 * y + b < 6);
            assert(2 * xVertex + a < DOFs);
            assert(2 * yVertex + b < DOFs);
            const REAL entry = H(2 * x + a, 2 * y + b);
            K(2 * xVertex + a, 2 * yVertex + b) += entry;
          }
      }
    }
  }

  return K;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR TRIANGLE_MESH::computeCollisionForces(const MATERIAL* material) const
{
  // you recomputed F, right?
  assert(!_staleFs);

  REAL eps = material->getEps();
  const int DOFs = _vertices.size() * 2;
  VECTOR forces(DOFs);
  forces.setZero();
  for (unsigned int i = 0; i < _boundaryVertices.size(); i++) {
    for (unsigned int j = 0; j < _boundaryEdges.size(); j++) {
      VECTOR2 x0 = _vertices[_boundaryVertices[i]];
      VECTOR2 x1 = _vertices[_boundaryEdges[i][0]];
      VECTOR2 x2 = _vertices[_boundaryEdges[i][1]];
      VECTOR2 normal = _R * (x2 - x1) / (x2 - x1).norm();
      // check for collision. if so, compute forces
      if (normal.transpose() * (x0 - x1) < eps) {
        VECTOR6 x;
        x << x0, x1, x2;
        REAL area = _boundaryVertexAreas[i] + _boundaryEdgeAreas[i];
        const VECTOR6 forceDensity = material->cPK1(x);
        const VECTOR6 force = -area * forceDensity;

        // fill in corresponding values
        forces[2 * _boundaryVertices[i]] += force[0];
        forces[2 * _boundaryVertices[i] + 1] += force[1];
        forces[2 * _boundaryEdges[i][0]] += force[2];
        forces[2 * _boundaryEdges[i][0] + 1] += force[3];
        forces[2 * _boundaryEdges[i][1]] += force[4];
        forces[2 * _boundaryEdges[i][1] + 1] += force[5];
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
        for (int b = 0; b < 2; b++)
          for (int a = 0; a < 2; a++)
          {
            assert(2 * x + a < 6);
            assert(2 * y + b < 6);
            assert(2 * xVertex + a < DOFs);
            assert(2 * yVertex + b < DOFs);
            const REAL entry = H(2 * x + a, 2 * y + b);
            K(2 * xVertex + a, 2 * yVertex + b) += entry;
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
