#include "SETTINGS.h"
#include <cmath>
#include <iostream>
#include <vector>

#include <float.h>
#include <random>

std::mt19937 gen(123);
std::uniform_real_distribution<REAL> dist(0.0, 1.0);

#include "Geometry/TRIANGLE_MESH.h"
#include "Geometry/SQUARE.h"
#include "Hyperelastic/STVK.h"
#include "Hyperelastic/SNH.h"
#include "Timestepper/TIMESTEPPER.h"
#include "Timestepper/FORWARD_EULER.h"

using namespace std;

//////////////////////////////////////////////////////////////////////////////
// do a convergence test on the normalJacobian
//////////////////////////////////////////////////////////////////////////////
bool convergenceTestNJ(const SNH* material, const VECTOR6 &x)
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING dndx for " << material->name().c_str() << endl;
  cout << "=============================================================== " << endl;

  VECTOR2 n0 = material->normal(x);
  MATRIX dndx = material->normalJacobian(x);

  double eps = 1e-4;
  int e = 0;
  double minSeen = FLT_MAX;
  while (eps > 1e-8)
  {
    MATRIX finiteDiffDndx = MATRIX(2, 6);

    // for each of the degrees of the freedom
    for (int i = 0; i < 6; i++)
    {
      VECTOR6 xnew = x;
      xnew[i] += eps;
      VECTOR2 n = material->normal(xnew);

      // store the finite difference
      VECTOR2 dn = (n - n0) / eps;
      finiteDiffDndx(0, i) = dn[0];
      finiteDiffDndx(1, i) = dn[1];
    }

    MATRIX diff = dndx - finiteDiffDndx;
    REAL diffNorm = (fabs(diff.norm() / dndx.norm())) / 12.0;
    if (diffNorm < minSeen)
      minSeen = diffNorm;
    cout << "eps: " << eps << " diff: " << diffNorm << endl;

    if (e == 4 && minSeen > 1e-6)
    {
      cout << " TEST FAILED!!!!!" << endl;
      cout << " dndx: " << endl << dndx << endl;
      cout << " finite diff: " << endl << finiteDiffDndx << endl;
      cout << " diff: " << endl << diff << endl;
      return false;
    }
    else
      eps *= 0.1;
    e++;
  }
  if (minSeen < 1e-6)
    cout << " TEST PASSED. " << endl;
  else
  {
    cout << " TEST FAILED. " << endl;
    return false;
  }
  return true;
}

//////////////////////////////////////////////////////////////////////////////
// do a convergence test on the normalHessian
//////////////////////////////////////////////////////////////////////////////
bool convergenceTestNH(const SNH* material, const VECTOR6 &x)
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING NH for " << material->name().c_str() << endl;
  cout << "=============================================================== " << endl;

  MATRIX J0 = material->normalJacobian(x);
  MATRIX H = material->normalHessian(x);

  J0 = material->reshape(J0, 12, 1);

  double eps = 1e-4;
  int e = 0;
  double minSeen = FLT_MAX;
  while (eps > 1e-8)
  {
    MATRIX finiteDiffH = MATRIX(12, 6);
    finiteDiffH.setZero();

    // for each of the degrees of the freedom
    for (int i = 0; i < 6; i++)
    {
      VECTOR6 xnew = x;
      xnew[i] += eps;
      MATRIX J = material->normalJacobian(xnew);
      J = material->reshape(J, 12, 1);

      // store the finite difference
      MATRIX dJ = (J - J0) / eps;
      for (int j = 0; j < 12; j++)
        finiteDiffH(j, i) = dJ(j, 0);
    }

    MATRIX diff = H - finiteDiffH;
    REAL diffNorm = (fabs(diff.norm() / H.norm())) / 72.0;
    if (diffNorm < minSeen)
      minSeen = diffNorm;
    cout << "eps: " << eps << " diff: " << diffNorm << endl;

    if (e == 4 && minSeen > 1e-6)
    {
      cout << " TEST FAILED!!!!!" << endl;
      cout << " H: " << endl << H << endl;
      cout << " finite diff: " << endl << finiteDiffH << endl;
      cout << " diff: " << endl << diff << endl;
      return false;
    }
    else
      eps *= 0.1;
    e++;
  }
  if (minSeen < 1e-6)
    cout << " TEST PASSED. " << endl;
  else
  {
    cout << " TEST FAILED. " << endl;
    return false;
  }
  return true;
}

//////////////////////////////////////////////////////////////////////////////
// do a convergence test on the collision Jacobian
//////////////////////////////////////////////////////////////////////////////
bool convergenceTestCJ(const SNH* material, const VECTOR6 &x)
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING Collision Jacobian for " << material->name().c_str() << endl;
  cout << "=============================================================== " << endl;

  REAL psi0 = material->cpsi(x);
  MATRIX P = material->cPK1(x); // shape (1, 6)

  double eps = 1e-4;
  int e = 0;
  double minSeen = FLT_MAX;
  while (eps > 1e-8)
  {
    MATRIX finiteDiffP = MATRIX(1, 6);
    finiteDiffP.setZero();

    // for each of the degrees of the freedom
    for (int i = 0; i < 6; i++)
    {
      VECTOR6 xnew = x;
      xnew(i) += eps;
      REAL psi = material->cpsi(xnew);

      // store the finite difference
      finiteDiffP(0, i) = (psi - psi0) / eps;
    }

    MATRIX diff = P - finiteDiffP;
    REAL diffNorm = (fabs(diff.norm() / P.norm())) / 6.0;
    if (diffNorm < minSeen)
      minSeen = diffNorm;
    cout << "eps: " << eps << " diff: " << diffNorm << endl;

    if (e == 4 && minSeen > 1e-6)
    {
      cout << " TEST FAILED!!!!!" << endl;
      cout << " P: " << endl << P << endl;
      cout << " finite diff: " << endl << finiteDiffP << endl;
      cout << " diff: " << endl << diff << endl;
      return false;
    }
    else
      eps *= 0.1;
    e++;
  }
  if (minSeen < 1e-6)
    cout << " TEST PASSED. " << endl;
  else
  {
    cout << " TEST FAILED. " << endl;
    return false;
  }
  return true;
}


//////////////////////////////////////////////////////////////////////////////
// do a convergence test on the collision Hessian
//////////////////////////////////////////////////////////////////////////////
bool convergenceTestCH(const SNH* material, const VECTOR6 &x)
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING Collision Hessian for " << material->name().c_str() << endl;
  cout << "=============================================================== " << endl;

  MATRIX H = material->cHessian(x);
  MATRIX J0 = material->cPK1(x);

  double eps = 1e-4;
  int e = 0;
  double minSeen = FLT_MAX;
  while (eps > 1e-8)
  {
    MATRIX finiteDiff = MATRIX(6, 6);
    finiteDiff.setZero();

    // for each of the degrees of the freedom
    for (int i = 0; i < 6; i++)
    {
      VECTOR6 xnew = x;
      xnew[i] += eps;

      // get the new psi
      MATRIX J = material->cPK1(xnew);

      // store the finite difference
      MATRIX diff = (J - J0) / eps;
      for (int j = 0; j < 6; j++)
        finiteDiff(j,i) = diff(0,j);
    }

    MATRIX diff = H - finiteDiff;
    REAL diffNorm = (fabs(diff.norm() / finiteDiff.norm())) / 36.0;
    if (diffNorm < minSeen)
      minSeen = diffNorm;
    cout << "eps: " << eps << " diff: " << diffNorm << endl;

    if (e == 4 && minSeen > 1e-6)
    {
      cout << " TEST FAILED!!!!!" << endl;
      cout << " H: " << endl << H << endl;
      cout << " finite diff: " << endl << finiteDiff << endl;
      cout << " diff: " << endl << diff << endl;
      return false;
    }
    eps *= 0.1;
    e++;
  }
  if (minSeen < 1e-6)
  {
    cout << " TEST PASSED. " << endl;
  }
  else
  {
    cout << " TEST FAILED. " << endl;
    return false;
  }
  return true;
}

//////////////////////////////////////////////////////////////////////////////
// do a convergence test on the PK1
//////////////////////////////////////////////////////////////////////////////
bool convergenceTestPK1(const MATERIAL* material, const MATRIX2 &F)
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING P for " << material->name().c_str() << endl;
  cout << "=============================================================== " << endl;

  REAL psi0 = material->psi(F);
  MATRIX2 P = material->PK1(F);

  double eps = 1e-4;
  int e = 0;
  double minSeen = FLT_MAX;
  while (eps > 1e-8)
  {
    MATRIX2 finiteDiffP;

    // for each of the degrees of the freedom
    for (int y = 0; y < 2; y++)
      for (int x = 0; x < 2; x++)
      {
        MATRIX2 Fnew = F;
        Fnew(x,y) += eps;
        double psi = material->psi(Fnew);

        // store the finite difference
        finiteDiffP(x, y) = (psi - psi0) / eps;
      }

    MATRIX2 diff = P - finiteDiffP;
    REAL diffNorm = (fabs(diff.norm() / P.norm())) / 4.0;
    if (diffNorm < minSeen)
      minSeen = diffNorm;
    cout << "eps: " << eps << " diff: " << diffNorm << endl;

    if (e == 4 && minSeen > 1e-6)
    {
      cout << " TEST FAILED!!!!!" << endl;
      cout << " P: " << endl << P << endl;
      cout << " finite diff: " << endl << finiteDiffP << endl;
      cout << " diff: " << endl << diff << endl;
      return false;
    }
    else
      eps *= 0.1;
    e++;
  }
  if (minSeen < 1e-6)
    cout << " TEST PASSED. " << endl;
  else
  {
    cout << " TEST FAILED. " << endl;
    return false;
  }
  return true;
}

///////////////////////////////////////////////////////////////////////
// convert a MATRIX2 to a VECTOR4 in a consistent way
///////////////////////////////////////////////////////////////////////
VECTOR4 flatten(const MATRIX2& A)
{
  VECTOR4 column;

  unsigned int index = 0;
  for (unsigned int j = 0; j < A.cols(); j++)
    for (unsigned int i = 0; i < A.rows(); i++, index++)
      column[index] = A(i,j);

  return column;
}

//////////////////////////////////////////////////////////////////////////////
// do a convergence test on Hessian
//////////////////////////////////////////////////////////////////////////////
bool convergenceTestHessian(const MATERIAL* material, const MATRIX2 &F)
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING Hessian for " << material->name().c_str() << endl;
  cout << "=============================================================== " << endl;

  MATRIX4 dPdF = material->hessian(F);
  MATRIX2 P = material->PK1(F);

  double eps = 1e-4;
  int e = 0;
  double minSeen = FLT_MAX;
  while (eps > 1e-8)
  {
    MATRIX4 finiteDiff;
    int column = 0;

    // for each of the degrees of the freedom
    for (int y = 0; y < 2; y++)
      for (int x = 0; x < 2; x++, column++)
      {
        MATRIX2 Fnew = F;
        Fnew(x,y) += eps;

        // get the new psi
        MATRIX2 Pnew = material->PK1(Fnew);

        // store the finite difference
        MATRIX2 diff = (Pnew - P) / eps;
        finiteDiff.col(column) = flatten(diff);
      }

    MATRIX4 diff = dPdF - finiteDiff;
    REAL diffNorm = (fabs(diff.norm() / P.norm())) / 16.0;
    if (diffNorm < minSeen)
      minSeen = diffNorm;
    cout << "eps: " << eps << " diff: " << diffNorm << endl;

    MATRIX4 div = finiteDiff;
    for (int y = 0; y < 4; y++)
      for (int x = 0; x < 4; x++)
        div(x,y) = div(x,y) / dPdF(x,y);

    if (e == 4 && minSeen > 1e-6)
    {
      cout << " TEST FAILED!!!!!" << endl;
      cout << " dPdF: " << endl << dPdF << endl;
      cout << " finite diff: " << endl << finiteDiff << endl;
      cout << " diff: " << endl << diff << endl;
      cout << " div: " << endl << div << endl;
      return false;
    }
    eps *= 0.1;
    e++;
  }
  if (minSeen < 1e-6)
  {
    cout << " TEST PASSED. " << endl;
  }
  else
  {
    cout << " TEST FAILED. " << endl;
    return false;
  }
  return true;
}

///////////////////////////////////////////////////////////////////////
// Let's make some random displacements
///////////////////////////////////////////////////////////////////////
VECTOR randomVector(const int DOFs, const REAL scaling)
{
  VECTOR result(DOFs);

  for (int x = 0; x < DOFs; x++)
  {
    result[x] = dist(gen);

    // randomize the sign
    if (dist(gen) < 0.5)
      result[x] *= -1.0;
  }
  result *= scaling;
  return result;
}

//////////////////////////////////////////////////////////////////////////////
// do a convergence test on Hessian
//////////////////////////////////////////////////////////////////////////////
bool convergenceTestK(const MATERIAL* material, TRIANGLE_MESH& mesh)
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING K for material " << material->name().c_str() << endl;
  cout << "=============================================================== " << endl;

  const int DOFs = mesh.DOFs();
  const VECTOR x0 = randomVector(DOFs, 1.0);
  // if things are really haywire, see if zero displacement works
  //const VECTOR x0 = randomVector(DOFs, 0.0);

  mesh.setDisplacement(x0);
  mesh.computeFs();
  const MATRIX K = mesh.computeStiffnessMatrix(material, true);
  const VECTOR f = mesh.computeMaterialForces(material);

  double eps = 1e-4;
  int e = 0;
  double minSeen = FLT_MAX;
  while (eps > 1e-8)
  {
    MATRIX finiteDiff(DOFs, DOFs);

    // for each of the degrees of the freedom
    for (int x = 0; x < DOFs; x++)
    {
      VECTOR xNew = x0;
      xNew[x] += eps;
      mesh.setDisplacement(xNew);

      mesh.computeFs();
      VECTOR fNew = mesh.computeMaterialForces(material);

      // store the finite difference
      VECTOR diff = (fNew - f) / eps;
      finiteDiff.col(x) = diff;
    }

    MATRIX diff = K - finiteDiff;
    REAL diffNorm = (fabs(diff.norm() / K.norm())) / (DOFs * DOFs);
    if (diffNorm < minSeen)
      minSeen = diffNorm;
    cout << "eps: " << eps << " diff: " << diffNorm << endl;

    if (e == 4 && minSeen > 1e-6)
    {
      cout << " TEST FAILED!!!!!" << endl;
      cout << " K: " << endl << K << endl;
      cout << " finite diff: " << endl << finiteDiff << endl;
      cout << " diff: " << endl << diff << endl;

      cout << " K norm: " << K.squaredNorm() << endl;
      cout << " finite diff  norm: " << finiteDiff.squaredNorm() << endl;
      return false;
    }
    eps *= 0.1;
    e++;
  }
  if (minSeen < 1e-6)
  {
    cout << " TEST PASSED. " << endl;
  }
  else
  {
    cout << " TEST FAILED. " << endl;
    return false;
  }
  return true;
}

///////////////////////////////////////////////////////////////////////
// Let's make some random deformation gradients
///////////////////////////////////////////////////////////////////////
MATRIX2 randomMatrix2(const REAL scaling)
{
  MATRIX2 F;
  F.setIdentity();

  for (int y = 0; y < 2; y++)
    for (int x = 0; x < 2; x++)
    {
      F(x,y) = dist(gen);

      // randomize the sign
      if (dist(gen) < 0.5)
        F(x,y) *= -1.0;
    }
  F *= scaling;
  return F;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
  MATRIX2 F = randomMatrix2(10.0);
  VECTOR6 x = randomVector(6, 1.0);
  //STVK stvk(1.0, 1.0);
  SNH material(1.0, 1.0);
  convergenceTestNJ(&material, x);
  convergenceTestNH(&material, x);
  convergenceTestCJ(&material, x);
  convergenceTestCH(&material, x);
  convergenceTestPK1(&material, F);
  convergenceTestHessian(&material, F);

  // build a square mesh
  vector<VECTOR2> nodes;
  vector<VECTOR3I> triangles;
  nodes.resize(4);
  nodes[0] = VECTOR2(0,0);
  nodes[1] = VECTOR2(1,0);
  nodes[2] = VECTOR2(0,1);
  nodes[3] = VECTOR2(1,1);
  triangles.resize(2);
  triangles[0] = VECTOR3I(0,1,2);
  triangles[1] = VECTOR3I(1,3,2);
  TRIANGLE_MESH triangleMesh(nodes, triangles);
  convergenceTestK(&material, triangleMesh);
  return 1;
}
