#include <cmath>
#include "Gilbert-Johnson-Keerthi Distance Algorithm.h"
#include "Collision Detection.h"

using namespace std;

bool DetectCollision2D(const AABB2D& rAABB2DDataA, const AABB2D& rAABB2DDataB)
{
  if (fabs(rAABB2DDataA.Center().fX - rAABB2DDataB.Center().fX) > rAABB2DDataA.Size().fWidthX + rAABB2DDataB.Size().fWidthX)
    return false;

  if (fabs(rAABB2DDataA.Center().fY - rAABB2DDataB.Center().fY) > rAABB2DDataA.Size().fWidthY + rAABB2DDataB.Size().fWidthY)
    return false;

  return true;
}

bool DetectCollision2D(const ConvexHull2D& rConvexHull2DDataA, const ConvexHull2D& rConvexHull2DDataB)
{
  return GJKIntersection2D(rConvexHull2DDataA, rConvexHull2DDataB);
}

bool DetectCollision2D(const Mesh2D& rMesh2DDataA, const Mesh2D& rMesh2DDataB)
{
  for (int i = 0; i < rMesh2DDataA.Size(); i++)
    for (int j = 0; j < rMesh2DDataB.Size(); j++)
      if (DetectCollision2D(rMesh2DDataA[i], rMesh2DDataB[j]))
        return true;

  return false;
}

bool DetectCollision2D(const vector<Vertex2D>& rvecVertexA, const vector<Vertex2D>& rvecVertexB, const PiecewisePolynomial2D& rPiecewisePolynomial2DDataA, const PiecewisePolynomial2D& rPiecewisePolynomial2DDataB, const float (*const& rpfTransformationMatrixA)[ELEMENT_COUNT_2D + 1], const float (*const& rpfTransformationMatrixB)[ELEMENT_COUNT_2D + 1])
{
  float fInversedTransformationMatrixA[ELEMENT_COUNT_2D + 1][ELEMENT_COUNT_2D + 1], fInversedTransformationMatrixB[ELEMENT_COUNT_2D + 1][ELEMENT_COUNT_2D + 1];

  for (int i = 0; i <= ELEMENT_COUNT_2D; i++)
    for (int j = 0; j <= ELEMENT_COUNT_2D; j++) {
      fInversedTransformationMatrixA[i][j] = rpfTransformationMatrixA[i][j];
      fInversedTransformationMatrixB[i][j] = rpfTransformationMatrixB[i][j];
    }

  InverseTransformation<ELEMENT_COUNT_2D + 1>(fInversedTransformationMatrixA);
  InverseTransformation<ELEMENT_COUNT_2D + 1>(fInversedTransformationMatrixB);

  for (int i = 0; i < rvecVertexA.size(); i++)
    if (rPiecewisePolynomial2DDataB.At(TransformCoordinate(fInversedTransformationMatrixB, rvecVertexA[i])) > 0.0F)
      return true;

  for (int i = 0; i < rvecVertexB.size(); i++)
    if (rPiecewisePolynomial2DDataA.At(TransformCoordinate(fInversedTransformationMatrixA, rvecVertexB[i])) > 0.0F)
      return true;

  return false;
}

bool DetectCollision3D(const AABB3D& rAABB3DDataA, const AABB3D& rAABB3DDataB)
{
  if (fabs(rAABB3DDataA.Center().fX - rAABB3DDataB.Center().fX) > rAABB3DDataA.Size().fWidthX + rAABB3DDataB.Size().fWidthX)
    return false;

  if (fabs(rAABB3DDataA.Center().fY - rAABB3DDataB.Center().fY) > rAABB3DDataA.Size().fWidthY + rAABB3DDataB.Size().fWidthY)
    return false;

  if (fabs(rAABB3DDataA.Center().fZ - rAABB3DDataB.Center().fZ) > rAABB3DDataA.Size().fWidthZ + rAABB3DDataB.Size().fWidthZ)
    return false;

  return true;
}

bool DetectCollision3D(const ConvexHull3D& rConvexHull3DDataA, const ConvexHull3D& rConvexHull3DDataB)
{
  return GJKIntersection3D(rConvexHull3DDataA, rConvexHull3DDataB);
}

bool DetectCollision3D(const Mesh3D& rMesh3DDataA, const Mesh3D& rMesh3DDataB)
{
  for (int i = 0; i < rMesh3DDataA.Size(); i++)
    for (int j = 0; j < rMesh3DDataB.Size(); j++)
      if (DetectCollision3D(rMesh3DDataA[i], rMesh3DDataB[j]))
        return true;

  return false;
}

bool DetectCollision3D(const vector<Vertex3D>& rvecVertexA, const vector<Vertex3D>& rvecVertexB, const PiecewisePolynomial3D& rPiecewisePolynomial3DDataA, const PiecewisePolynomial3D& rPiecewisePolynomial3DDataB, const float (*const& rpfTransformationMatrixA)[ELEMENT_COUNT_3D + 1], const float (*const& rpfTransformationMatrixB)[ELEMENT_COUNT_3D + 1])
{
  float fInversedTransformationMatrixA[ELEMENT_COUNT_3D + 1][ELEMENT_COUNT_3D + 1], fInversedTransformationMatrixB[ELEMENT_COUNT_3D + 1][ELEMENT_COUNT_3D + 1];

  for (int i = 0; i <= ELEMENT_COUNT_3D; i++)
    for (int j = 0; j <= ELEMENT_COUNT_3D; j++) {
      fInversedTransformationMatrixA[i][j] = rpfTransformationMatrixA[i][j];
      fInversedTransformationMatrixB[i][j] = rpfTransformationMatrixB[i][j];
    }

  InverseTransformation<ELEMENT_COUNT_3D + 1>(fInversedTransformationMatrixA);
  InverseTransformation<ELEMENT_COUNT_3D + 1>(fInversedTransformationMatrixB);

  for (int i = 0; i < rvecVertexA.size(); i++)
    if (rPiecewisePolynomial3DDataB.At(TransformCoordinate(fInversedTransformationMatrixB, rvecVertexA[i])) > 0.0F)
      return true;

  for (int i = 0; i < rvecVertexB.size(); i++)
    if (rPiecewisePolynomial3DDataA.At(TransformCoordinate(fInversedTransformationMatrixA, rvecVertexB[i])) > 0.0F)
      return true;

  return false;
}