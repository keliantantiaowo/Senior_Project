#define _USE_MATH_DEFINES
#include <random>
#undef _USE_MATH_DEFINES
#include "Geometric Operation.h"

using namespace std;

float Norm(const Vertex2D& rVertex2DVector)
{
  return sqrt(pow(rVertex2DVector.fX, 2.0F) + pow(rVertex2DVector.fY, 2.0F));
}

float Norm(const Vertex3D& rVertex3DVector)
{
  return sqrt(pow(rVertex3DVector.fX, 2.0F) + pow(rVertex3DVector.fY, 2.0F) + pow(rVertex3DVector.fZ, 2.0F));
}

float DotProduct(const Vertex2D& rVertex2DVectorA, const Vertex2D& rVertex2DVectorB)
{
  return rVertex2DVectorA.fX * rVertex2DVectorB.fX + rVertex2DVectorA.fY * rVertex2DVectorB.fY;
}

float DotProduct(const Vertex3D& rVertex3DVectorA, const Vertex3D& rVertex3DVectorB)
{
  return rVertex3DVectorA.fX * rVertex3DVectorB.fX + rVertex3DVectorA.fY * rVertex3DVectorB.fY + rVertex3DVectorA.fZ * rVertex3DVectorB.fZ;
}

Vertex2D CrossProduct(const Vertex2D& rVertex2DVectorA, const Vertex2D& rVertex2DVectorB)
{
  throw length_error("外積演算に使われるベクトルの次元数が不正です。");
}

Vertex3D CrossProduct(const Vertex3D& rVertex3DVectorA, const Vertex3D& rVertex3DVectorB)
{
  return Vertex3D{rVertex3DVectorA.fY * rVertex3DVectorB.fZ - rVertex3DVectorA.fZ * rVertex3DVectorB.fY, rVertex3DVectorA.fZ * rVertex3DVectorB.fX - rVertex3DVectorA.fX * rVertex3DVectorB.fZ, rVertex3DVectorA.fX * rVertex3DVectorB.fY - rVertex3DVectorA.fY * rVertex3DVectorB.fX};
}

Vertex2D Normalize(const Vertex2D& rVertex2DVector)
{
  float fNorm = Norm(rVertex2DVector);

  return Vertex2D{rVertex2DVector.fX / fNorm, rVertex2DVector.fY / fNorm};
}

Vertex3D Normalize(const Vertex3D& rVertex3DVector)
{
  float fNorm = Norm(rVertex3DVector);

  return Vertex3D{rVertex3DVector.fX / fNorm, rVertex3DVector.fY / fNorm, rVertex3DVector.fZ / fNorm};
}

Vertex2D TransformCoordinate(const float (*const& rpfTransformationMatrix)[ELEMENT_COUNT_2D + 1], const Vertex2D& rVertex2DData)
{
  float Vertex2D::*pfVertex2DMember[ELEMENT_COUNT_2D] = {&Vertex2D::fX, &Vertex2D::fY};
  Vertex2D Vertex2DResult{0.0F, 0.0F};

  for (int i = 0; i < ELEMENT_COUNT_2D; i++) {
    for (int j = 0; j < ELEMENT_COUNT_2D; j++)
      Vertex2DResult.*pfVertex2DMember[i] += rpfTransformationMatrix[i][j] * rVertex2DData.*pfVertex2DMember[j];

    Vertex2DResult.*pfVertex2DMember[i] += rpfTransformationMatrix[i][ELEMENT_COUNT_2D];
  }

  return Vertex2DResult;
}

Vertex3D TransformCoordinate(const float (*const& rpfTransformationMatrix)[ELEMENT_COUNT_3D + 1], const Vertex3D& rVertex3DData)
{
  float Vertex3D::*pfVertex3DMember[ELEMENT_COUNT_3D] = {&Vertex3D::fX, &Vertex3D::fY, &Vertex3D::fZ};
  Vertex3D Vertex3DResult{0.0F, 0.0F, 0.0F};

  for (int i = 0; i < ELEMENT_COUNT_3D; i++) {
    for (int j = 0; j < ELEMENT_COUNT_3D; j++)
      Vertex3DResult.*pfVertex3DMember[i] += rpfTransformationMatrix[i][j] * rVertex3DData.*pfVertex3DMember[j];

    Vertex3DResult.*pfVertex3DMember[i] += rpfTransformationMatrix[i][ELEMENT_COUNT_3D];
  }

  return Vertex3DResult;
}

void TransformCoordinate(const float (*const& rpfTransformationMatrix)[ELEMENT_COUNT_2D + 1], const vector<Vertex2D>& rvecVertex, vector<Vertex2D>& rvecTransformedVertex)
{
  rvecTransformedVertex.resize(rvecVertex.size());

  for (int i = 0; i < rvecVertex.size(); i++)
    rvecTransformedVertex[i] = TransformCoordinate(rpfTransformationMatrix, rvecVertex[i]);
}

void TransformCoordinate(const float (*const& rpfTransformationMatrix)[ELEMENT_COUNT_3D + 1], const vector<Vertex3D>& rvecVertex, vector<Vertex3D>& rvecTransformedVertex)
{
  rvecTransformedVertex.resize(rvecVertex.size());

  for (int i = 0; i < rvecVertex.size(); i++)
    rvecTransformedVertex[i] = TransformCoordinate(rpfTransformationMatrix, rvecVertex[i]);
}

void RandomizeTransformation(float (*const& rpfTransformationMatrix)[ELEMENT_COUNT_2D + 1], const Vertex2D& rVertex2DMinimumCorner, const Vertex2D& rVertex2DMaximumCorner)
{
  mt19937 othRandomNumberGenerator((random_device())());
  uniform_real_distribution<float> othDistributor;
  float (Vertex2D::*pfVertex2DMember[ELEMENT_COUNT_2D]) = {&Vertex2D::fX, &Vertex2D::fY};
  float (*pfnFunctionMatrix[ELEMENT_COUNT_2D][ELEMENT_COUNT_2D])(float fIdentifier) = {
    {cos, sin},
    {sin, cos}
  };
  Sign SignMatrix[ELEMENT_COUNT_2D][ELEMENT_COUNT_2D] = {
    {Plus, Minus},
    {Plus, Plus}
  };
  float fTheta = (othDistributor.param(uniform_real_distribution<float>::param_type(0.0F, M_PI * 2.0F)), othDistributor(othRandomNumberGenerator));

  for (int i = 0; i < ELEMENT_COUNT_2D; i++) {
    float fFactor = (othDistributor.param(uniform_real_distribution<float>::param_type(0.0F, nextafterf(1.0F, numeric_limits<float>::max()))), othDistributor(othRandomNumberGenerator));

    for (int j = 0; j < ELEMENT_COUNT_2D; j++)
      rpfTransformationMatrix[i][j] = pfnFunctionMatrix[i][j](fTheta) * (SignMatrix[i][j] == Plus ? 1.0F : SignMatrix[i][j] == Minus ? -1.0F : 0.0F);

    rpfTransformationMatrix[i][ELEMENT_COUNT_2D] = rVertex2DMinimumCorner.*pfVertex2DMember[i] * (1.0F - fFactor) + rVertex2DMaximumCorner.*pfVertex2DMember[i] * fFactor;
    rpfTransformationMatrix[ELEMENT_COUNT_2D][i] = 0.0F;
  }

  rpfTransformationMatrix[ELEMENT_COUNT_2D][ELEMENT_COUNT_2D] = 1.0F;
}

void RandomizeTransformation(float (*const& rpfTransformationMatrix)[ELEMENT_COUNT_3D + 1], const Vertex3D& rVertex3DMinimumCorner, const Vertex3D& rVertex3DMaximumCorner)
{
  mt19937 othRandomNumberGenerator((random_device())());
  uniform_real_distribution<float> othDistributor;
  float (Vertex3D::*pfVertex3DMember[ELEMENT_COUNT_3D]) = {&Vertex3D::fX, &Vertex3D::fY, &Vertex3D::fZ};
  float (*pfnFunctionMatrix[ELEMENT_COUNT_3D][ELEMENT_COUNT_3D][ELEMENT_COUNT_3D])(float fIdentifier) = {
    {
      {nullptr, nullptr, nullptr},
      {nullptr, cos, sin},
      {nullptr, sin, cos}
    },
    {
      {cos, nullptr, sin},
      {nullptr, nullptr, nullptr},
      {sin, nullptr, cos}
    },
    {
      {cos, sin, nullptr},
      {sin, cos, nullptr},
      {nullptr, nullptr, nullptr}
    }
  };
  Sign SignMatrix[ELEMENT_COUNT_3D][ELEMENT_COUNT_3D] = {
    {Plus, Minus, Plus},
    {Plus, Plus, Minus},
    {Minus, Plus, Plus}
  };
  float fRotationMatrix[ELEMENT_COUNT_3D + 1][ELEMENT_COUNT_3D][ELEMENT_COUNT_3D];
  float fComposedRotationMatrix[ELEMENT_COUNT_3D][ELEMENT_COUNT_3D];
  float fTheta = (othDistributor.param(uniform_real_distribution<float>::param_type(0.0F, M_PI * 2.0F)), othDistributor(othRandomNumberGenerator)), fPhi = (othDistributor.param(uniform_real_distribution<float>::param_type(0.0F, M_PI * 2.0F)), othDistributor(othRandomNumberGenerator)), fPsi = (othDistributor.param(uniform_real_distribution<float>::param_type(0.0F, M_PI * 2.0F)), othDistributor(othRandomNumberGenerator));

  for (int i = 0; i < ELEMENT_COUNT_3D; i++)
    for (int j = 0; j < ELEMENT_COUNT_3D; j++) {
      fRotationMatrix[INDICATOR_X][i][j] = (pfnFunctionMatrix[INDICATOR_X][i][j] == nullptr ? 0.0F : pfnFunctionMatrix[INDICATOR_X][i][j](fTheta)) * (SignMatrix[i][j] == Plus ? 1.0F : SignMatrix[i][j] == Minus ? -1.0F : 0.0F);
      fRotationMatrix[INDICATOR_Y][i][j] = (pfnFunctionMatrix[INDICATOR_Y][i][j] == nullptr ? 0.0F : pfnFunctionMatrix[INDICATOR_Y][i][j](fPhi)) * (SignMatrix[i][j] == Plus ? 1.0F : SignMatrix[i][j] == Minus ? -1.0F : 0.0F);
      fRotationMatrix[INDICATOR_Z][i][j] = (pfnFunctionMatrix[INDICATOR_Z][i][j] == nullptr ? 0.0F : pfnFunctionMatrix[INDICATOR_Z][i][j](fPsi)) * (SignMatrix[i][j] == Plus ? 1.0F : SignMatrix[i][j] == Minus ? -1.0F : 0.0F);
    }

  fRotationMatrix[INDICATOR_X][INDICATOR_X][INDICATOR_X] = fRotationMatrix[INDICATOR_Y][INDICATOR_Y][INDICATOR_Y] = fRotationMatrix[INDICATOR_Z][INDICATOR_Z][INDICATOR_Z] = 1.0F;

  MultiplyMatrix<float, ELEMENT_COUNT_3D, ELEMENT_COUNT_3D, ELEMENT_COUNT_3D>(fRotationMatrix[INDICATOR_X], fRotationMatrix[INDICATOR_Y], fComposedRotationMatrix);
  MultiplyMatrix<float, ELEMENT_COUNT_3D, ELEMENT_COUNT_3D, ELEMENT_COUNT_3D>(fComposedRotationMatrix, fRotationMatrix[INDICATOR_Z], fRotationMatrix[ELEMENT_COUNT_3D]);

  for (int i = 0; i < ELEMENT_COUNT_3D; i++) {
    float fFactor = (othDistributor.param(uniform_real_distribution<float>::param_type(0.0F, nextafterf(1.0F, numeric_limits<float>::max()))), othDistributor(othRandomNumberGenerator));

    for (int j = 0; j < ELEMENT_COUNT_3D; j++)
      rpfTransformationMatrix[i][j] = fRotationMatrix[ELEMENT_COUNT_3D][i][j];

    rpfTransformationMatrix[i][ELEMENT_COUNT_3D] = rVertex3DMinimumCorner.*pfVertex3DMember[i] * (1.0F - fFactor) + rVertex3DMaximumCorner.*pfVertex3DMember[i] * fFactor;
    rpfTransformationMatrix[ELEMENT_COUNT_3D][i] = 0.0F;
  }

  rpfTransformationMatrix[ELEMENT_COUNT_3D][ELEMENT_COUNT_3D] = 1.0F;
}

float Distance(const Vertex2D& rVertex2DDataA, const Vertex2D& rVertex2DDataB)
{
  return Norm(Vertex2D{rVertex2DDataA.fX - rVertex2DDataB.fX, rVertex2DDataA.fY - rVertex2DDataB.fY});
}

float Distance(const Vertex2D& rVertex2DData, const Vertex2D& rVertex2DEndpointA, const Vertex2D& rVertex2DEndpointB)
{
  Vertex2D Vertex2DEdgeNormal = EdgeNormal(Vertex2D{rVertex2DEndpointB.fX - rVertex2DEndpointA.fX, rVertex2DEndpointB.fY - rVertex2DEndpointA.fY});

  return DotProduct(Vertex2D{rVertex2DData.fX - rVertex2DEndpointA.fX, rVertex2DData.fY - rVertex2DEndpointA.fY}, Normalize(Vertex2DEdgeNormal));
}

float Distance(const vector<Vertex2D>& rvecVertex, const vector<Vertex2D>& rvecVertexNormal, const Vertex2D& rVertex2DData, const int& riIndex)
{
  Vertex2D Vertex2DEdge{rvecVertex[(riIndex + 1) % rvecVertex.size()].fX - rvecVertex[riIndex].fX, rvecVertex[(riIndex + 1) % rvecVertex.size()].fY - rvecVertex[riIndex].fY};
  float fFactor = DotProduct(Vertex2D{rVertex2DData.fX - rvecVertex[riIndex].fX, rVertex2DData.fY - rvecVertex[riIndex].fY}, Vertex2DEdge) / DotProduct(Vertex2DEdge, Vertex2DEdge);

  if (fFactor <= 0.0F)
    return copysign(Distance(rvecVertex[riIndex], rVertex2DData), DotProduct(Vertex2D{rVertex2DData.fX - rvecVertex[riIndex].fX, rVertex2DData.fY - rvecVertex[riIndex].fY}, rvecVertexNormal[riIndex]));
  else if (fFactor >= 1.0F)
    return copysign(Distance(rvecVertex[(riIndex + 1) % rvecVertex.size()], rVertex2DData), DotProduct(Vertex2D{rVertex2DData.fX - rvecVertex[(riIndex + 1) % rvecVertex.size()].fX, rVertex2DData.fY - rvecVertex[(riIndex + 1) % rvecVertex.size()].fY}, rvecVertexNormal[(riIndex + 1) % rvecVertexNormal.size()]));
  else {
    Vertex2D Vertex2DClosestPoint{rvecVertex[riIndex].fX + Vertex2DEdge.fX * fFactor, rvecVertex[riIndex].fY + Vertex2DEdge.fY * fFactor};

    return copysign(Distance(rVertex2DData, Vertex2DClosestPoint), DotProduct(Vertex2D{rVertex2DData.fX - Vertex2DClosestPoint.fX, rVertex2DData.fY - Vertex2DClosestPoint.fY}, EdgeNormal(Vertex2DEdge)));
  }
}

float Distance(const Vertex3D& rVertex3DDataA, const Vertex3D& rVertex3DDataB)
{
  return Norm(Vertex3D{rVertex3DDataA.fX - rVertex3DDataB.fX, rVertex3DDataA.fY - rVertex3DDataB.fY, rVertex3DDataA.fZ - rVertex3DDataB.fZ});
}

float Distance(const Vertex3D& rVertex3DData, const Vertex3D& rVertex3DEndpointA, const Vertex3D& rVertex3DEndpointB)
{
  Vertex3D Vertex3DVectorA{rVertex3DData.fX - rVertex3DEndpointA.fX, rVertex3DData.fY - rVertex3DEndpointA.fY, rVertex3DData.fZ - rVertex3DEndpointA.fZ}, Vertex3DVectorB{rVertex3DEndpointB.fX - rVertex3DEndpointA.fX, rVertex3DEndpointB.fY - rVertex3DEndpointA.fY, rVertex3DEndpointB.fZ - rVertex3DEndpointA.fZ};

  return sqrt(DotProduct(Vertex3DVectorA, Vertex3DVectorA) - pow(DotProduct(Vertex3DVectorA, Vertex3DVectorB), 2.0F) / DotProduct(Vertex3DVectorB, Vertex3DVectorB));
}

float Distance(const vector<Vertex3D>& rvecVertex, const vector<Vertex3D>& rvecVertexNormal, const vector<Vertex3D>& rvecSurfaceNormal, const vector<vector<unsigned>>& rvecTriangleVertexIndex, const vector<vector<unsigned>>& rvecTriangleNormalIndex, const Vertex3D& rVertex3DData, const int& riIndex)
{
  Vertex3D Vertex3DEdge[TRIANGLE_ELEMENT_COUNT];
  float fFactor[TRIANGLE_ELEMENT_COUNT];

  for (int i = 0; i < TRIANGLE_ELEMENT_COUNT; i++) {
    Vertex3DEdge[i] = Vertex3D{rvecVertex[rvecTriangleVertexIndex[riIndex][(i + 1) % TRIANGLE_ELEMENT_COUNT]].fX - rvecVertex[rvecTriangleVertexIndex[riIndex][i]].fX, rvecVertex[rvecTriangleVertexIndex[riIndex][(i + 1) % TRIANGLE_ELEMENT_COUNT]].fY - rvecVertex[rvecTriangleVertexIndex[riIndex][i]].fY, rvecVertex[rvecTriangleVertexIndex[riIndex][(i + 1) % TRIANGLE_ELEMENT_COUNT]].fZ - rvecVertex[rvecTriangleVertexIndex[riIndex][i]].fZ};
    fFactor[i] = DotProduct(Vertex3D{rVertex3DData.fX - rvecVertex[rvecTriangleVertexIndex[riIndex][i]].fX, rVertex3DData.fY - rvecVertex[rvecTriangleVertexIndex[riIndex][i]].fY, rVertex3DData.fZ - rvecVertex[rvecTriangleVertexIndex[riIndex][i]].fZ}, Vertex3DEdge[i]) / DotProduct(Vertex3DEdge[i], Vertex3DEdge[i]);
  }

  float fSignedArea[TRIANGLE_ELEMENT_COUNT];
  float fGrossArea = 0.0F;

  for (int i = 0; i < TRIANGLE_ELEMENT_COUNT; i++)
    if (fFactor[i] <= 0.0F && fFactor[(i + TRIANGLE_ELEMENT_COUNT - 1) % TRIANGLE_ELEMENT_COUNT] >= 1.0F)
      return copysign(Distance(rvecVertex[rvecTriangleVertexIndex[riIndex][i]], rVertex3DData), DotProduct(Vertex3D{rVertex3DData.fX - rvecVertex[rvecTriangleVertexIndex[riIndex][i]].fX, rVertex3DData.fY - rvecVertex[rvecTriangleVertexIndex[riIndex][i]].fY, rVertex3DData.fZ - rvecVertex[rvecTriangleVertexIndex[riIndex][i]].fZ}, rvecVertexNormal[rvecTriangleVertexIndex[riIndex][i]]));
    else {
      fGrossArea += fSignedArea[i] = DotProduct(CrossProduct(Vertex3D{rvecVertex[rvecTriangleVertexIndex[riIndex][i]].fX - rVertex3DData.fX, rvecVertex[rvecTriangleVertexIndex[riIndex][i]].fY - rVertex3DData.fY, rvecVertex[rvecTriangleVertexIndex[riIndex][i]].fZ - rVertex3DData.fZ}, Vertex3D{rvecVertex[rvecTriangleVertexIndex[riIndex][(i + 1) % TRIANGLE_ELEMENT_COUNT]].fX - rVertex3DData.fX, rvecVertex[rvecTriangleVertexIndex[riIndex][(i + 1) % TRIANGLE_ELEMENT_COUNT]].fY - rVertex3DData.fY, rvecVertex[rvecTriangleVertexIndex[riIndex][(i + 1) % TRIANGLE_ELEMENT_COUNT]].fZ - rVertex3DData.fZ}), rvecSurfaceNormal[rvecTriangleNormalIndex[riIndex][i]]);

      if (fFactor[i] > 0.0F && fFactor[i] < 1.0F && fSignedArea[i] <= 0.0F) {
        Vertex3D Vertex3DClosestPoint{rvecVertex[rvecTriangleVertexIndex[riIndex][i]].fX + Vertex3DEdge[i].fX * fFactor[i], rvecVertex[rvecTriangleVertexIndex[riIndex][i]].fY + Vertex3DEdge[i].fY * fFactor[i], rvecVertex[rvecTriangleVertexIndex[riIndex][i]].fZ + Vertex3DEdge[i].fZ * fFactor[i]}, Vertex3DEdgeNormal{rvecVertexNormal[rvecTriangleVertexIndex[riIndex][i]].fX + rvecVertexNormal[rvecTriangleVertexIndex[riIndex][(i + 1) % TRIANGLE_ELEMENT_COUNT]].fX, rvecVertexNormal[rvecTriangleVertexIndex[riIndex][i]].fY + rvecVertexNormal[rvecTriangleVertexIndex[riIndex][(i + 1) % TRIANGLE_ELEMENT_COUNT]].fY, rvecVertexNormal[rvecTriangleVertexIndex[riIndex][i]].fZ + rvecVertexNormal[rvecTriangleVertexIndex[riIndex][(i + 1) % TRIANGLE_ELEMENT_COUNT]].fZ};

        return copysign(Distance(rVertex3DData, Vertex3DClosestPoint), DotProduct(Vertex3D{rVertex3DData.fX - Vertex3DClosestPoint.fX, rVertex3DData.fY - Vertex3DClosestPoint.fY, rVertex3DData.fZ - Vertex3DClosestPoint.fZ}, Vertex3DEdgeNormal));
      }
    }

  Vertex3D Vertex3DClosestPoint{0.0F, 0.0F, 0.0F};

  for (int i = 0; i < TRIANGLE_ELEMENT_COUNT; i++) {
    Vertex3DClosestPoint.fX += rvecVertex[rvecTriangleVertexIndex[riIndex][i]].fX * fSignedArea[(i + 1) % TRIANGLE_ELEMENT_COUNT] / fGrossArea;
    Vertex3DClosestPoint.fY += rvecVertex[rvecTriangleVertexIndex[riIndex][i]].fY * fSignedArea[(i + 1) % TRIANGLE_ELEMENT_COUNT] / fGrossArea;
    Vertex3DClosestPoint.fZ += rvecVertex[rvecTriangleVertexIndex[riIndex][i]].fZ * fSignedArea[(i + 1) % TRIANGLE_ELEMENT_COUNT] / fGrossArea;
  }

  return copysign(Distance(rVertex3DData, Vertex3DClosestPoint), DotProduct(Vertex3D{rVertex3DData.fX - Vertex3DClosestPoint.fX, rVertex3DData.fY - Vertex3DClosestPoint.fY, rVertex3DData.fZ - Vertex3DClosestPoint.fZ}, rvecSurfaceNormal[*rvecTriangleNormalIndex[riIndex].begin()]));
}

Vertex2D EdgeNormal(const Vertex2D& rVertex2DVector)
{
  return Normalize(Vertex2D{rVertex2DVector.fY, -rVertex2DVector.fX});
}

Vertex3D SurfaceNormal(const Surface3D& rSurface3DData, const float& rfEpsilon)
{
  Vertex3D Vertex3DResult{0.0F, 0.0F, 0.0F};
  HalfEdge3D* pHalfEdge3DAnchor = rSurface3DData.pHalfEdge3DPointer;
  int iCount = 0;

  if (pHalfEdge3DAnchor->pHalfEdge3DNext == pHalfEdge3DAnchor || pHalfEdge3DAnchor->pHalfEdge3DNext == pHalfEdge3DAnchor->pHalfEdge3DPrevious)
    throw invalid_argument("ポリゴンを構成するのに必要な頂点数が不足します。");
  else {
    HalfEdge3D* pHalfEdge3DCurrent = pHalfEdge3DAnchor;

    do {
      Vertex3D Vertex3DPartialSurfaceNormal = CrossProduct(Vertex3D{pHalfEdge3DCurrent->pHalfEdge3DNext->pVertex3DPointer->fX - pHalfEdge3DCurrent->pVertex3DPointer->fX, pHalfEdge3DCurrent->pHalfEdge3DNext->pVertex3DPointer->fY - pHalfEdge3DCurrent->pVertex3DPointer->fY, pHalfEdge3DCurrent->pHalfEdge3DNext->pVertex3DPointer->fZ - pHalfEdge3DCurrent->pVertex3DPointer->fZ}, Vertex3D{pHalfEdge3DCurrent->pHalfEdge3DPrevious->pVertex3DPointer->fX - pHalfEdge3DCurrent->pVertex3DPointer->fX, pHalfEdge3DCurrent->pHalfEdge3DPrevious->pVertex3DPointer->fY - pHalfEdge3DCurrent->pVertex3DPointer->fY, pHalfEdge3DCurrent->pHalfEdge3DPrevious->pVertex3DPointer->fZ - pHalfEdge3DCurrent->pVertex3DPointer->fZ});
      float fNorm = Norm(Vertex3DPartialSurfaceNormal);

      if (fNorm >= rfEpsilon) {
        Vertex3DResult.fX += Vertex3DPartialSurfaceNormal.fX / fNorm;
        Vertex3DResult.fY += Vertex3DPartialSurfaceNormal.fY / fNorm;
        Vertex3DResult.fZ += Vertex3DPartialSurfaceNormal.fZ / fNorm;
        iCount++;
      }

      pHalfEdge3DCurrent = pHalfEdge3DCurrent->pHalfEdge3DNext;
    } while (pHalfEdge3DCurrent != pHalfEdge3DAnchor);
  }

  if (!iCount)
    throw invalid_argument("ポリゴンが点あるいは線分に縮退しています。");

  return Normalize(Vertex3DResult);
}