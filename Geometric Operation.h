#ifndef GEOMETRIC_OPERATION
#define GEOMETRIC_OPERATION

#include <vector>
#include <limits>
#include <stdexcept>
#include <cmath>
#include "Definition.h"
#include "Data Structure.h"

#define GET_ARRAY_SIZE(varArray) sizeof(*GetArraySize(&(varArray)))

template <typename Type, std::size_t szSize> inline char (*GetArraySize(Type (*pvarArray)[szSize]))[szSize];

template <typename Type> inline Type Norm(const Type* const& rpvarArray, const std::size_t& rszSize)
{
  Type varResult = static_cast<Type>(0);

  for (int i = 0; i < rszSize; i++)
    varResult += static_cast<Type>(std::pow(static_cast<float>(rpvarArray[i]), 2.0F));

  return static_cast<Type>(std::sqrt(static_cast<float>(varResult)));
}

template <typename Type> inline Type Norm(const std::vector<Type>& rvecData)
{
  Type varResult = static_cast<Type>(0);

  for (int i = 0; i < rvecData.size(); i++)
    varResult += static_cast<Type>(std::pow(static_cast<float>(rvecData[i]), 2.0F));

  return static_cast<Type>(std::sqrt(static_cast<float>(varResult)));
}

template <typename Type> inline Type DotProduct(const Type* const& rpvarArrayA, const Type* const& rpvarArrayB, const std::size_t& rszSize)
{
  Type varResult = static_cast<Type>(0);

  for (int i = 0; i < rszSize; i++)
    varResult += rpvarArrayA[i] * rpvarArrayB[i];

  return varResult;
}

template <typename Type> inline Type DotProduct(const std::vector<Type>& rvecDataA, const std::vector<Type>& rvecDataB)
{
  Type varResult = static_cast<Type>(0);

  if (rvecDataA.size() == rvecDataB.size())
    for (int i = 0; i < rvecDataA.size() && i < rvecDataB.size(); i++)
      varResult += rvecDataA[i] * rvecDataB[i];
  else
    throw std::length_error("内積演算に使われるベクトルの次元数が異なります。");

  return varResult;
}

template <typename Type> inline void CrossProduct(const Type* const& rpvarArrayA, const Type* const& rpvarArrayB, Type* const& rpvarResult, const std::size_t& rszSize)
{
  if (rszSize == ELEMENT_COUNT_3D) {
    rpvarResult[INDICATOR_X] = rpvarArrayA[INDICATOR_Y] * rpvarArrayB[INDICATOR_Z] - rpvarArrayA[INDICATOR_Z] * rpvarArrayB[INDICATOR_Y];
    rpvarResult[INDICATOR_Y] = rpvarArrayA[INDICATOR_Z] * rpvarArrayB[INDICATOR_X] - rpvarArrayA[INDICATOR_X] * rpvarArrayB[INDICATOR_Z];
    rpvarResult[INDICATOR_Z] = rpvarArrayA[INDICATOR_X] * rpvarArrayB[INDICATOR_Y] - rpvarArrayA[INDICATOR_Y] * rpvarArrayB[INDICATOR_X];
  } else
    throw std::length_error("外積演算に使われるベクトルの次元数が不正です。");
}

template <typename Type> inline void CrossProduct(const std::vector<Type>& rvecDataA, const std::vector<Type>& rvecDataB, std::vector<Type>& rvecResult)
{
  if (rvecDataA.size() == ELEMENT_COUNT_3D && rvecDataB.size() == ELEMENT_COUNT_3D) {
    rvecResult.resize(ELEMENT_COUNT_3D);

    rvecResult[INDICATOR_X] = rvecDataA[INDICATOR_Y] * rvecDataB[INDICATOR_Z] - rvecDataA[INDICATOR_Z] * rvecDataB[INDICATOR_Y];
    rvecResult[INDICATOR_Y] = rvecDataA[INDICATOR_Z] * rvecDataB[INDICATOR_X] - rvecDataA[INDICATOR_X] * rvecDataB[INDICATOR_Z];
    rvecResult[INDICATOR_Z] = rvecDataA[INDICATOR_X] * rvecDataB[INDICATOR_Y] - rvecDataA[INDICATOR_Y] * rvecDataB[INDICATOR_X];
  } else
    throw std::length_error("外積演算に使われるベクトルの次元数が不正です。");
}

template <typename Type> inline void Normalize(const Type* const& rpvarArray, Type* const& rpvarResult, const std::size_t& rszSize)
{
  Type varNorm = Norm(rpvarArray, rszSize);

  for (int i = 0; i < rszSize; i++)
    rpvarResult[i] = rpvarArray[i] / varNorm;
}

template <typename Type> inline std::vector<Type> Normalize(const std::vector<Type>& rvecData)
{
  std::vector<Type> vecResult(rvecData.size());
  Type varNorm = Norm(rvecData);

  for (int i = 0; i < rvecData.size(); i++)
    vecResult[i] = rvecData[i] / varNorm;

  return vecResult;
}

template <typename Type, std::size_t szSize> inline void TransformCoordinate(const Type (*const& rpvarTransformationMatrix)[szSize], const Type* const& rpvarCoordinate, Type* const& rpvarResult)
{
  for (int i = 0; i < szSize - 1; i++) {
    rpvarResult[i] = 0.0F;

    for (int j = 0; j < szSize - 1; j++)
      rpvarResult[i] += rpvarTransformationMatrix[i][j] * rpvarCoordinate[j];

    rpvarResult[i] += rpvarTransformationMatrix[i][szSize - 1];
  }
}

template <typename Type, std::size_t szSize> inline void TransformCoordinate(const Type (*const& rpvarTransformationMatrix)[szSize], const std::vector<Type>& rvecVertex, std::vector<Type>& rvecTransformedVertex)
{
  rvecTransformedVertex.resize(rvecVertex.size());

  for (int i = 0; i < rvecVertex.size(); i++)
    TransformCoordinate(rpvarTransformationMatrix, &rvecVertex[i], &rvecTransformedVertex[i]);
}

template <typename Type, std::size_t szSizeA, std::size_t szSizeB, std::size_t szSizeC> inline void MultiplyMatrix(const Type varMatrixA[szSizeA][szSizeB], const Type varMatrixB[szSizeB][szSizeC], Type varProduct[szSizeA][szSizeC])
{
  for (int i = 0; i < szSizeA; i++)
    for (int j = 0; j < szSizeC; j++) {
      varProduct[i][j] = 0.0F;

      for (int k = 0; k < szSizeB; k++)
        varProduct[i][j] += varMatrixA[i][k] * varMatrixB[k][j];
    }
}

template <std::size_t szSize> inline void InverseTransformation(float (*const& rpfTransformationMatrix)[szSize])
{
  float fTranslationData[szSize - 1];

  for (int i = 0; i < szSize - 1; i++) {
    fTranslationData[i] = rpfTransformationMatrix[i][szSize - 1];
    rpfTransformationMatrix[i][szSize - 1] = 0.0F;
  }

  for (int i = 0; i < szSize - 1; i++) {
    for (int j = i + 1; j < szSize - 1; j++) {
      float fTemp;

      fTemp = rpfTransformationMatrix[i][j];
      rpfTransformationMatrix[i][j] = rpfTransformationMatrix[j][i];
      rpfTransformationMatrix[j][i] = fTemp;
    }

    for (int j = 0; j < szSize - 1; j++)
      rpfTransformationMatrix[i][szSize - 1] -= rpfTransformationMatrix[i][j] * fTranslationData[j];
  }
}

float Norm(const Vertex2D& rVertex2DVector);

float Norm(const Vertex3D& rVertex3DVector);

float DotProduct(const Vertex2D& rVertex2DVectorA, const Vertex2D& rVertex2DVectorB);

float DotProduct(const Vertex3D& rVertex3DVectorA, const Vertex3D& rVertex3DVectorB);

Vertex2D CrossProduct(const Vertex2D& rVertex2DVectorA, const Vertex2D& rVertex2DVectorB);

Vertex3D CrossProduct(const Vertex3D& rVertex3DVectorA, const Vertex3D& rVertex3DVectorB);

Vertex2D Normalize(const Vertex2D& rVertex2DVector);

Vertex3D Normalize(const Vertex3D& rVertex3DVector);

Vertex2D TransformCoordinate(const float (*const& rpfTransformationMatrix)[ELEMENT_COUNT_2D + 1], const Vertex2D& rVertex2DData);

Vertex3D TransformCoordinate(const float (*const& rpfTransformationMatrix)[ELEMENT_COUNT_3D + 1], const Vertex3D& rVertex3DData);

void TransformCoordinate(const float (*const& rpfTransformationMatrix)[ELEMENT_COUNT_2D + 1], const std::vector<Vertex2D>& rvecVertex, std::vector<Vertex2D>& rvecTransformedVertex);

void TransformCoordinate(const float (*const& rpfTransformationMatrix)[ELEMENT_COUNT_3D + 1], const std::vector<Vertex3D>& rvecVertex, std::vector<Vertex3D>& rvecTransformedVertex);

void RandomizeTransformation(float (*const& rpfTransformationMatrix)[ELEMENT_COUNT_2D + 1], const Vertex2D& rVertex2DMinimumCorner, const Vertex2D& rVertex2DMaximumCorner);

void RandomizeTransformation(float (*const& rpfTransformationMatrix)[ELEMENT_COUNT_3D + 1], const Vertex3D& rVertex3DMinimumCorner, const Vertex3D& rVertex3DMaximumCorner);

float Distance(const Vertex2D& rVertex2DDataA, const Vertex2D& rVertex2DDataB);

float Distance(const Vertex2D& rVertex2DData, const Vertex2D& rVertex2DEndpointA, const Vertex2D& rVertex2DEndpointB);

float Distance(const std::vector<Vertex2D>& rvecVertex, const std::vector<Vertex2D>& rvecVertexNormal, const Vertex2D& rVertex2DData, const int& riIndex);

float Distance(const Vertex3D& rVertex3DDataA, const Vertex3D& rVertex3DDataB);

float Distance(const Vertex3D& rVertex3DData, const Vertex3D& rVertex3DEndpointA, const Vertex3D& rVertex3DEndpointB);

float Distance(const std::vector<Vertex3D>& rvecVertex, const std::vector<Vertex3D>& rvecVertexNormal, const std::vector<Vertex3D>& rvecSurfaceNormal, const std::vector<std::vector<unsigned>>& rvecTriangleVertexIndex, const std::vector<std::vector<unsigned>>& rvecTriangleNormalIndex, const Vertex3D& rVertex3DData, const int& riIndex);

Vertex2D EdgeNormal(const Vertex2D& rVertex2DVector);

Vertex3D SurfaceNormal(const Surface3D& rSurface3DData, const float& rfEpsilon = std::numeric_limits<float>::epsilon());

#endif