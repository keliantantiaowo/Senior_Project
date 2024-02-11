#include <iterator>
#include <limits>
#include "Display Result.h"
#include "Geometric Operation.h"
#include "Equation Solver.h"
#include "Generate Implicit Surface.h"

using namespace std;

void CalculateVertexNormal2D(const vector<Vertex2D>& rvecVertex, vector<Vertex2D>& rvecVertexNormal, const NormalStyle& rNormalStyleSetting)
{
  rvecVertexNormal.clear();

  for (int i = 0; i < rvecVertex.size(); i++) {
    Vertex2D Vertex2DEdgeNormalA{Normalize(Vertex2D{rvecVertex[(i + 1) % rvecVertex.size()].fX - rvecVertex[i].fX, rvecVertex[(i + 1) % rvecVertex.size()].fY - rvecVertex[i].fY})}, Vertex2DEdgeNormalB{Normalize(Vertex2D{rvecVertex[(i + rvecVertex.size() - 1) % rvecVertex.size()].fX - rvecVertex[i].fX, rvecVertex[(i + rvecVertex.size() - 1) % rvecVertex.size()].fY - rvecVertex[i].fY})};

    rvecVertexNormal.push_back(EdgeNormal(Vertex2D{Vertex2DEdgeNormalA.fX - Vertex2DEdgeNormalB.fX, Vertex2DEdgeNormalA.fY - Vertex2DEdgeNormalB.fY}));

    switch (rNormalStyleSetting) {
    case Normalized : break;
    case Weighted :
      do {
        float fFactor = 2.0F / (DotProduct(rvecVertexNormal[i], Vertex2DEdgeNormalA) + DotProduct(rvecVertexNormal[i], Vertex2DEdgeNormalB));

        rvecVertexNormal[i].fX *= fFactor;
        rvecVertexNormal[i].fY *= fFactor;
      } while (false);

      break;
    default : break;
    }
  }

  ShowElement(rvecVertexNormal, "rvecVertexNormal", Newline);
}

void CalculateVertexNormal3D(const vector<Vertex3D>& rvecVertex, vector<Vertex3D>& rvecVertexNormal, const vector<Vertex3D>& rvecSurfaceNormal, const vector<vector<unsigned>>& rvecTriangleVertexIndex, const vector<vector<unsigned>>& rvecTriangleNormalIndex, const NormalStyle& rNormalStyleSetting)
{
  vector<vector<WeightedNormal>> vecSharedNormalIndex(rvecVertex.size(), vector<WeightedNormal>(0));

  for (int i = 0; i < rvecTriangleVertexIndex.size(); i++)
    for (int j = 0; j < TRIANGLE_ELEMENT_COUNT; j++) {
      Vertex3D Vertex3DNormalizedEdgeA{Normalize(Vertex3D{rvecVertex[rvecTriangleVertexIndex[i][(j + 1) % TRIANGLE_ELEMENT_COUNT]].fX - rvecVertex[rvecTriangleVertexIndex[i][j]].fX, rvecVertex[rvecTriangleVertexIndex[i][(j + 1) % TRIANGLE_ELEMENT_COUNT]].fY - rvecVertex[rvecTriangleVertexIndex[i][j]].fY, rvecVertex[rvecTriangleVertexIndex[i][(j + 1) % TRIANGLE_ELEMENT_COUNT]].fZ - rvecVertex[rvecTriangleVertexIndex[i][j]].fZ})}, Vertex3DNormalizedEdgeB{Normalize(Vertex3D{rvecVertex[rvecTriangleVertexIndex[i][(j + TRIANGLE_ELEMENT_COUNT - 1) % TRIANGLE_ELEMENT_COUNT]].fX - rvecVertex[rvecTriangleVertexIndex[i][j]].fX, rvecVertex[rvecTriangleVertexIndex[i][(j + TRIANGLE_ELEMENT_COUNT - 1) % TRIANGLE_ELEMENT_COUNT]].fY - rvecVertex[rvecTriangleVertexIndex[i][j]].fY, rvecVertex[rvecTriangleVertexIndex[i][(j + TRIANGLE_ELEMENT_COUNT - 1) % TRIANGLE_ELEMENT_COUNT]].fZ - rvecVertex[rvecTriangleVertexIndex[i][j]].fZ})};

      vecSharedNormalIndex[rvecTriangleVertexIndex[i][j]].push_back(WeightedNormal{rvecTriangleNormalIndex[i][j], acosf(DotProduct(Vertex3DNormalizedEdgeA, Vertex3DNormalizedEdgeB))});
    }

  rvecVertexNormal.assign(rvecVertex.size(), Vertex3D{0.0F, 0.0F, 0.0F});

  for (int i = 0; i < rvecVertexNormal.size(); i++) {
    for (int j = 0; j < vecSharedNormalIndex[i].size(); j++) {
      rvecVertexNormal[i].fX += rvecSurfaceNormal[vecSharedNormalIndex[i][j].uIndex].fX * vecSharedNormalIndex[i][j].fAngle;
      rvecVertexNormal[i].fY += rvecSurfaceNormal[vecSharedNormalIndex[i][j].uIndex].fY * vecSharedNormalIndex[i][j].fAngle;
      rvecVertexNormal[i].fZ += rvecSurfaceNormal[vecSharedNormalIndex[i][j].uIndex].fZ * vecSharedNormalIndex[i][j].fAngle;
    }
    
    rvecVertexNormal[i] = Normalize(rvecVertexNormal[i]);

    switch (rNormalStyleSetting) {
    case Normalized : break;
    case Weighted :
      do {
        float fExtremeValue = 0.0F;

        for (int j = 0; j < vecSharedNormalIndex[i].size(); j++) {
          float fFactor = 1.0F / DotProduct(rvecSurfaceNormal[vecSharedNormalIndex[i][j].uIndex], rvecVertexNormal[i]);

          if (fExtremeValue < fFactor)
            fExtremeValue = fFactor;
        }

        rvecVertexNormal[i].fX *= fExtremeValue;
        rvecVertexNormal[i].fY *= fExtremeValue;
        rvecVertexNormal[i].fZ *= fExtremeValue;
      } while (false);

      break;
    default : break;
    }
  }

  ShowElement(vecSharedNormalIndex, "vecSharedNormalIndex", Newline);
  ShowElement(rvecVertexNormal, "rvecVertexNormal", Newline);
}

void InsertVertex2D(vector<Vertex2D>& rvecVertex, vector<Vertex2D>& rvecVertexNormal, const int& riInsertionCount, const int& riBias)
{
  if (!riInsertionCount)
    return ;

  const int iOriginalSize = rvecVertex.size();

  if (riBias) {
    const int iAdditionalSize = (riInsertionCount + EDGE_ELEMENT_COUNT) * iOriginalSize, iAdditionalSizePerEdge = riInsertionCount + EDGE_ELEMENT_COUNT;

    rvecVertex.resize(iOriginalSize + iAdditionalSize);
    rvecVertexNormal.resize(iOriginalSize + iAdditionalSize);

    for (int i = 0; i < iOriginalSize; i++) {
      rvecVertex[rvecVertex.size() - (iAdditionalSizePerEdge + 1) * (i + 1)] = rvecVertex[iOriginalSize - i - 1];
      rvecVertexNormal[rvecVertexNormal.size() - (iAdditionalSizePerEdge + 1) * (i + 1)] = rvecVertexNormal[iOriginalSize - i - 1];
    }

    float fFactor = 1.0F / ((static_cast<float>(riInsertionCount) + 1.0F) * (static_cast<float>(riBias) + 1.0F));

    for (int i = 0, j = iOriginalSize - 1; i < iOriginalSize; j = i++) {
      Vertex2D Vertex2DDifference{(rvecVertex[(iAdditionalSizePerEdge + 1) * i].fX - rvecVertex[(iAdditionalSizePerEdge + 1) * j].fX) / (static_cast<float>(riInsertionCount) + 1.0F), (rvecVertex[(iAdditionalSizePerEdge + 1) * i].fY - rvecVertex[(iAdditionalSizePerEdge + 1) * j].fY) / (static_cast<float>(riInsertionCount) + 1.0F)}, Vertex2DNormal{EdgeNormal(Vertex2DDifference)};

      for (int k = 0; k < iAdditionalSizePerEdge; k++) {
        if (!k)
          rvecVertex[(iAdditionalSizePerEdge + 1) * j + k + 1] = Vertex2D{rvecVertex[(iAdditionalSizePerEdge + 1) * j].fX + Vertex2DDifference.fX * fFactor, rvecVertex[(iAdditionalSizePerEdge + 1) * j].fY + Vertex2DDifference.fY * fFactor};
        else if (!(iAdditionalSizePerEdge - k - 1))
          rvecVertex[(iAdditionalSizePerEdge + 1) * j + k + 1] = Vertex2D{rvecVertex[(iAdditionalSizePerEdge + 1) * i].fX - Vertex2DDifference.fX * fFactor, rvecVertex[(iAdditionalSizePerEdge + 1) * i].fY - Vertex2DDifference.fY * fFactor};
        else
          rvecVertex[(iAdditionalSizePerEdge + 1) * j + k + 1] = Vertex2D{rvecVertex[(iAdditionalSizePerEdge + 1) * j].fX + Vertex2DDifference.fX * static_cast<float>(k), rvecVertex[(iAdditionalSizePerEdge + 1) * j].fY + Vertex2DDifference.fY * static_cast<float>(k)};

        rvecVertexNormal[(iAdditionalSizePerEdge + 1) * j + k + 1] = Vertex2DNormal;
      }
    }
  } else {
    rvecVertex.resize((riInsertionCount + 1) * iOriginalSize);
    rvecVertexNormal.resize((riInsertionCount + 1) * iOriginalSize);

    for (int i = 0; i < iOriginalSize; i++) {
      rvecVertex[rvecVertex.size() - (riInsertionCount + 1) * (i + 1)] = rvecVertex[iOriginalSize - i - 1];
      rvecVertexNormal[rvecVertexNormal.size() - (riInsertionCount + 1) * (i + 1)] = rvecVertexNormal[iOriginalSize - i - 1];
    }

    for (int i = 0, j = iOriginalSize - 1; i < iOriginalSize; j = i++) {
      Vertex2D Vertex2DDifference{(rvecVertex[(riInsertionCount + 1) * i].fX - rvecVertex[(riInsertionCount + 1) * j].fX) / (static_cast<float>(riInsertionCount) + 1.0F), (rvecVertex[(riInsertionCount + 1) * i].fY - rvecVertex[(riInsertionCount + 1) * j].fY) / (static_cast<float>(riInsertionCount) + 1.0F)}, Vertex2DNormal{EdgeNormal(Vertex2DDifference)};

      for (int k = 0; k < riInsertionCount; k++) {
        rvecVertex[(riInsertionCount + 1) * j + k + 1] = Vertex2D{rvecVertex[(riInsertionCount + 1) * j].fX + Vertex2DDifference.fX * (static_cast<float>(k) + 1.0F), rvecVertex[(riInsertionCount + 1) * j].fY + Vertex2DDifference.fY * (static_cast<float>(k) + 1.0F)};
        rvecVertexNormal[(riInsertionCount + 1) * j + k + 1] = Vertex2DNormal;
      }
    }
  }

  ShowElement(rvecVertex, "rvecVertex", Newline);
  ShowElement(rvecVertexNormal, "rvecVertexNormal", Newline);
}

void InsertVertex3D(vector<Vertex3D>& rvecVertex, vector<Vertex3D>& rvecVertexNormal, const vector<Vertex3D>& rvecSurfaceNormal, const vector<vector<unsigned>>& rvecTriangleVertexIndex, const vector<vector<unsigned>>& rvecTriangleNormalIndex, const int& riInsertionCount, const int& riBias)
{
  if (!riInsertionCount)
    return ;

  list<unsigned> listEndpointIndexA(0), listEndpointIndexB(0);

  for (int i = 0; i < rvecTriangleVertexIndex.size(); i++)
    for (int j = 0; j < TRIANGLE_ELEMENT_COUNT; j++) {
      bool bDuplication = false;

      for (list<unsigned>::iterator k = listEndpointIndexA.begin(), l = listEndpointIndexB.begin(); k != listEndpointIndexA.end() && l != listEndpointIndexB.end(); k++, l++)
        if (rvecTriangleVertexIndex[i][(j + 1) % TRIANGLE_ELEMENT_COUNT] == *k && rvecTriangleVertexIndex[i][j] == *l) {
          bDuplication = true;

          break;
        }

      if (!bDuplication) {
        listEndpointIndexA.push_back(rvecTriangleVertexIndex[i][j]);
        listEndpointIndexB.push_back(rvecTriangleVertexIndex[i][(j + 1) % TRIANGLE_ELEMENT_COUNT]);
      }
    }

  list<unsigned>::iterator iterStart = listEndpointIndexA.begin(), iterEnd = listEndpointIndexB.begin();
  const int iOriginalSize = rvecVertex.size();

  if (riBias) {
    const int iAdditionalSize = rvecTriangleVertexIndex.size() * riInsertionCount * TRIANGLE_ELEMENT_COUNT / 2 + rvecTriangleVertexIndex.size() * (riInsertionCount + 1) * TRIANGLE_ELEMENT_COUNT * 2, iAdditionalSizePerEdge = (riInsertionCount + 1) * 2, iAdditionalSizePerSurface = iAdditionalSizePerEdge * TRIANGLE_ELEMENT_COUNT;

    rvecVertex.resize(iOriginalSize + iAdditionalSize);
    rvecVertexNormal.resize(iOriginalSize + iAdditionalSize);

    float fExtremeValue = numeric_limits<float>::max();

    for (int i = 0; i < listEndpointIndexA.size() && i < listEndpointIndexB.size(); i++) {
      Vertex3D Vertex3DDifference{(rvecVertex[*iterEnd].fX - rvecVertex[*iterStart].fX) / (static_cast<float>(riInsertionCount) + 1.0F), (rvecVertex[*iterEnd].fY - rvecVertex[*iterStart].fY) / (static_cast<float>(riInsertionCount) + 1.0F), (rvecVertex[*iterEnd].fZ - rvecVertex[*iterStart].fZ) / (static_cast<float>(riInsertionCount) + 1.0F)}, Vertex3DNormal{Normalize(Vertex3D{rvecVertexNormal[*iterStart].fX + rvecVertexNormal[*iterEnd].fX, rvecVertexNormal[*iterStart].fY + rvecVertexNormal[*iterEnd].fY, rvecVertexNormal[*iterStart].fZ + rvecVertexNormal[*iterEnd].fZ})};

      for (int j = 0; j < riInsertionCount; j++) {
        rvecVertex[riInsertionCount * i + iOriginalSize + j] = Vertex3D{rvecVertex[*iterStart].fX + Vertex3DDifference.fX * (static_cast<float>(j) + 1.0F), rvecVertex[*iterStart].fY + Vertex3DDifference.fY * (static_cast<float>(j) + 1.0F), rvecVertex[*iterStart].fZ + Vertex3DDifference.fZ * (static_cast<float>(j) + 1.0F)};
        rvecVertexNormal[riInsertionCount * i + iOriginalSize + j] = Vertex3DNormal;
      }

      float fDistance = Distance(rvecVertex[*iterStart], rvecVertex[*iterEnd]);

      if (fExtremeValue > fDistance)
        fExtremeValue = fDistance;

      iterStart++;
      iterEnd++;
    }

    for (int i = 0; i < rvecTriangleVertexIndex.size(); i++) {
      unsigned uIndex = rvecVertex.size() - (rvecTriangleVertexIndex.size() - static_cast<unsigned>(i)) * static_cast<unsigned>(iAdditionalSizePerSurface);

      FillUpTriangle(rvecVertex, rvecVertexNormal, *rvecTriangleVertexIndex[i].begin(), *next(rvecTriangleVertexIndex[i].begin()), *prev(rvecTriangleVertexIndex[i].end()), uIndex, rvecSurfaceNormal[*rvecTriangleNormalIndex[i].begin()], riInsertionCount, riBias, fExtremeValue <= 1.0F ? fExtremeValue : 1.0F);
    }
  } else {
    const int iAdditionalSize = rvecTriangleVertexIndex.size() * riInsertionCount * TRIANGLE_ELEMENT_COUNT / 2 + rvecTriangleVertexIndex.size() * (static_cast<int>(pow(static_cast<float>(TRIANGLE_ELEMENT_COUNT), static_cast<float>(riInsertionCount))) - 1) / 2, iAdditionalSizePerEdge = 0, iAdditionalSizePerSurface = (static_cast<int>(pow(static_cast<float>(TRIANGLE_ELEMENT_COUNT), static_cast<float>(riInsertionCount))) - 1) / 2;

    rvecVertex.resize(iOriginalSize + iAdditionalSize);
    rvecVertexNormal.resize(iOriginalSize + iAdditionalSize);

    for (int i = 0; i < listEndpointIndexA.size() && i < listEndpointIndexB.size(); i++) {
      Vertex3D Vertex3DDifference{(rvecVertex[*iterEnd].fX - rvecVertex[*iterStart].fX) / (static_cast<float>(riInsertionCount) + 1.0F), (rvecVertex[*iterEnd].fY - rvecVertex[*iterStart].fY) / (static_cast<float>(riInsertionCount) + 1.0F), (rvecVertex[*iterEnd].fZ - rvecVertex[*iterStart].fZ) / (static_cast<float>(riInsertionCount) + 1.0F)}, Vertex3DNormal{Normalize(Vertex3D{rvecVertexNormal[*iterStart].fX + rvecVertexNormal[*iterEnd].fX, rvecVertexNormal[*iterStart].fY + rvecVertexNormal[*iterEnd].fY, rvecVertexNormal[*iterStart].fZ + rvecVertexNormal[*iterEnd].fZ})};

      for (int j = 0; j < riInsertionCount; j++) {
        rvecVertex[riInsertionCount * i + iOriginalSize + j] = Vertex3D{rvecVertex[*iterStart].fX + Vertex3DDifference.fX * (static_cast<float>(j) + 1.0F), rvecVertex[*iterStart].fY + Vertex3DDifference.fY * (static_cast<float>(j) + 1.0F), rvecVertex[*iterStart].fZ + Vertex3DDifference.fZ * (static_cast<float>(j) + 1.0F)};
        rvecVertexNormal[riInsertionCount * i + iOriginalSize + j] = Vertex3DNormal;
      }

      iterStart++;
      iterEnd++;
    }

    for (int i = 0; i < rvecTriangleVertexIndex.size() && i < rvecTriangleNormalIndex.size(); i++) {
      unsigned uIndex = rvecVertex.size() - (rvecTriangleVertexIndex.size() - static_cast<unsigned>(i)) * static_cast<unsigned>(iAdditionalSizePerSurface);

      FillUpTriangle(rvecVertex, rvecVertexNormal, *rvecTriangleVertexIndex[i].begin(), *next(rvecTriangleVertexIndex[i].begin()), *prev(rvecTriangleVertexIndex[i].end()), unsigned(uIndex), uIndex, rvecSurfaceNormal[*rvecTriangleNormalIndex[i].begin()], riInsertionCount);
    }
  }

  ShowElement(rvecVertex, "rvecVertex", Newline);
  ShowElement(rvecVertexNormal, "rvecVertexNormal", Newline);
}

float FunctionValue2D(const float* const& rpfCoordinate, const vector<Vertex2D>& rvecVertex, const vector<float>& rvecCoefficient)
{
  float fResult = 0.0F;

  for (int i = 0; i < rvecVertex.size(); i++) {
    float fDifference[ELEMENT_COUNT_2D] = {rpfCoordinate[INDICATOR_X] - rvecVertex[i].fX, rpfCoordinate[INDICATOR_Y] - rvecVertex[i].fY};

    fResult += pow(Norm(reinterpret_cast<const float* const>(fDifference), ELEMENT_COUNT_2D), 3.0F) * rvecCoefficient[i];
  }

  for (int i = 0; i < ELEMENT_COUNT_2D; i++)
    fResult += rpfCoordinate[i] * rvecCoefficient[rvecVertex.size() + i];

  fResult += rvecCoefficient[rvecCoefficient.size() - 1];
}

float FunctionValue3D(const float* const& rpfCoordinate, const vector<Vertex3D>& rvecVertex, const vector<float>& rvecCoefficient)
{
  float fResult = 0.0F;

  for (int i = 0; i < rvecVertex.size(); i++) {
    float fDifference[ELEMENT_COUNT_3D] = {rpfCoordinate[INDICATOR_X] - rvecVertex[i].fX, rpfCoordinate[INDICATOR_Y] - rvecVertex[i].fY, rpfCoordinate[INDICATOR_Z] - rvecVertex[i].fZ};

    fResult += pow(Norm(reinterpret_cast<const float* const>(fDifference), ELEMENT_COUNT_3D), 3.0F) * rvecCoefficient[i];
  }

  for (int i = 0; i < ELEMENT_COUNT_3D; i++)
    fResult += rpfCoordinate[i] * rvecCoefficient[rvecVertex.size() + i];

  fResult += rvecCoefficient[rvecCoefficient.size() - 1];
}

void GenerateImplicitSurface2D(vector<Vertex2D>& rvecVertex, const vector<Vertex2D>& rvecVertexNormal, vector<float>& rvecCoefficient, const float& rfEpsilon)
{
  float** pfMatrix = new float*[rvecVertex.size() * 2 + EXTRA_COEFFICIENT_COUNT_2D];
  const int iMatrixSize = rvecVertex.size() * 2 + EXTRA_COEFFICIENT_COUNT_2D;

  rvecCoefficient.resize(iMatrixSize);

  for (int i = 0; i < rvecVertexNormal.size(); i++)
    rvecVertex.push_back(Vertex2D{rvecVertex[i].fX - rvecVertexNormal[i].fX * rfEpsilon, rvecVertex[i].fY - rvecVertexNormal[i].fY * rfEpsilon});

  for (int i = 0; i < iMatrixSize; i++)
    pfMatrix[i] = new float[iMatrixSize];

  for (int i = 0; i < rvecVertex.size(); i++) {
    for (int j = 0; j < rvecVertex.size(); j++) {
      float fDifference[ELEMENT_COUNT_2D] = {rvecVertex[i].fX - rvecVertex[j].fX, rvecVertex[i].fY - rvecVertex[j].fY};

      pfMatrix[i][j] = pow(Norm(reinterpret_cast<const float* const>(fDifference), ELEMENT_COUNT_2D), 3.0F);
    }

    pfMatrix[i][rvecVertex.size() + INDICATOR_X] = pfMatrix[rvecVertex.size() + INDICATOR_X][i] = rvecVertex[i].fX;
    pfMatrix[i][rvecVertex.size() + INDICATOR_Y] = pfMatrix[rvecVertex.size() + INDICATOR_Y][i] = rvecVertex[i].fY;
    pfMatrix[i][iMatrixSize - 1] = pfMatrix[iMatrixSize - 1][i] = 1.0F;
    rvecCoefficient[i] = rvecVertexNormal.size() > i ? 0.0F : rfEpsilon;
  }

  for (int i = 0; i < EXTRA_COEFFICIENT_COUNT_2D; i++) {
    for (int j = 0; j < EXTRA_COEFFICIENT_COUNT_2D; j++)
      pfMatrix[rvecVertex.size() + i][rvecVertex.size() + j] = 0.0F;

    rvecCoefficient[rvecVertex.size() + i] = 0.0F;
  }

  unsigned* puPivot = new unsigned[iMatrixSize];

  ShowElement(rvecCoefficient, "rvecCoefficient", Newline);

  SolveEquation(pfMatrix, puPivot, rvecCoefficient, iMatrixSize);

  ConfirmSolution(pfMatrix, rvecCoefficient, iMatrixSize, NoLineBreak);

  delete[] puPivot;

  for (int i = 0; i < iMatrixSize; i++)
    delete[] pfMatrix[i];

  delete[] pfMatrix;
}

void GenerateImplicitSurface3D(vector<Vertex3D>& rvecVertex, const vector<Vertex3D>& rvecVertexNormal, vector<float>& rvecCoefficient, const float& rfEpsilon)
{
  float** pfMatrix = new float*[rvecVertex.size() * 2 + EXTRA_COEFFICIENT_COUNT_3D];
  const int iMatrixSize = rvecVertex.size() * 2 + EXTRA_COEFFICIENT_COUNT_3D;

  rvecCoefficient.resize(iMatrixSize);

  for (int i = 0; i < rvecVertexNormal.size(); i++)
    rvecVertex.push_back(Vertex3D{rvecVertex[i].fX - rvecVertexNormal[i].fX * rfEpsilon, rvecVertex[i].fY - rvecVertexNormal[i].fY * rfEpsilon, rvecVertex[i].fZ - rvecVertexNormal[i].fZ * rfEpsilon});

  for (int i = 0; i < iMatrixSize; i++)
    pfMatrix[i] = new float[iMatrixSize];

  for (int i = 0; i < rvecVertex.size(); i++) {
    for (int j = 0; j < rvecVertex.size(); j++) {
      float fDifference[ELEMENT_COUNT_3D] = {rvecVertex[i].fX - rvecVertex[j].fX, rvecVertex[i].fY - rvecVertex[j].fY, rvecVertex[i].fZ - rvecVertex[j].fZ};

      pfMatrix[i][j] = pow(Norm(reinterpret_cast<const float* const>(fDifference), ELEMENT_COUNT_3D), 3.0F);
    }

    pfMatrix[i][rvecVertex.size() + INDICATOR_X] = pfMatrix[rvecVertex.size() + INDICATOR_X][i] = rvecVertex[i].fX;
    pfMatrix[i][rvecVertex.size() + INDICATOR_Y] = pfMatrix[rvecVertex.size() + INDICATOR_Y][i] = rvecVertex[i].fY;
    pfMatrix[i][rvecVertex.size() + INDICATOR_Z] = pfMatrix[rvecVertex.size() + INDICATOR_Z][i] = rvecVertex[i].fZ;
    pfMatrix[i][iMatrixSize - 1] = pfMatrix[iMatrixSize - 1][i] = 1.0F;
    rvecCoefficient[i] = rvecVertexNormal.size() > i ? 0.0F : rfEpsilon;
  }
  
  for (int i = 0; i < EXTRA_COEFFICIENT_COUNT_3D; i++) {
    for (int j = 0; j < EXTRA_COEFFICIENT_COUNT_3D; j++)
      pfMatrix[rvecVertex.size() + i][rvecVertex.size() + j] = 0.0F;

    rvecCoefficient[rvecVertex.size() + i] = 0.0F;
  }

  unsigned* puPivot = new unsigned[iMatrixSize];

  ShowElement(rvecCoefficient, "rvecCoefficient", Newline);

  SolveEquation(pfMatrix, puPivot, rvecCoefficient, iMatrixSize);

  ConfirmSolution(pfMatrix, rvecCoefficient, iMatrixSize, NoLineBreak);

  delete[] puPivot;

  for (int i = 0; i < iMatrixSize; i++)
    delete[] pfMatrix[i];

  delete[] pfMatrix;
}

void FillUpTriangle(vector<Vertex3D>& rvecVertex, vector<Vertex3D>& rvecVertexNormal, const unsigned& ruStart, const unsigned& ruIntermediate, const unsigned& ruEnd, unsigned& ruIndex, const Vertex3D& rVertex3DNormal, const int& riInsertionCount, const int& riBias, const float& rfFixedLength)
{
  Vertex3D Vertex3DNormalizedEdge{0.0F, 0.0F, 0.0F}, Vertex3DEdgeNormal{0.0F, 0.0F, 0.0F};
  float fGlobalFactor = 1.0F / ((static_cast<float>(riInsertionCount) + 1.0F) * (static_cast<float>(riBias) + 1.0F));

  Vertex3DNormalizedEdge = Normalize(Vertex3D{rvecVertex[ruIntermediate].fX - rvecVertex[ruStart].fX, rvecVertex[ruIntermediate].fY - rvecVertex[ruStart].fY, rvecVertex[ruIntermediate].fZ - rvecVertex[ruStart].fZ});
  Vertex3DEdgeNormal = CrossProduct(Vertex3DNormalizedEdge, rVertex3DNormal);
  rvecVertex[ruIndex] = Vertex3D{rvecVertex[ruStart].fX + Vertex3DNormalizedEdge.fX * rfFixedLength * fGlobalFactor / 2.0F - Vertex3DEdgeNormal.fX * rfFixedLength * pow(fGlobalFactor, 2.0F), rvecVertex[ruStart].fY + Vertex3DNormalizedEdge.fY * rfFixedLength * fGlobalFactor / 2.0F - Vertex3DEdgeNormal.fY * rfFixedLength * pow(fGlobalFactor, 2.0F), rvecVertex[ruStart].fZ + Vertex3DNormalizedEdge.fZ * rfFixedLength * fGlobalFactor / 2.0F - Vertex3DEdgeNormal.fZ * rfFixedLength * pow(fGlobalFactor, 2.0F)};
  rvecVertexNormal[ruIndex] = rVertex3DNormal;
  ruIndex++;

  for (int i = 0; i < riInsertionCount; i++) {
    float fLocalFactor = (static_cast<float>(i) + 1.0F) / (static_cast<float>(riInsertionCount) + 1.0F);
    Vertex3D Vertex3DEndpoint{rvecVertex[ruStart].fX * (1.0F - fLocalFactor) + rvecVertex[ruIntermediate].fX * fLocalFactor, rvecVertex[ruStart].fY * (1.0F - fLocalFactor) + rvecVertex[ruIntermediate].fY * fLocalFactor, rvecVertex[ruStart].fZ * (1.0F - fLocalFactor) + rvecVertex[ruIntermediate].fZ * fLocalFactor};

    rvecVertex[ruIndex] = Vertex3D{Vertex3DEndpoint.fX - Vertex3DNormalizedEdge.fX * rfFixedLength * fGlobalFactor - Vertex3DEdgeNormal.fX * rfFixedLength * fGlobalFactor, Vertex3DEndpoint.fY - Vertex3DNormalizedEdge.fY * rfFixedLength * fGlobalFactor - Vertex3DEdgeNormal.fY * rfFixedLength * fGlobalFactor, Vertex3DEndpoint.fZ - Vertex3DNormalizedEdge.fZ * rfFixedLength * fGlobalFactor - Vertex3DEdgeNormal.fZ * rfFixedLength * fGlobalFactor};
    rvecVertexNormal[ruIndex] = rVertex3DNormal;
    ruIndex++;
    rvecVertex[ruIndex] = Vertex3D{Vertex3DEndpoint.fX + Vertex3DNormalizedEdge.fX * rfFixedLength * fGlobalFactor - Vertex3DEdgeNormal.fX * rfFixedLength * fGlobalFactor, Vertex3DEndpoint.fY + Vertex3DNormalizedEdge.fY * rfFixedLength * fGlobalFactor - Vertex3DEdgeNormal.fY * rfFixedLength * fGlobalFactor, Vertex3DEndpoint.fZ + Vertex3DNormalizedEdge.fZ * rfFixedLength * fGlobalFactor - Vertex3DEdgeNormal.fZ * rfFixedLength * fGlobalFactor};
    rvecVertexNormal[ruIndex] = rVertex3DNormal;
    ruIndex++;
  }

  rvecVertex[ruIndex] = Vertex3D{rvecVertex[ruIntermediate].fX - Vertex3DNormalizedEdge.fX * rfFixedLength * fGlobalFactor / 2.0F - Vertex3DEdgeNormal.fX * rfFixedLength * pow(fGlobalFactor, 2.0F), rvecVertex[ruIntermediate].fY - Vertex3DNormalizedEdge.fY * rfFixedLength * fGlobalFactor / 2.0F - Vertex3DEdgeNormal.fY * rfFixedLength * pow(fGlobalFactor, 2.0F), rvecVertex[ruIntermediate].fZ - Vertex3DNormalizedEdge.fZ * rfFixedLength * fGlobalFactor / 2.0F - Vertex3DEdgeNormal.fZ * rfFixedLength * pow(fGlobalFactor, 2.0F)};
  rvecVertexNormal[ruIndex] = rVertex3DNormal;
  ruIndex++;
  Vertex3DNormalizedEdge = Normalize(Vertex3D{rvecVertex[ruEnd].fX - rvecVertex[ruIntermediate].fX, rvecVertex[ruEnd].fY - rvecVertex[ruIntermediate].fY, rvecVertex[ruEnd].fZ - rvecVertex[ruIntermediate].fZ});
  Vertex3DEdgeNormal = CrossProduct(Vertex3DNormalizedEdge, rVertex3DNormal);
  rvecVertex[ruIndex] = Vertex3D{rvecVertex[ruIntermediate].fX + Vertex3DNormalizedEdge.fX * rfFixedLength * fGlobalFactor / 2.0F - Vertex3DEdgeNormal.fX * rfFixedLength * pow(fGlobalFactor, 2.0F), rvecVertex[ruIntermediate].fY + Vertex3DNormalizedEdge.fY * rfFixedLength * fGlobalFactor / 2.0F - Vertex3DEdgeNormal.fY * rfFixedLength * pow(fGlobalFactor, 2.0F), rvecVertex[ruIntermediate].fZ + Vertex3DNormalizedEdge.fZ * rfFixedLength * fGlobalFactor / 2.0F - Vertex3DEdgeNormal.fZ * rfFixedLength * pow(fGlobalFactor, 2.0F)};
  rvecVertexNormal[ruIndex] = rVertex3DNormal;
  ruIndex++;

  for (int i = 0; i < riInsertionCount; i++) {
    float fLocalFactor = (static_cast<float>(i) + 1.0F) / (static_cast<float>(riInsertionCount) + 1.0F);
    Vertex3D Vertex3DEndpoint{rvecVertex[ruIntermediate].fX * (1.0F - fLocalFactor) + rvecVertex[ruEnd].fX * fLocalFactor, rvecVertex[ruIntermediate].fY * (1.0F - fLocalFactor) + rvecVertex[ruEnd].fY * fLocalFactor, rvecVertex[ruIntermediate].fZ * (1.0F - fLocalFactor) + rvecVertex[ruEnd].fZ * fLocalFactor};

    rvecVertex[ruIndex] = Vertex3D{Vertex3DEndpoint.fX - Vertex3DNormalizedEdge.fX * rfFixedLength * fGlobalFactor - Vertex3DEdgeNormal.fX * rfFixedLength * fGlobalFactor, Vertex3DEndpoint.fY - Vertex3DNormalizedEdge.fY * rfFixedLength * fGlobalFactor - Vertex3DEdgeNormal.fY * rfFixedLength * fGlobalFactor, Vertex3DEndpoint.fZ - Vertex3DNormalizedEdge.fZ * rfFixedLength * fGlobalFactor - Vertex3DEdgeNormal.fZ * rfFixedLength * fGlobalFactor};
    rvecVertexNormal[ruIndex] = rVertex3DNormal;
    ruIndex++;
    rvecVertex[ruIndex] = Vertex3D{Vertex3DEndpoint.fX + Vertex3DNormalizedEdge.fX * rfFixedLength * fGlobalFactor - Vertex3DEdgeNormal.fX * rfFixedLength * fGlobalFactor, Vertex3DEndpoint.fY + Vertex3DNormalizedEdge.fY * rfFixedLength * fGlobalFactor - Vertex3DEdgeNormal.fY * rfFixedLength * fGlobalFactor, Vertex3DEndpoint.fZ + Vertex3DNormalizedEdge.fZ * rfFixedLength * fGlobalFactor - Vertex3DEdgeNormal.fZ * rfFixedLength * fGlobalFactor};
    rvecVertexNormal[ruIndex] = rVertex3DNormal;
    ruIndex++;
  }

  rvecVertex[ruIndex] = Vertex3D{rvecVertex[ruEnd].fX - Vertex3DNormalizedEdge.fX * rfFixedLength * fGlobalFactor / 2.0F - Vertex3DEdgeNormal.fX * rfFixedLength * pow(fGlobalFactor, 2.0F), rvecVertex[ruEnd].fY - Vertex3DNormalizedEdge.fY * rfFixedLength * fGlobalFactor / 2.0F - Vertex3DEdgeNormal.fY * rfFixedLength * pow(fGlobalFactor, 2.0F), rvecVertex[ruEnd].fZ - Vertex3DNormalizedEdge.fZ * rfFixedLength * fGlobalFactor / 2.0F - Vertex3DEdgeNormal.fZ * rfFixedLength * pow(fGlobalFactor, 2.0F)};
  rvecVertexNormal[ruIndex] = rVertex3DNormal;
  ruIndex++;
  Vertex3DNormalizedEdge = Normalize(Vertex3D{rvecVertex[ruStart].fX - rvecVertex[ruEnd].fX, rvecVertex[ruStart].fY - rvecVertex[ruEnd].fY, rvecVertex[ruStart].fZ - rvecVertex[ruEnd].fZ});
  Vertex3DEdgeNormal = CrossProduct(Vertex3DNormalizedEdge, rVertex3DNormal);
  rvecVertex[ruIndex] = Vertex3D{rvecVertex[ruEnd].fX + Vertex3DNormalizedEdge.fX * rfFixedLength * fGlobalFactor / 2.0F - Vertex3DEdgeNormal.fX * rfFixedLength * pow(fGlobalFactor, 2.0F), rvecVertex[ruEnd].fY + Vertex3DNormalizedEdge.fY * rfFixedLength * fGlobalFactor / 2.0F - Vertex3DEdgeNormal.fY * rfFixedLength * pow(fGlobalFactor, 2.0F), rvecVertex[ruEnd].fZ + Vertex3DNormalizedEdge.fZ * rfFixedLength * fGlobalFactor / 2.0F - Vertex3DEdgeNormal.fZ * rfFixedLength * pow(fGlobalFactor, 2.0F)};
  rvecVertexNormal[ruIndex] = rVertex3DNormal;
  ruIndex++;

  for (int i = 0; i < riInsertionCount; i++) {
    float fLocalFactor = (static_cast<float>(i) + 1.0F) / (static_cast<float>(riInsertionCount) + 1.0F);
    Vertex3D Vertex3DEndpoint{rvecVertex[ruEnd].fX * (1.0F - fLocalFactor) + rvecVertex[ruStart].fX * fLocalFactor, rvecVertex[ruEnd].fY * (1.0F - fLocalFactor) + rvecVertex[ruStart].fY * fLocalFactor, rvecVertex[ruEnd].fZ * (1.0F - fLocalFactor) + rvecVertex[ruStart].fZ * fLocalFactor};

    rvecVertex[ruIndex] = Vertex3D{Vertex3DEndpoint.fX - Vertex3DNormalizedEdge.fX * rfFixedLength * fGlobalFactor - Vertex3DEdgeNormal.fX * rfFixedLength * fGlobalFactor, Vertex3DEndpoint.fY - Vertex3DNormalizedEdge.fY * rfFixedLength * fGlobalFactor - Vertex3DEdgeNormal.fY * rfFixedLength * fGlobalFactor, Vertex3DEndpoint.fZ - Vertex3DNormalizedEdge.fZ * rfFixedLength * fGlobalFactor - Vertex3DEdgeNormal.fZ * rfFixedLength * fGlobalFactor};
    rvecVertexNormal[ruIndex] = rVertex3DNormal;
    ruIndex++;
    rvecVertex[ruIndex] = Vertex3D{Vertex3DEndpoint.fX + Vertex3DNormalizedEdge.fX * rfFixedLength * fGlobalFactor - Vertex3DEdgeNormal.fX * rfFixedLength * fGlobalFactor, Vertex3DEndpoint.fY + Vertex3DNormalizedEdge.fY * rfFixedLength * fGlobalFactor - Vertex3DEdgeNormal.fY * rfFixedLength * fGlobalFactor, Vertex3DEndpoint.fZ + Vertex3DNormalizedEdge.fZ * rfFixedLength * fGlobalFactor - Vertex3DEdgeNormal.fZ * rfFixedLength * fGlobalFactor};
    rvecVertexNormal[ruIndex] = rVertex3DNormal;
    ruIndex++;
  }

  rvecVertex[ruIndex] = Vertex3D{rvecVertex[ruStart].fX - Vertex3DNormalizedEdge.fX * rfFixedLength * fGlobalFactor / 2.0F - Vertex3DEdgeNormal.fX * rfFixedLength * pow(fGlobalFactor, 2.0F), rvecVertex[ruStart].fY - Vertex3DNormalizedEdge.fY * rfFixedLength * fGlobalFactor / 2.0F - Vertex3DEdgeNormal.fY * rfFixedLength * pow(fGlobalFactor, 2.0F), rvecVertex[ruStart].fZ - Vertex3DNormalizedEdge.fZ * rfFixedLength * fGlobalFactor / 2.0F - Vertex3DEdgeNormal.fZ * rfFixedLength * pow(fGlobalFactor, 2.0F)};
  rvecVertexNormal[ruIndex] = rVertex3DNormal;
  ruIndex++;
}

void FillUpTriangle(vector<Vertex3D>& rvecVertex, vector<Vertex3D>& rvecVertexNormal, const unsigned& ruStart, const unsigned& ruIntermediate, const unsigned& ruEnd, const unsigned& ruMark, unsigned& ruIndex, const Vertex3D& rVertex3DNormal, const int& riInsertionCount)
{
  rvecVertex[ruMark] = Vertex3D{(rvecVertex[ruStart].fX + rvecVertex[ruIntermediate].fX + rvecVertex[ruEnd].fX) / static_cast<float>(TRIANGLE_ELEMENT_COUNT), (rvecVertex[ruStart].fY + rvecVertex[ruIntermediate].fY + rvecVertex[ruEnd].fY) / static_cast<float>(TRIANGLE_ELEMENT_COUNT), (rvecVertex[ruStart].fZ + rvecVertex[ruIntermediate].fZ + rvecVertex[ruEnd].fZ) / static_cast<float>(TRIANGLE_ELEMENT_COUNT)};
  rvecVertexNormal[ruMark] = rVertex3DNormal;
  ruIndex++;

  if (riInsertionCount - 1) {
    FillUpTriangle(rvecVertex, rvecVertexNormal, ruStart, ruIntermediate, ruMark, unsigned(ruIndex), ruIndex, rVertex3DNormal, riInsertionCount - 1);
    FillUpTriangle(rvecVertex, rvecVertexNormal, ruIntermediate, ruEnd, ruMark, unsigned(ruIndex), ruIndex, rVertex3DNormal, riInsertionCount - 1);
    FillUpTriangle(rvecVertex, rvecVertexNormal, ruStart, ruEnd, ruMark, unsigned(ruIndex), ruIndex, rVertex3DNormal, riInsertionCount - 1);
  }
}

void ConfirmSolution(const float* const* const& rpfMatrix, const vector<float>& rvecCoefficient, const size_t& rszSize, const DisplayMode& rDisplayModeSetting)
{
  float* pfResult = new float[rszSize]();

  for (int i = 0; i < rszSize; i++)
    for (int j = 0; j < rszSize; j++)
      pfResult[i] += rpfMatrix[i][j] * rvecCoefficient[j];

  ShowElement(pfResult, rszSize, "pfResult", rDisplayModeSetting);

  delete[] pfResult;
}