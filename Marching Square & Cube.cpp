#include "Definition.h"
#include "Geometric Operation.h"
#include "Marching Square & Cube.h"

using namespace std;

void MarchingSquare(ofstream& rofsFile, const PiecewisePolynomial2D& rPiecewisePolynomial2DData, const float (*const& rpfTransformationMatrix)[ELEMENT_COUNT_2D + 1], const int& riDivisionCountX, const int& riDivisionCountY)
{
  static const vector<vector<unsigned char>> vecTrianglePattern{
    vector<unsigned char>{},
    vector<unsigned char>{0U, 4U, 7U},
    vector<unsigned char>{1U, 5U, 4U},
    vector<unsigned char>{0U, 1U, 7U, 1U, 5U, 7U},
    vector<unsigned char>{2U, 6U, 5U},
    vector<unsigned char>{0U, 4U, 7U, 2U, 6U, 5U},
    vector<unsigned char>{1U, 2U, 4U, 2U, 6U, 4U},
    vector<unsigned char>{0U, 1U, 7U, 1U, 6U, 7U, 1U, 2U, 6U},
    vector<unsigned char>{3U, 7U, 6U},
    vector<unsigned char>{0U, 4U, 3U, 3U, 4U, 6U},
    vector<unsigned char>{1U, 5U, 4U, 3U, 7U, 6U},
    vector<unsigned char>{0U, 1U, 5U, 0U, 5U, 6U, 0U, 6U, 3U},
    vector<unsigned char>{2U, 3U, 5U, 3U, 7U, 5U},
    vector<unsigned char>{0U, 4U, 3U, 3U, 4U, 5U, 2U, 3U, 5U},
    vector<unsigned char>{1U, 2U, 4U, 2U, 7U, 4U, 2U, 3U, 7U},
    vector<unsigned char>{0U, 1U, 2U, 0U, 2U, 3U}
  }, vecSpecialPattern{
    vector<unsigned char>{0U, 4U, 7U, 2U, 6U, 5U, 4U, 5U, 7U, 6U, 7U, 5U},
    vector<unsigned char>{1U, 5U, 4U, 3U, 7U, 6U, 4U, 5U, 7U, 6U, 7U, 5U}
  };
  static const unsigned char uchEdgePattern[1U << GRID_VERTEX_COUNT_2D] = {0x0U, 0x91U, 0x32U, 0xA3U, 0x64U, 0xF5U, 0x56U, 0xC7U, 0xC8U, 0x59U, 0xFAU, 0x6BU, 0xACU, 0x3DU, 0x9EU, 0xFU};

  vector<Vertex2D> vecTransformedCorner{rPiecewisePolynomial2DData.MinimumCorner(), rPiecewisePolynomial2DData.MaximumCorner()};
  float** pfGridValueMatrix = new float*[riDivisionCountX + 1];
  float fInversedTransformationMatrix[ELEMENT_COUNT_2D + 1][ELEMENT_COUNT_2D + 1];
  int iTriangleCount = 0;

  for (int i = 0; i <= riDivisionCountX; i++)
    pfGridValueMatrix[i] = new float[riDivisionCountY + 1];

  if (rpfTransformationMatrix == nullptr)
    ;
  else {
    for (int i = 0; i <= ELEMENT_COUNT_2D; i++)
      for (int j = 0; j <= ELEMENT_COUNT_2D; j++)
        fInversedTransformationMatrix[i][j] = rpfTransformationMatrix[i][j];

    InverseTransformation<ELEMENT_COUNT_2D + 1>(fInversedTransformationMatrix);

    float Vertex2D::*pfVertex2DMember[ELEMENT_COUNT_2D] = {&Vertex2D::fX, &Vertex2D::fY};
    vector<Vertex2D> vecOriginalCorner{rPiecewisePolynomial2DData.MinimumCorner(), rPiecewisePolynomial2DData.MaximumCorner()};

    for (int i = 0; i < ELEMENT_COUNT_2D; i++) {
      vecTransformedCorner[INDICATOR_MINIMUM].*pfVertex2DMember[i] = vecTransformedCorner[INDICATOR_MAXIMUM].*pfVertex2DMember[i] = 0.0F;

      for (int j = 0; j < ELEMENT_COUNT_2D; j++)
        if (rpfTransformationMatrix[i][j] * vecOriginalCorner[INDICATOR_MINIMUM].*pfVertex2DMember[j] < rpfTransformationMatrix[i][j] * vecOriginalCorner[INDICATOR_MAXIMUM].*pfVertex2DMember[j]) {
          vecTransformedCorner[INDICATOR_MINIMUM].*pfVertex2DMember[i] += rpfTransformationMatrix[i][j] * vecOriginalCorner[INDICATOR_MINIMUM].*pfVertex2DMember[j];
          vecTransformedCorner[INDICATOR_MAXIMUM].*pfVertex2DMember[i] += rpfTransformationMatrix[i][j] * vecOriginalCorner[INDICATOR_MAXIMUM].*pfVertex2DMember[j];
        } else {
          vecTransformedCorner[INDICATOR_MINIMUM].*pfVertex2DMember[i] += rpfTransformationMatrix[i][j] * vecOriginalCorner[INDICATOR_MAXIMUM].*pfVertex2DMember[j];
          vecTransformedCorner[INDICATOR_MAXIMUM].*pfVertex2DMember[i] += rpfTransformationMatrix[i][j] * vecOriginalCorner[INDICATOR_MINIMUM].*pfVertex2DMember[j];
        }

      vecTransformedCorner[INDICATOR_MINIMUM].*pfVertex2DMember[i] += rpfTransformationMatrix[i][ELEMENT_COUNT_2D];
      vecTransformedCorner[INDICATOR_MAXIMUM].*pfVertex2DMember[i] += rpfTransformationMatrix[i][ELEMENT_COUNT_2D];
    }
  }

  float fStep[ELEMENT_COUNT_2D] = {(vecTransformedCorner[INDICATOR_MAXIMUM].fX - vecTransformedCorner[INDICATOR_MINIMUM].fX) / static_cast<float>(riDivisionCountX), (vecTransformedCorner[INDICATOR_MAXIMUM].fY - vecTransformedCorner[INDICATOR_MINIMUM].fY) / static_cast<float>(riDivisionCountY)};

  for (int i = 0; i <= riDivisionCountX; i++)
    for (int j = 0; j <= riDivisionCountY; j++) {
      Vertex2D Vertex2DSamplingPoint{vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * static_cast<float>(i), vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * static_cast<float>(j)};

      pfGridValueMatrix[i][j] = rPiecewisePolynomial2DData.At(rpfTransformationMatrix == nullptr ? Vertex2DSamplingPoint : TransformCoordinate(fInversedTransformationMatrix, Vertex2DSamplingPoint));
    }

  unsigned char uchOffset[GRID_VERTEX_COUNT_2D][ELEMENT_COUNT_2D] = {
    {0U, 0U},
    {0U, 1U},
    {1U, 1U},
    {1U, 0U}
  };

  for (int i = 0; i < riDivisionCountX; i++)
    for (int j = 0; j < riDivisionCountY; j++) {
      unsigned char uchFlag = 0U;

      do {
        unsigned char uchDigit = 1U;

        for (const unsigned char *const& rpuchPointer : uchOffset) {
          if (pfGridValueMatrix[rpuchPointer[INDICATOR_X] + i][rpuchPointer[INDICATOR_Y] + j] >= 0.0F)
            uchFlag |= uchDigit;

          uchDigit <<= 1U;
        }
      } while (false);

      if (!uchEdgePattern[uchFlag])
        continue;

      vector<Vertex2D> vecInterpolationPoint(GRID_VERTEX_COUNT_2D + GRID_EDGE_COUNT_2D, Vertex2D{0.0F, 0.0F});

      do {
        unsigned char uchFactor = 0U;

        if (uchEdgePattern[uchFlag] & 1U << uchFactor++)
          vecInterpolationPoint[uchFactor - 1] = Vertex2D{vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * static_cast<float>(i), vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * static_cast<float>(j)};

        if (uchEdgePattern[uchFlag] & 1U << uchFactor++)
          vecInterpolationPoint[uchFactor - 1] = Vertex2D{vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * static_cast<float>(i), vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * (static_cast<float>(j) + 1.0F)};

        if (uchEdgePattern[uchFlag] & 1U << uchFactor++)
          vecInterpolationPoint[uchFactor - 1] = Vertex2D{vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * (static_cast<float>(i) + 1.0F), vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * (static_cast<float>(j) + 1.0F)};

        if (uchEdgePattern[uchFlag] & 1U << uchFactor++)
          vecInterpolationPoint[uchFactor - 1] = Vertex2D{vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * (static_cast<float>(i) + 1.0F), vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * static_cast<float>(j)};

        if (uchEdgePattern[uchFlag] & 1U << uchFactor++)
          vecInterpolationPoint[uchFactor - 1] = Vertex2D{vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * static_cast<float>(i), (pfGridValueMatrix[i][j + 1] * (vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * static_cast<float>(j)) - pfGridValueMatrix[i][j] * (vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * (static_cast<float>(j) + 1.0F))) / (pfGridValueMatrix[i][j + 1] - pfGridValueMatrix[i][j])};

        if (uchEdgePattern[uchFlag] & 1U << uchFactor++)
          vecInterpolationPoint[uchFactor - 1] = Vertex2D{(pfGridValueMatrix[i + 1][j + 1] * (vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * static_cast<float>(i)) - pfGridValueMatrix[i][j + 1] * (vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * (static_cast<float>(i) + 1.0F))) / (pfGridValueMatrix[i + 1][j + 1] - pfGridValueMatrix[i][j + 1]), vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * (static_cast<float>(j) + 1.0F)};

        if (uchEdgePattern[uchFlag] & 1U << uchFactor++)
          vecInterpolationPoint[uchFactor - 1] = Vertex2D{vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * (static_cast<float>(i) + 1.0F), (pfGridValueMatrix[i + 1][j + 1] * (vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * static_cast<float>(j)) - pfGridValueMatrix[i + 1][j] * (vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * (static_cast<float>(j) + 1.0F))) / (pfGridValueMatrix[i + 1][j + 1] - pfGridValueMatrix[i + 1][j])};

        if (uchEdgePattern[uchFlag] & 1U << uchFactor++)
          vecInterpolationPoint[uchFactor - 1] = Vertex2D{(pfGridValueMatrix[i + 1][j] * (vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * static_cast<float>(i)) - pfGridValueMatrix[i][j] * (vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * (static_cast<float>(i) + 1.0F))) / (pfGridValueMatrix[i + 1][j] - pfGridValueMatrix[i][j]), vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * static_cast<float>(j)};
      } while (false);

      Vertex2D Vertex2DCenter{vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * (static_cast<float>(i) + 0.5F), vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * (static_cast<float>(j) + 0.5F)};

      if (uchFlag == SPECIAL_PATTERN_A && rPiecewisePolynomial2DData.At(rpfTransformationMatrix == nullptr ? Vertex2DCenter : TransformCoordinate(fInversedTransformationMatrix, Vertex2DCenter)) >= 0.0F)
        for (int k = 0; k * TRIANGLE_ELEMENT_COUNT < vecSpecialPattern.front().size(); k++) {
          for (int l = 0; l < TRIANGLE_ELEMENT_COUNT; l++)
            rofsFile << VERTEX_CHARACTER << ' ' << vecInterpolationPoint[vecSpecialPattern.front()[k * TRIANGLE_ELEMENT_COUNT + l]].fX << ' ' << vecInterpolationPoint[vecSpecialPattern.front()[k * TRIANGLE_ELEMENT_COUNT + l]].fY << ' ' << 0.0F << endl;

          iTriangleCount++;
        }
      else if (uchFlag == SPECIAL_PATTERN_B && rPiecewisePolynomial2DData.At(rpfTransformationMatrix == nullptr ? Vertex2DCenter : TransformCoordinate(fInversedTransformationMatrix, Vertex2DCenter)) >= 0.0F)
        for (int k = 0; k * TRIANGLE_ELEMENT_COUNT < vecSpecialPattern.back().size(); k++) {
          for (int l = 0; l < TRIANGLE_ELEMENT_COUNT; l++)
            rofsFile << VERTEX_CHARACTER << ' ' << vecInterpolationPoint[vecSpecialPattern.back()[k * TRIANGLE_ELEMENT_COUNT + l]].fX << ' ' << vecInterpolationPoint[vecSpecialPattern.back()[k * TRIANGLE_ELEMENT_COUNT + l]].fY << ' ' << 0.0F << endl;

          iTriangleCount++;
        }
      else
        for (int k = 0; k * TRIANGLE_ELEMENT_COUNT < vecTrianglePattern[uchFlag].size(); k++) {
          for (int l = 0; l < TRIANGLE_ELEMENT_COUNT; l++)
            rofsFile << VERTEX_CHARACTER << ' ' << vecInterpolationPoint[vecTrianglePattern[uchFlag][k * TRIANGLE_ELEMENT_COUNT + l]].fX << ' ' << vecInterpolationPoint[vecTrianglePattern[uchFlag][k * TRIANGLE_ELEMENT_COUNT + l]].fY << ' ' << 0.0F << endl;

          iTriangleCount++;
        }
    }

  int iVertexCount = 0;

  for (int i = 0; i < iTriangleCount; i++) {
    rofsFile << FACE_ELEMENT_CHARACTER << ' ';

    for (int j = 0; j < TRIANGLE_ELEMENT_COUNT; j++) {
      rofsFile << iVertexCount + j + 1 << DELIMITER_CHARACTER << DELIMITER_CHARACTER;

      if (TRIANGLE_ELEMENT_COUNT - j - 1)
        rofsFile << ' ';
      else
        rofsFile << endl;
    }

    iVertexCount += TRIANGLE_ELEMENT_COUNT;
  }
}

void MarchingSquare(ofstream& rofsFile, const vector<PiecewisePolynomial2D>& rvecImplicitSurface, const float (*const& rpfTransformationMatrix)[ELEMENT_COUNT_2D + 1][ELEMENT_COUNT_2D + 1], const int& riDivisionCountX, const int& riDivisionCountY)
{
  static const vector<vector<unsigned char>> vecTrianglePattern{
    vector<unsigned char>{},
    vector<unsigned char>{0U, 4U, 7U},
    vector<unsigned char>{1U, 5U, 4U},
    vector<unsigned char>{0U, 1U, 7U, 1U, 5U, 7U},
    vector<unsigned char>{2U, 6U, 5U},
    vector<unsigned char>{0U, 4U, 7U, 2U, 6U, 5U},
    vector<unsigned char>{1U, 2U, 4U, 2U, 6U, 4U},
    vector<unsigned char>{0U, 1U, 7U, 1U, 6U, 7U, 1U, 2U, 6U},
    vector<unsigned char>{3U, 7U, 6U},
    vector<unsigned char>{0U, 4U, 3U, 3U, 4U, 6U},
    vector<unsigned char>{1U, 5U, 4U, 3U, 7U, 6U},
    vector<unsigned char>{0U, 1U, 5U, 0U, 5U, 6U, 0U, 6U, 3U},
    vector<unsigned char>{2U, 3U, 5U, 3U, 7U, 5U},
    vector<unsigned char>{0U, 4U, 3U, 3U, 4U, 5U, 2U, 3U, 5U},
    vector<unsigned char>{1U, 2U, 4U, 2U, 7U, 4U, 2U, 3U, 7U},
    vector<unsigned char>{0U, 1U, 2U, 0U, 2U, 3U}
  }, vecSpecialPattern{
    vector<unsigned char>{0U, 4U, 7U, 2U, 6U, 5U, 4U, 5U, 7U, 6U, 7U, 5U},
    vector<unsigned char>{1U, 5U, 4U, 3U, 7U, 6U, 4U, 5U, 7U, 6U, 7U, 5U}
  };
  static const unsigned char uchEdgePattern[1U << GRID_VERTEX_COUNT_2D] = {0x0U, 0x91U, 0x32U, 0xA3U, 0x64U, 0xF5U, 0x56U, 0xC7U, 0xC8U, 0x59U, 0xFAU, 0x6BU, 0xACU, 0x3DU, 0x9EU, 0xFU};

  float** pfGridValueMatrix = new float*[riDivisionCountX + 1];
  int* piTriangleCountArray = new int[rvecImplicitSurface.size()]();

  for (int i = 0; i <= riDivisionCountX; i++)
    pfGridValueMatrix[i] = new float[riDivisionCountY + 1];

  for (int i = 0; i < rvecImplicitSurface.size(); i++) {
    vector<Vertex2D> vecTransformedCorner{rvecImplicitSurface[i].MinimumCorner(), rvecImplicitSurface[i].MaximumCorner()};
    float fInversedTransformationMatrix[ELEMENT_COUNT_2D + 1][ELEMENT_COUNT_2D + 1];

    for (int j = 0; j <= ELEMENT_COUNT_2D; j++)
      for (int k = 0; k <= ELEMENT_COUNT_2D; k++)
        fInversedTransformationMatrix[j][k] = rpfTransformationMatrix[i][j][k];

    InverseTransformation<ELEMENT_COUNT_2D + 1>(fInversedTransformationMatrix);

    float Vertex2D::*pfVertex2DMember[ELEMENT_COUNT_2D] = {&Vertex2D::fX, &Vertex2D::fY};
    vector<Vertex2D> vecOriginalCorner{rvecImplicitSurface[i].MinimumCorner(), rvecImplicitSurface[i].MaximumCorner()};

    for (int j = 0; j < ELEMENT_COUNT_2D; j++) {
      vecTransformedCorner[INDICATOR_MINIMUM].*pfVertex2DMember[j] = vecTransformedCorner[INDICATOR_MAXIMUM].*pfVertex2DMember[j] = 0.0F;

      for (int k = 0; k < ELEMENT_COUNT_2D; k++)
        if (rpfTransformationMatrix[i][j][k] * vecOriginalCorner[INDICATOR_MINIMUM].*pfVertex2DMember[k] < rpfTransformationMatrix[i][j][k] * vecOriginalCorner[INDICATOR_MAXIMUM].*pfVertex2DMember[k]) {
          vecTransformedCorner[INDICATOR_MINIMUM].*pfVertex2DMember[j] += rpfTransformationMatrix[i][j][k] * vecOriginalCorner[INDICATOR_MINIMUM].*pfVertex2DMember[k];
          vecTransformedCorner[INDICATOR_MAXIMUM].*pfVertex2DMember[j] += rpfTransformationMatrix[i][j][k] * vecOriginalCorner[INDICATOR_MAXIMUM].*pfVertex2DMember[k];
        } else {
          vecTransformedCorner[INDICATOR_MINIMUM].*pfVertex2DMember[j] += rpfTransformationMatrix[i][j][k] * vecOriginalCorner[INDICATOR_MAXIMUM].*pfVertex2DMember[k];
          vecTransformedCorner[INDICATOR_MAXIMUM].*pfVertex2DMember[j] += rpfTransformationMatrix[i][j][k] * vecOriginalCorner[INDICATOR_MINIMUM].*pfVertex2DMember[k];
        }

      vecTransformedCorner[INDICATOR_MINIMUM].*pfVertex2DMember[j] += rpfTransformationMatrix[i][j][ELEMENT_COUNT_2D];
      vecTransformedCorner[INDICATOR_MAXIMUM].*pfVertex2DMember[j] += rpfTransformationMatrix[i][j][ELEMENT_COUNT_2D];
    }

    float fStep[ELEMENT_COUNT_2D] = {(vecTransformedCorner[INDICATOR_MAXIMUM].fX - vecTransformedCorner[INDICATOR_MINIMUM].fX) / static_cast<float>(riDivisionCountX), (vecTransformedCorner[INDICATOR_MAXIMUM].fY - vecTransformedCorner[INDICATOR_MINIMUM].fY) / static_cast<float>(riDivisionCountY)};

    for (int j = 0; j <= riDivisionCountX; j++)
      for (int k = 0; k <= riDivisionCountY; k++) {
        Vertex2D Vertex2DSamplingPoint{vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * static_cast<float>(j), vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * static_cast<float>(k)};

        pfGridValueMatrix[j][k] = rvecImplicitSurface[i].At(TransformCoordinate(fInversedTransformationMatrix, Vertex2DSamplingPoint));
      }

    unsigned char uchOffset[GRID_VERTEX_COUNT_2D][ELEMENT_COUNT_2D] = {
      {0U, 0U},
      {0U, 1U},
      {1U, 1U},
      {1U, 0U}
    };

    for (int j = 0; j < riDivisionCountX; j++)
      for (int k = 0; k < riDivisionCountY; k++) {
        unsigned char uchFlag = 0U;

        do {
          unsigned char uchDigit = 1U;

          for (const unsigned char* const& rpuchPointer : uchOffset) {
            if (pfGridValueMatrix[rpuchPointer[INDICATOR_X] + i][rpuchPointer[INDICATOR_Y] + j] >= 0.0F)
              uchFlag |= uchDigit;

            uchDigit <<= 1U;
          }
        } while (false);

        if (!uchEdgePattern[uchFlag])
          continue;

        vector<Vertex2D> vecInterpolationPoint(GRID_VERTEX_COUNT_2D + GRID_EDGE_COUNT_2D, Vertex2D{0.0F, 0.0F});

        do {
          unsigned char uchFactor = 0U;

          if (uchEdgePattern[uchFlag] & 1U << uchFactor++)
            vecInterpolationPoint[uchFactor - 1] = Vertex2D{vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * static_cast<float>(i), vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * static_cast<float>(j)};

          if (uchEdgePattern[uchFlag] & 1U << uchFactor++)
            vecInterpolationPoint[uchFactor - 1] = Vertex2D{vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * static_cast<float>(i), vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * (static_cast<float>(j) + 1.0F)};

          if (uchEdgePattern[uchFlag] & 1U << uchFactor++)
            vecInterpolationPoint[uchFactor - 1] = Vertex2D{vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * (static_cast<float>(i) + 1.0F), vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * (static_cast<float>(j) + 1.0F)};

          if (uchEdgePattern[uchFlag] & 1U << uchFactor++)
            vecInterpolationPoint[uchFactor - 1] = Vertex2D{vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * (static_cast<float>(i) + 1.0F), vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * static_cast<float>(j)};

          if (uchEdgePattern[uchFlag] & 1U << uchFactor++)
            vecInterpolationPoint[uchFactor - 1] = Vertex2D{vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * static_cast<float>(i), (pfGridValueMatrix[i][j + 1] * (vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * static_cast<float>(j)) - pfGridValueMatrix[i][j] * (vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * (static_cast<float>(j) + 1.0F))) / (pfGridValueMatrix[i][j + 1] - pfGridValueMatrix[i][j])};

          if (uchEdgePattern[uchFlag] & 1U << uchFactor++)
            vecInterpolationPoint[uchFactor - 1] = Vertex2D{(pfGridValueMatrix[i + 1][j + 1] * (vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * static_cast<float>(i)) - pfGridValueMatrix[i][j + 1] * (vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * (static_cast<float>(i) + 1.0F))) / (pfGridValueMatrix[i + 1][j + 1] - pfGridValueMatrix[i][j + 1]), vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * (static_cast<float>(j) + 1.0F)};

          if (uchEdgePattern[uchFlag] & 1U << uchFactor++)
            vecInterpolationPoint[uchFactor - 1] = Vertex2D{vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * (static_cast<float>(i) + 1.0F), (pfGridValueMatrix[i + 1][j + 1] * (vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * static_cast<float>(j)) - pfGridValueMatrix[i + 1][j] * (vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * (static_cast<float>(j) + 1.0F))) / (pfGridValueMatrix[i + 1][j + 1] - pfGridValueMatrix[i + 1][j])};

          if (uchEdgePattern[uchFlag] & 1U << uchFactor++)
            vecInterpolationPoint[uchFactor - 1] = Vertex2D{(pfGridValueMatrix[i + 1][j] * (vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * static_cast<float>(i)) - pfGridValueMatrix[i][j] * (vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * (static_cast<float>(i) + 1.0F))) / (pfGridValueMatrix[i + 1][j] - pfGridValueMatrix[i][j]), vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * static_cast<float>(j)};
        } while (false);

        Vertex2D Vertex2DCenter{vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * (static_cast<float>(j) + 0.5F), vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * (static_cast<float>(k) + 0.5F)};

        if (uchFlag == SPECIAL_PATTERN_A && rvecImplicitSurface[i].At(TransformCoordinate(fInversedTransformationMatrix, Vertex2DCenter)) >= 0.0F)
          for (int l = 0; l * TRIANGLE_ELEMENT_COUNT < vecSpecialPattern.front().size(); l++) {
            for (int m = 0; m < TRIANGLE_ELEMENT_COUNT; m++)
              rofsFile << VERTEX_CHARACTER << ' ' << vecInterpolationPoint[vecSpecialPattern.front()[l * TRIANGLE_ELEMENT_COUNT + m]].fX << ' ' << vecInterpolationPoint[vecSpecialPattern.front()[l * TRIANGLE_ELEMENT_COUNT + m]].fY << ' ' << 0.0F << endl;

            piTriangleCountArray[i]++;
          }
        else if (uchFlag == SPECIAL_PATTERN_B && rvecImplicitSurface[i].At(TransformCoordinate(fInversedTransformationMatrix, Vertex2DCenter)) >= 0.0F)
          for (int l = 0; l * TRIANGLE_ELEMENT_COUNT < vecSpecialPattern.back().size(); l++) {
            for (int m = 0; m < TRIANGLE_ELEMENT_COUNT; m++)
              rofsFile << VERTEX_CHARACTER << ' ' << vecInterpolationPoint[vecSpecialPattern.back()[l * TRIANGLE_ELEMENT_COUNT + m]].fX << ' ' << vecInterpolationPoint[vecSpecialPattern.back()[l * TRIANGLE_ELEMENT_COUNT + m]].fY << ' ' << 0.0F << endl;

            piTriangleCountArray[i]++;
          }
        else
          for (int l = 0; l * TRIANGLE_ELEMENT_COUNT < vecTrianglePattern[uchFlag].size(); l++) {
            for (int m = 0; m < TRIANGLE_ELEMENT_COUNT; m++)
              rofsFile << VERTEX_CHARACTER << ' ' << vecInterpolationPoint[vecTrianglePattern[uchFlag][l * TRIANGLE_ELEMENT_COUNT + m]].fX << ' ' << vecInterpolationPoint[vecTrianglePattern[uchFlag][l * TRIANGLE_ELEMENT_COUNT + m]].fY << ' ' << 0.0F << endl;

            piTriangleCountArray[i]++;
          }
      }

    if (rvecImplicitSurface.size() - i - 1)
      rofsFile << COMMENT_CHARACTER << endl;
  }

  for (int i = 0; i <= riDivisionCountX; i++)
    delete[] pfGridValueMatrix[i];

  delete[] pfGridValueMatrix;

  int iVertexCount = 0;

  for (int i = 0; i < rvecImplicitSurface.size(); i++) {
    for (int j = 0; j < piTriangleCountArray[i]; j++) {
      rofsFile << FACE_ELEMENT_CHARACTER << ' ';

      for (int k = 0; k < TRIANGLE_ELEMENT_COUNT; k++) {
        rofsFile << iVertexCount + k + 1 << DELIMITER_CHARACTER << DELIMITER_CHARACTER;

        if (TRIANGLE_ELEMENT_COUNT - k - 1)
          rofsFile << ' ';
        else
          rofsFile << endl;
      }

      iVertexCount += TRIANGLE_ELEMENT_COUNT;
    }

    if (rvecImplicitSurface.size() - i - 1)
      rofsFile << COMMENT_CHARACTER << endl;
  }

  delete[] piTriangleCountArray;
}

void MarchingCube(ofstream& rofsFile, const PiecewisePolynomial3D& rPiecewisePolynomial3DData, const float (*const& rpfTransformationMatrix)[ELEMENT_COUNT_3D + 1], const int& riDivisionCountX, const int& riDivisionCountY, const int& riDivisionCountZ)
{
  static const vector<vector<unsigned char>> vecTrianglePattern{
    vector<unsigned char>{},
    vector<unsigned char>{0U, 8U, 3U},
    vector<unsigned char>{0U, 1U, 9U},
    vector<unsigned char>{1U, 8U, 3U, 9U, 8U, 1U},
    vector<unsigned char>{1U, 2U, 10U},
    vector<unsigned char>{0U, 8U, 3U, 1U, 2U, 10U},
    vector<unsigned char>{9U, 2U, 10U, 0U, 2U, 9U},
    vector<unsigned char>{2U, 8U, 3U, 2U, 10U, 8U, 10U, 9U, 8U},
    vector<unsigned char>{3U, 11U, 2U},
    vector<unsigned char>{0U, 11U, 2U, 8U, 11U, 0U},
    vector<unsigned char>{1U, 9U, 0U, 2U, 3U, 11U},
    vector<unsigned char>{1U, 11U, 2U, 1U, 9U, 11U, 9U, 8U, 11U},
    vector<unsigned char>{3U, 10U, 1U, 11U, 10U, 3U},
    vector<unsigned char>{0U, 10U, 1U, 0U, 8U, 10U, 8U, 11U, 10U},
    vector<unsigned char>{3U, 9U, 0U, 3U, 11U, 9U, 11U, 10U, 9U},
    vector<unsigned char>{9U, 8U, 10U, 10U, 8U, 11U},
    vector<unsigned char>{4U, 7U, 8U},
    vector<unsigned char>{4U, 3U, 0U, 7U, 3U, 4U},
    vector<unsigned char>{0U, 1U, 9U, 8U, 4U, 7U},
    vector<unsigned char>{4U, 1U, 9U, 4U, 7U, 1U, 7U, 3U, 1U},
    vector<unsigned char>{1U, 2U, 10U, 8U, 4U, 7U},
    vector<unsigned char>{3U, 4U, 7U, 3U, 0U, 4U, 1U, 2U, 10U},
    vector<unsigned char>{9U, 2U, 10U, 9U, 0U, 2U, 8U, 4U, 7U},
    vector<unsigned char>{2U, 10U, 9U, 2U, 9U, 7U, 2U, 7U, 3U, 7U, 9U, 4U},
    vector<unsigned char>{8U, 4U, 7U, 3U, 11U, 2U},
    vector<unsigned char>{11U, 4U, 7U, 11U, 2U, 4U, 2U, 0U, 4U},
    vector<unsigned char>{9U, 0U, 1U, 8U, 4U, 7U, 2U, 3U, 11U},
    vector<unsigned char>{4U, 7U, 11U, 9U, 4U, 11U, 9U, 11U, 2U, 9U, 2U, 1U},
    vector<unsigned char>{3U, 10U, 1U, 3U, 11U, 10U, 7U, 8U, 4U},
    vector<unsigned char>{1U, 11U, 10U, 1U, 4U, 11U, 1U, 0U, 4U, 7U, 11U, 4U},
    vector<unsigned char>{4U, 7U, 8U, 9U, 0U, 11U, 9U, 11U, 10U, 11U, 0U, 3U},
    vector<unsigned char>{4U, 7U, 11U, 4U, 11U, 9U, 9U, 11U, 10U},
    vector<unsigned char>{9U, 5U, 4U},
    vector<unsigned char>{9U, 5U, 4U, 0U, 8U, 3U},
    vector<unsigned char>{0U, 5U, 4U, 1U, 5U, 0U},
    vector<unsigned char>{8U, 5U, 4U, 8U, 3U, 5U, 3U, 1U, 5U},
    vector<unsigned char>{1U, 2U, 10U, 9U, 5U, 4U},
    vector<unsigned char>{3U, 0U, 8U, 1U, 2U, 10U, 4U, 9U, 5U},
    vector<unsigned char>{5U, 2U, 10U, 5U, 4U, 2U, 4U, 0U, 2U},
    vector<unsigned char>{2U, 10U, 5U, 3U, 2U, 5U, 3U, 5U, 4U, 3U, 4U, 8U},
    vector<unsigned char>{9U, 5U, 4U, 2U, 3U, 11U},
    vector<unsigned char>{0U, 11U, 2U, 0U, 8U, 11U, 4U, 9U, 5U},
    vector<unsigned char>{0U, 5U, 4U, 0U, 1U, 5U, 2U, 3U, 11U},
    vector<unsigned char>{2U, 1U, 5U, 2U, 5U, 8U, 2U, 8U, 11U, 4U, 8U, 5U},
    vector<unsigned char>{10U, 3U, 11U, 10U, 1U, 3U, 9U, 5U, 4U},
    vector<unsigned char>{4U, 9U, 5U, 0U, 8U, 1U, 8U, 10U, 1U, 8U, 11U, 10U},
    vector<unsigned char>{5U, 4U, 0U, 5U, 0U, 11U, 5U, 11U, 10U, 11U, 0U, 3U},
    vector<unsigned char>{5U, 4U, 8U, 5U, 8U, 10U, 10U, 8U, 11U},
    vector<unsigned char>{9U, 7U, 8U, 5U, 7U, 9U},
    vector<unsigned char>{9U, 3U, 0U, 9U, 5U, 3U, 5U, 7U, 3U},
    vector<unsigned char>{0U, 7U, 8U, 0U, 1U, 7U, 1U, 5U, 7U},
    vector<unsigned char>{1U, 5U, 3U, 3U, 5U, 7U},
    vector<unsigned char>{9U, 7U, 8U, 9U, 5U, 7U, 10U, 1U, 2U},
    vector<unsigned char>{10U, 1U, 2U, 9U, 5U, 0U, 5U, 3U, 0U, 5U, 7U, 3U},
    vector<unsigned char>{8U, 0U, 2U, 8U, 2U, 5U, 8U, 5U, 7U, 10U, 5U, 2U},
    vector<unsigned char>{2U, 10U, 5U, 2U, 5U, 3U, 3U, 5U, 7U},
    vector<unsigned char>{7U, 9U, 5U, 7U, 8U, 9U, 3U, 11U, 2U},
    vector<unsigned char>{9U, 5U, 7U, 9U, 7U, 2U, 9U, 2U, 0U, 2U, 7U, 11U},
    vector<unsigned char>{2U, 3U, 11U, 0U, 1U, 8U, 1U, 7U, 8U, 1U, 5U, 7U},
    vector<unsigned char>{11U, 2U, 1U, 11U, 1U, 7U, 7U, 1U, 5U},
    vector<unsigned char>{9U, 5U, 8U, 8U, 5U, 7U, 10U, 1U, 3U, 10U, 3U, 11U},
    vector<unsigned char>{5U, 7U, 0U, 5U, 0U, 9U, 7U, 11U, 0U, 1U, 0U, 10U, 11U, 10U, 0U},
    vector<unsigned char>{11U, 10U, 0U, 11U, 0U, 3U, 10U, 5U, 0U, 8U, 0U, 7U, 5U, 7U, 0U},
    vector<unsigned char>{11U, 10U, 5U, 7U, 11U, 5U},
    vector<unsigned char>{10U, 6U, 5U},
    vector<unsigned char>{0U, 8U, 3U, 5U, 10U, 6U},
    vector<unsigned char>{9U, 0U, 1U, 5U, 10U, 6U},
    vector<unsigned char>{1U, 8U, 3U, 1U, 9U, 8U, 5U, 10U, 6U},
    vector<unsigned char>{1U, 6U, 5U, 2U, 6U, 1U},
    vector<unsigned char>{1U, 6U, 5U, 1U, 2U, 6U, 3U, 0U, 8U},
    vector<unsigned char>{9U, 6U, 5U, 9U, 0U, 6U, 0U, 2U, 6U},
    vector<unsigned char>{5U, 9U, 8U, 5U, 8U, 2U, 5U, 2U, 6U, 3U, 2U, 8U},
    vector<unsigned char>{2U, 3U, 11U, 10U, 6U, 5U},
    vector<unsigned char>{11U, 0U, 8U, 11U, 2U, 0U, 10U, 6U, 5U},
    vector<unsigned char>{0U, 1U, 9U, 2U, 3U, 11U, 5U, 10U, 6U},
    vector<unsigned char>{5U, 10U, 6U, 1U, 9U, 2U, 9U, 11U, 2U, 9U, 8U, 11U},
    vector<unsigned char>{6U, 3U, 11U, 6U, 5U, 3U, 5U, 1U, 3U},
    vector<unsigned char>{0U, 8U, 11U, 0U, 11U, 5U, 0U, 5U, 1U, 5U, 11U, 6U},
    vector<unsigned char>{3U, 11U, 6U, 0U, 3U, 6U, 0U, 6U, 5U, 0U, 5U, 9U},
    vector<unsigned char>{6U, 5U, 9U, 6U, 9U, 11U, 11U, 9U, 8U},
    vector<unsigned char>{5U, 10U, 6U, 4U, 7U, 8U},
    vector<unsigned char>{4U, 3U, 0U, 4U, 7U, 3U, 6U, 5U, 10U},
    vector<unsigned char>{1U, 9U, 0U, 5U, 10U, 6U, 8U, 4U, 7U},
    vector<unsigned char>{10U, 6U, 5U, 1U, 9U, 7U, 1U, 7U, 3U, 7U, 9U, 4U},
    vector<unsigned char>{6U, 1U, 2U, 6U, 5U, 1U, 4U, 7U, 8U},
    vector<unsigned char>{1U, 2U, 5U, 5U, 2U, 6U, 3U, 0U, 4U, 3U, 4U, 7U},
    vector<unsigned char>{8U, 4U, 7U, 9U, 0U, 5U, 0U, 6U, 5U, 0U, 2U, 6U},
    vector<unsigned char>{7U, 3U, 9U, 7U, 9U, 4U, 3U, 2U, 9U, 5U, 9U, 6U, 2U, 6U, 9U},
    vector<unsigned char>{3U, 11U, 2U, 7U, 8U, 4U, 10U, 6U, 5U},
    vector<unsigned char>{5U, 10U, 6U, 4U, 7U, 2U, 4U, 2U, 0U, 2U, 7U, 11U},
    vector<unsigned char>{0U, 1U, 9U, 4U, 7U, 8U, 2U, 3U, 11U, 5U, 10U, 6U},
    vector<unsigned char>{9U, 2U, 1U, 9U, 11U, 2U, 9U, 4U, 11U, 7U, 11U, 4U, 5U, 10U, 6U},
    vector<unsigned char>{8U, 4U, 7U, 3U, 11U, 5U, 3U, 5U, 1U, 5U, 11U, 6U},
    vector<unsigned char>{5U, 1U, 11U, 5U, 11U, 6U, 1U, 0U, 11U, 7U, 11U, 4U, 0U, 4U, 11U},
    vector<unsigned char>{0U, 5U, 9U, 0U, 6U, 5U, 0U, 3U, 6U, 11U, 6U, 3U, 8U, 4U, 7U},
    vector<unsigned char>{6U, 5U, 9U, 6U, 9U, 11U, 4U, 7U, 9U, 7U, 11U, 9U},
    vector<unsigned char>{10U, 4U, 9U, 6U, 4U, 10U},
    vector<unsigned char>{4U, 10U, 6U, 4U, 9U, 10U, 0U, 8U, 3U},
    vector<unsigned char>{10U, 0U, 1U, 10U, 6U, 0U, 6U, 4U, 0U},
    vector<unsigned char>{8U, 3U, 1U, 8U, 1U, 6U, 8U, 6U, 4U, 6U, 1U, 10U},
    vector<unsigned char>{1U, 4U, 9U, 1U, 2U, 4U, 2U, 6U, 4U},
    vector<unsigned char>{3U, 0U, 8U, 1U, 2U, 9U, 2U, 4U, 9U, 2U, 6U, 4U},
    vector<unsigned char>{0U, 2U, 4U, 4U, 2U, 6U},
    vector<unsigned char>{8U, 3U, 2U, 8U, 2U, 4U, 4U, 2U, 6U},
    vector<unsigned char>{10U, 4U, 9U, 10U, 6U, 4U, 11U, 2U, 3U},
    vector<unsigned char>{0U, 8U, 2U, 2U, 8U, 11U, 4U, 9U, 10U, 4U, 10U, 6U},
    vector<unsigned char>{3U, 11U, 2U, 0U, 1U, 6U, 0U, 6U, 4U, 6U, 1U, 10U},
    vector<unsigned char>{6U, 4U, 1U, 6U, 1U, 10U, 4U, 8U, 1U, 2U, 1U, 11U, 8U, 11U, 1U},
    vector<unsigned char>{9U, 6U, 4U, 9U, 3U, 6U, 9U, 1U, 3U, 11U, 6U, 3U},
    vector<unsigned char>{8U, 11U, 1U, 8U, 1U, 0U, 11U, 6U, 1U, 9U, 1U, 4U, 6U, 4U, 1U},
    vector<unsigned char>{3U, 11U, 6U, 3U, 6U, 0U, 0U, 6U, 4U},
    vector<unsigned char>{6U, 4U, 8U, 11U, 6U, 8U},
    vector<unsigned char>{7U, 10U, 6U, 7U, 8U, 10U, 8U, 9U, 10U},
    vector<unsigned char>{0U, 7U, 3U, 0U, 10U, 7U, 0U, 9U, 10U, 6U, 7U, 10U},
    vector<unsigned char>{10U, 6U, 7U, 1U, 10U, 7U, 1U, 7U, 8U, 1U, 8U, 0U},
    vector<unsigned char>{10U, 6U, 7U, 10U, 7U, 1U, 1U, 7U, 3U},
    vector<unsigned char>{1U, 2U, 6U, 1U, 6U, 8U, 1U, 8U, 9U, 8U, 6U, 7U},
    vector<unsigned char>{2U, 6U, 9U, 2U, 9U, 1U, 6U, 7U, 9U, 0U, 9U, 3U, 7U, 3U, 9U},
    vector<unsigned char>{7U, 8U, 0U, 7U, 0U, 6U, 6U, 0U, 2U},
    vector<unsigned char>{7U, 3U, 2U, 6U, 7U, 2U},
    vector<unsigned char>{2U, 3U, 11U, 10U, 6U, 8U, 10U, 8U, 9U, 8U, 6U, 7U},
    vector<unsigned char>{2U, 0U, 7U, 2U, 7U, 11U, 0U, 9U, 7U, 6U, 7U, 10U, 9U, 10U, 7U},
    vector<unsigned char>{1U, 8U, 0U, 1U, 7U, 8U, 1U, 10U, 7U, 6U, 7U, 10U, 2U, 3U, 11U},
    vector<unsigned char>{11U, 2U, 1U, 11U, 1U, 7U, 10U, 6U, 1U, 6U, 7U, 1U},
    vector<unsigned char>{8U, 9U, 6U, 8U, 6U, 7U, 9U, 1U, 6U, 11U, 6U, 3U, 1U, 3U, 6U},
    vector<unsigned char>{0U, 9U, 1U, 11U, 6U, 7U},
    vector<unsigned char>{7U, 8U, 0U, 7U, 0U, 6U, 3U, 11U, 0U, 11U, 6U, 0U},
    vector<unsigned char>{7U, 11U, 6U},
    vector<unsigned char>{7U, 6U, 11U},
    vector<unsigned char>{3U, 0U, 8U, 11U, 7U, 6U},
    vector<unsigned char>{0U, 1U, 9U, 11U, 7U, 6U},
    vector<unsigned char>{8U, 1U, 9U, 8U, 3U, 1U, 11U, 7U, 6U},
    vector<unsigned char>{10U, 1U, 2U, 6U, 11U, 7U},
    vector<unsigned char>{1U, 2U, 10U, 3U, 0U, 8U, 6U, 11U, 7U},
    vector<unsigned char>{2U, 9U, 0U, 2U, 10U, 9U, 6U, 11U, 7U},
    vector<unsigned char>{6U, 11U, 7U, 2U, 10U, 3U, 10U, 8U, 3U, 10U, 9U, 8U},
    vector<unsigned char>{7U, 2U, 3U, 6U, 2U, 7U},
    vector<unsigned char>{7U, 0U, 8U, 7U, 6U, 0U, 6U, 2U, 0U},
    vector<unsigned char>{2U, 7U, 6U, 2U, 3U, 7U, 0U, 1U, 9U},
    vector<unsigned char>{1U, 6U, 2U, 1U, 8U, 6U, 1U, 9U, 8U, 8U, 7U, 6U},
    vector<unsigned char>{10U, 7U, 6U, 10U, 1U, 7U, 1U, 3U, 7U},
    vector<unsigned char>{10U, 7U, 6U, 1U, 7U, 10U, 1U, 8U, 7U, 1U, 0U, 8U},
    vector<unsigned char>{0U, 3U, 7U, 0U, 7U, 10U, 0U, 10U, 9U, 6U, 10U, 7U},
    vector<unsigned char>{7U, 6U, 10U, 7U, 10U, 8U, 8U, 10U, 9U},
    vector<unsigned char>{6U, 8U, 4U, 11U, 8U, 6U},
    vector<unsigned char>{3U, 6U, 11U, 3U, 0U, 6U, 0U, 4U, 6U},
    vector<unsigned char>{8U, 6U, 11U, 8U, 4U, 6U, 9U, 0U, 1U},
    vector<unsigned char>{9U, 4U, 6U, 9U, 6U, 3U, 9U, 3U, 1U, 11U, 3U, 6U},
    vector<unsigned char>{6U, 8U, 4U, 6U, 11U, 8U, 2U, 10U, 1U},
    vector<unsigned char>{1U, 2U, 10U, 3U, 0U, 11U, 0U, 6U, 11U, 0U, 4U, 6U},
    vector<unsigned char>{4U, 11U, 8U, 4U, 6U, 11U, 0U, 2U, 9U, 2U, 10U, 9U},
    vector<unsigned char>{10U, 9U, 3U, 10U, 3U, 2U, 9U, 4U, 3U, 11U, 3U, 6U, 4U, 6U, 3U},
    vector<unsigned char>{8U, 2U, 3U, 8U, 4U, 2U, 4U, 6U, 2U},
    vector<unsigned char>{0U, 4U, 2U, 4U, 6U, 2U},
    vector<unsigned char>{1U, 9U, 0U, 2U, 3U, 4U, 2U, 4U, 6U, 4U, 3U, 8U},
    vector<unsigned char>{1U, 9U, 4U, 1U, 4U, 2U, 2U, 4U, 6U},
    vector<unsigned char>{8U, 1U, 3U, 8U, 6U, 1U, 8U, 4U, 6U, 6U, 10U, 1U},
    vector<unsigned char>{10U, 1U, 0U, 10U, 0U, 6U, 6U, 0U, 4U},
    vector<unsigned char>{4U, 6U, 3U, 4U, 3U, 8U, 6U, 10U, 3U, 0U, 3U, 9U, 10U, 9U, 3U},
    vector<unsigned char>{10U, 9U, 4U, 6U, 10U, 4U},
    vector<unsigned char>{4U, 9U, 5U, 7U, 6U, 11U},
    vector<unsigned char>{0U, 8U, 3U, 4U, 9U, 5U, 11U, 7U, 6U},
    vector<unsigned char>{5U, 0U, 1U, 5U, 4U, 0U, 7U, 6U, 11U},
    vector<unsigned char>{11U, 7U, 6U, 8U, 3U, 4U, 3U, 5U, 4U, 3U, 1U, 5U},
    vector<unsigned char>{9U, 5U, 4U, 10U, 1U, 2U, 7U, 6U, 11U},
    vector<unsigned char>{6U, 11U, 7U, 1U, 2U, 10U, 0U, 8U, 3U, 4U, 9U, 5U},
    vector<unsigned char>{7U, 6U, 11U, 5U, 4U, 10U, 4U, 2U, 10U, 4U, 0U, 2U},
    vector<unsigned char>{3U, 4U, 8U, 3U, 5U, 4U, 3U, 2U, 5U, 10U, 5U, 2U, 11U, 7U, 6U},
    vector<unsigned char>{7U, 2U, 3U, 7U, 6U, 2U, 5U, 4U, 9U},
    vector<unsigned char>{9U, 5U, 4U, 0U, 8U, 6U, 0U, 6U, 2U, 6U, 8U, 7U},
    vector<unsigned char>{3U, 6U, 2U, 3U, 7U, 6U, 1U, 5U, 0U, 5U, 4U, 0U},
    vector<unsigned char>{6U, 2U, 8U, 6U, 8U, 7U, 2U, 1U, 8U, 4U, 8U, 5U, 1U, 5U, 8U},
    vector<unsigned char>{9U, 5U, 4U, 10U, 1U, 6U, 1U, 7U, 6U, 1U, 3U, 7U},
    vector<unsigned char>{1U, 6U, 10U, 1U, 7U, 6U, 1U, 0U, 7U, 8U, 7U, 0U, 9U, 5U, 4U},
    vector<unsigned char>{4U, 0U, 10U, 4U, 10U, 5U, 0U, 3U, 10U, 6U, 10U, 7U, 3U, 7U, 10U},
    vector<unsigned char>{7U, 6U, 10U, 7U, 10U, 8U, 5U, 4U, 10U, 4U, 8U, 10U},
    vector<unsigned char>{6U, 9U, 5U, 6U, 11U, 9U, 11U, 8U, 9U},
    vector<unsigned char>{3U, 6U, 11U, 0U, 6U, 3U, 0U, 5U, 6U, 0U, 9U, 5U},
    vector<unsigned char>{0U, 11U, 8U, 0U, 5U, 11U, 0U, 1U, 5U, 5U, 6U, 11U},
    vector<unsigned char>{6U, 11U, 3U, 6U, 3U, 5U, 5U, 3U, 1U},
    vector<unsigned char>{1U, 2U, 10U, 9U, 5U, 11U, 9U, 11U, 8U, 11U, 5U, 6U},
    vector<unsigned char>{0U, 11U, 3U, 0U, 6U, 11U, 0U, 9U, 6U, 5U, 6U, 9U, 1U, 2U, 10U},
    vector<unsigned char>{11U, 8U, 5U, 11U, 5U, 6U, 8U, 0U, 5U, 10U, 5U, 2U, 0U, 2U, 5U},
    vector<unsigned char>{6U, 11U, 3U, 6U, 3U, 5U, 2U, 10U, 3U, 10U, 5U, 3U},
    vector<unsigned char>{5U, 8U, 9U, 5U, 2U, 8U, 5U, 6U, 2U, 3U, 8U, 2U},
    vector<unsigned char>{9U, 5U, 6U, 9U, 6U, 0U, 0U, 6U, 2U},
    vector<unsigned char>{1U, 5U, 8U, 1U, 8U, 0U, 5U, 6U, 8U, 3U, 8U, 2U, 6U, 2U, 8U},
    vector<unsigned char>{1U, 5U, 6U, 2U, 1U, 6U},
    vector<unsigned char>{1U, 3U, 6U, 1U, 6U, 10U, 3U, 8U, 6U, 5U, 6U, 9U, 8U, 9U, 6U},
    vector<unsigned char>{10U, 1U, 0U, 10U, 0U, 6U, 9U, 5U, 0U, 5U, 6U, 0U},
    vector<unsigned char>{0U, 3U, 8U, 5U, 6U, 10U},
    vector<unsigned char>{10U, 5U, 6U},
    vector<unsigned char>{11U, 5U, 10U, 7U, 5U, 11U},
    vector<unsigned char>{11U, 5U, 10U, 11U, 7U, 5U, 8U, 3U, 0U},
    vector<unsigned char>{5U, 11U, 7U, 5U, 10U, 11U, 1U, 9U, 0U},
    vector<unsigned char>{10U, 7U, 5U, 10U, 11U, 7U, 9U, 8U, 1U, 8U, 3U, 1U},
    vector<unsigned char>{11U, 1U, 2U, 11U, 7U, 1U, 7U, 5U, 1U},
    vector<unsigned char>{0U, 8U, 3U, 1U, 2U, 7U, 1U, 7U, 5U, 7U, 2U, 11U},
    vector<unsigned char>{9U, 7U, 5U, 9U, 2U, 7U, 9U, 0U, 2U, 2U, 11U, 7U},
    vector<unsigned char>{7U, 5U, 2U, 7U, 2U, 11U, 5U, 9U, 2U, 3U, 2U, 8U, 9U, 8U, 2U},
    vector<unsigned char>{2U, 5U, 10U, 2U, 3U, 5U, 3U, 7U, 5U},
    vector<unsigned char>{8U, 2U, 0U, 8U, 5U, 2U, 8U, 7U, 5U, 10U, 2U, 5U},
    vector<unsigned char>{9U, 0U, 1U, 5U, 10U, 3U, 5U, 3U, 7U, 3U, 10U, 2U},
    vector<unsigned char>{9U, 8U, 2U, 9U, 2U, 1U, 8U, 7U, 2U, 10U, 2U, 5U, 7U, 5U, 2U},
    vector<unsigned char>{1U, 3U, 5U, 3U, 7U, 5U},
    vector<unsigned char>{0U, 8U, 7U, 0U, 7U, 1U, 1U, 7U, 5U},
    vector<unsigned char>{9U, 0U, 3U, 9U, 3U, 5U, 5U, 3U, 7U},
    vector<unsigned char>{9U, 8U, 7U, 5U, 9U, 7U},
    vector<unsigned char>{5U, 8U, 4U, 5U, 10U, 8U, 10U, 11U, 8U},
    vector<unsigned char>{5U, 0U, 4U, 5U, 11U, 0U, 5U, 10U, 11U, 11U, 3U, 0U},
    vector<unsigned char>{0U, 1U, 9U, 8U, 4U, 10U, 8U, 10U, 11U, 10U, 4U, 5U},
    vector<unsigned char>{10U, 11U, 4U, 10U, 4U, 5U, 11U, 3U, 4U, 9U, 4U, 1U, 3U, 1U, 4U},
    vector<unsigned char>{2U, 5U, 1U, 2U, 8U, 5U, 2U, 11U, 8U, 4U, 5U, 8U},
    vector<unsigned char>{0U, 4U, 11U, 0U, 11U, 3U, 4U, 5U, 11U, 2U, 11U, 1U, 5U, 1U, 11U},
    vector<unsigned char>{0U, 2U, 5U, 0U, 5U, 9U, 2U, 11U, 5U, 4U, 5U, 8U, 11U, 8U, 5U},
    vector<unsigned char>{9U, 4U, 5U, 2U, 11U, 3U},
    vector<unsigned char>{2U, 5U, 10U, 3U, 5U, 2U, 3U, 4U, 5U, 3U, 8U, 4U},
    vector<unsigned char>{5U, 10U, 2U, 5U, 2U, 4U, 4U, 2U, 0U},
    vector<unsigned char>{3U, 10U, 2U, 3U, 5U, 10U, 3U, 8U, 5U, 4U, 5U, 8U, 0U, 1U, 9U},
    vector<unsigned char>{5U, 10U, 2U, 5U, 2U, 4U, 1U, 9U, 2U, 9U, 4U, 2U},
    vector<unsigned char>{8U, 4U, 5U, 8U, 5U, 3U, 3U, 5U, 1U},
    vector<unsigned char>{0U, 4U, 5U, 1U, 0U, 5U},
    vector<unsigned char>{8U, 4U, 5U, 8U, 5U, 3U, 9U, 0U, 5U, 0U, 3U, 5U},
    vector<unsigned char>{9U, 4U, 5U},
    vector<unsigned char>{4U, 11U, 7U, 4U, 9U, 11U, 9U, 10U, 11U},
    vector<unsigned char>{0U, 8U, 3U, 4U, 9U, 7U, 9U, 11U, 7U, 9U, 10U, 11U},
    vector<unsigned char>{1U, 10U, 11U, 1U, 11U, 4U, 1U, 4U, 0U, 7U, 4U, 11U},
    vector<unsigned char>{3U, 1U, 4U, 3U, 4U, 8U, 1U, 10U, 4U, 7U, 4U, 11U, 10U, 11U, 4U},
    vector<unsigned char>{4U, 11U, 7U, 9U, 11U, 4U, 9U, 2U, 11U, 9U, 1U, 2U},
    vector<unsigned char>{9U, 7U, 4U, 9U, 11U, 7U, 9U, 1U, 11U, 2U, 11U, 1U, 0U, 8U, 3U},
    vector<unsigned char>{11U, 7U, 4U, 11U, 4U, 2U, 2U, 4U, 0U},
    vector<unsigned char>{11U, 7U, 4U, 11U, 4U, 2U, 8U, 3U, 4U, 3U, 2U, 4U},
    vector<unsigned char>{2U, 9U, 10U, 2U, 7U, 9U, 2U, 3U, 7U, 7U, 4U, 9U},
    vector<unsigned char>{9U, 10U, 7U, 9U, 7U, 4U, 10U, 2U, 7U, 8U, 7U, 0U, 2U, 0U, 7U},
    vector<unsigned char>{3U, 7U, 10U, 3U, 10U, 2U, 7U, 4U, 10U, 1U, 10U, 0U, 4U, 0U, 10U},
    vector<unsigned char>{1U, 10U, 2U, 8U, 7U, 4U},
    vector<unsigned char>{4U, 9U, 1U, 4U, 1U, 7U, 7U, 1U, 3U},
    vector<unsigned char>{4U, 9U, 1U, 4U, 1U, 7U, 0U, 8U, 1U, 8U, 7U, 1U},
    vector<unsigned char>{4U, 0U, 3U, 7U, 4U, 3U},
    vector<unsigned char>{4U, 8U, 7U},
    vector<unsigned char>{9U, 10U, 8U, 10U, 11U, 8U},
    vector<unsigned char>{3U, 0U, 9U, 3U, 9U, 11U, 11U, 9U, 10U},
    vector<unsigned char>{0U, 1U, 10U, 0U, 10U, 8U, 8U, 10U, 11U},
    vector<unsigned char>{3U, 1U, 10U, 11U, 3U, 10U},
    vector<unsigned char>{1U, 2U, 11U, 1U, 11U, 9U, 9U, 11U, 8U},
    vector<unsigned char>{3U, 0U, 9U, 3U, 9U, 11U, 1U, 2U, 9U, 2U, 11U, 9U},
    vector<unsigned char>{0U, 2U, 11U, 8U, 0U, 11U},
    vector<unsigned char>{3U, 2U, 11U},
    vector<unsigned char>{2U, 3U, 8U, 2U, 8U, 10U, 10U, 8U, 9U},
    vector<unsigned char>{9U, 10U, 2U, 0U, 9U, 2U},
    vector<unsigned char>{2U, 3U, 8U, 2U, 8U, 10U, 0U, 1U, 8U, 1U, 10U, 8U},
    vector<unsigned char>{1U, 10U, 2U},
    vector<unsigned char>{1U, 3U, 8U, 9U, 1U, 8U},
    vector<unsigned char>{0U, 9U, 1U},
    vector<unsigned char>{0U, 3U, 8U},
    vector<unsigned char>{}
  };
  static const unsigned short ushEdgePattern[1U << GRID_VERTEX_COUNT_3D] = {0x0U, 0x109U, 0x203U, 0x30AU, 0x406U, 0x50FU, 0x605U, 0x70CU, 0x80CU, 0x905U, 0xA0FU, 0xB06U, 0xC0AU, 0xD03U, 0xE09U, 0xF00U, 0x190U, 0x99U, 0x393U, 0x29AU, 0x596U, 0x49FU, 0x795U, 0x69CU, 0x99CU, 0x895U, 0xB9FU, 0xA96U, 0xD9AU, 0xC93U, 0xF99U, 0xE90U, 0x230U, 0x339U, 0x33U, 0x13AU, 0x636U, 0x73FU, 0x435U, 0x53CU, 0xA3CU, 0xB35U, 0x83FU, 0x936U, 0xE3AU, 0xF33U, 0xC39U, 0xD30U, 0x3A0U, 0x2A9U, 0x1A3U, 0xAAU, 0x7A6U, 0x6AFU, 0x5A5U, 0x4ACU, 0xBACU, 0xAA5U, 0x9AFU, 0x8A6U, 0xFAAU, 0xEA3U, 0xDA9U, 0xCA0U, 0x460U, 0x569U, 0x663U, 0x76AU, 0x66U, 0x16FU, 0x265U, 0x36CU, 0xC6CU, 0xD65U, 0xE6FU, 0xF66U, 0x86AU, 0x963U, 0xA69U, 0xB60U, 0x5F0U, 0x4F9U, 0x7F3U, 0x6FAU, 0x1F6U, 0xFFU, 0x3F5U, 0x2FCU, 0xDFCU, 0xCF5U, 0xFFFU, 0xEF6U, 0x9FAU, 0x8F3U, 0xBF9U, 0xAF0U, 0x650U, 0x759U, 0x453U, 0x55AU, 0x256U, 0x35FU, 0x55U, 0x15CU, 0xE5CU, 0xF55U, 0xC5FU, 0xD56U, 0xA5AU, 0xB53U, 0x859U, 0x950U, 0x7C0U, 0x6C9U, 0x5C3U, 0x4CAU, 0x3C6U, 0x2CFU, 0x1C5U, 0xCCU, 0xFCCU, 0xEC5U, 0xDCFU, 0xCC6U, 0xBCAU, 0xAC3U, 0x9C9U, 0x8C0U, 0x8C0U, 0x9C9U, 0xAC3U, 0xBCAU, 0xCC6U, 0xDCFU, 0xEC5U, 0xFCCU, 0xCCU, 0x1C5U, 0x2CFU, 0x3C6U, 0x4CAU, 0x5C3U, 0x6C9U, 0x7C0U, 0x950U, 0x859U, 0xB53U, 0xA5AU, 0xD56U, 0xC5FU, 0xF55U, 0xE5CU, 0x15CU, 0x55U, 0x35FU, 0x256U, 0x55AU, 0x453U, 0x759U, 0x650U, 0xAF0U, 0xBF9U, 0x8F3U, 0x9FAU, 0xEF6U, 0xFFFU, 0xCF5U, 0xDFCU, 0x2FCU, 0x3F5U, 0xFFU, 0x1F6U, 0x6FAU, 0x7F3U, 0x4F9U, 0x5F0U, 0xB60U, 0xA69U, 0x963U, 0x86AU, 0xF66U, 0xE6FU, 0xD65U, 0xC6CU, 0x36CU, 0x265U, 0x16FU, 0x66U, 0x76AU, 0x663U, 0x569U, 0x460U, 0xCA0U, 0xDA9U, 0xEA3U, 0xFAAU, 0x8A6U, 0x9AFU, 0xAA5U, 0xBACU, 0x4ACU, 0x5A5U, 0x6AFU, 0x7A6U, 0xAAU, 0x1A3U, 0x2A9U, 0x3A0U, 0xD30U, 0xC39U, 0xF33U, 0xE3AU, 0x936U, 0x83FU, 0xB35U, 0xA3CU, 0x53CU, 0x435U, 0x73FU, 0x636U, 0x13AU, 0x33U, 0x339U, 0x230U, 0xE90U, 0xF99U, 0xC93U, 0xD9AU, 0xA96U, 0xB9FU, 0x895U, 0x99CU, 0x69CU, 0x795U, 0x49FU, 0x596U, 0x29AU, 0x393U, 0x99U, 0x190U, 0xF00U, 0xE09U, 0xD03U, 0xC0AU, 0xB06U, 0xA0FU, 0x905U, 0x80CU, 0x70CU, 0x605U, 0x50FU, 0x406U, 0x30AU, 0x203U, 0x109U, 0x0U};

  vector<Vertex3D> vecTransformedCorner{rPiecewisePolynomial3DData.MinimumCorner(), rPiecewisePolynomial3DData.MaximumCorner()};
  float*** pfGridValueMatrix = new float**[riDivisionCountX + 1];
  float fInversedTransformationMatrix[ELEMENT_COUNT_3D + 1][ELEMENT_COUNT_3D + 1];
  int iTriangleCount = 0;

  for (int i = 0; i <= riDivisionCountX; i++) {
    pfGridValueMatrix[i] = new float*[riDivisionCountY + 1];

    for (int j = 0; j <= riDivisionCountY; j++)
      pfGridValueMatrix[i][j] = new float[riDivisionCountZ + 1];
  }

  if (rpfTransformationMatrix == nullptr)
    ;
  else {
    for (int i = 0; i <= ELEMENT_COUNT_2D; i++)
      for (int j = 0; j <= ELEMENT_COUNT_2D; j++)
        fInversedTransformationMatrix[i][j] = rpfTransformationMatrix[i][j];

    InverseTransformation<ELEMENT_COUNT_3D + 1>(fInversedTransformationMatrix);

    float Vertex3D::*pfVertex3DMember[ELEMENT_COUNT_3D] = {&Vertex3D::fX, &Vertex3D::fY, &Vertex3D::fZ};
    vector<Vertex3D> vecOriginalCorner{rPiecewisePolynomial3DData.MinimumCorner(), rPiecewisePolynomial3DData.MaximumCorner()};

    for (int i = 0; i < ELEMENT_COUNT_3D; i++) {
      vecTransformedCorner[INDICATOR_MINIMUM].*pfVertex3DMember[i] = vecTransformedCorner[INDICATOR_MAXIMUM].*pfVertex3DMember[i] = 0.0F;

      for (int j = 0; j < ELEMENT_COUNT_3D; j++)
        if (rpfTransformationMatrix[i][j] * vecOriginalCorner[INDICATOR_MINIMUM].*pfVertex3DMember[j] < rpfTransformationMatrix[i][j] * vecOriginalCorner[INDICATOR_MAXIMUM].*pfVertex3DMember[j]) {
          vecTransformedCorner[INDICATOR_MINIMUM].*pfVertex3DMember[i] += rpfTransformationMatrix[i][j] * vecOriginalCorner[INDICATOR_MINIMUM].*pfVertex3DMember[j];
          vecTransformedCorner[INDICATOR_MAXIMUM].*pfVertex3DMember[i] += rpfTransformationMatrix[i][j] * vecOriginalCorner[INDICATOR_MAXIMUM].*pfVertex3DMember[j];
        } else {
          vecTransformedCorner[INDICATOR_MINIMUM].*pfVertex3DMember[i] += rpfTransformationMatrix[i][j] * vecOriginalCorner[INDICATOR_MAXIMUM].*pfVertex3DMember[j];
          vecTransformedCorner[INDICATOR_MAXIMUM].*pfVertex3DMember[i] += rpfTransformationMatrix[i][j] * vecOriginalCorner[INDICATOR_MINIMUM].*pfVertex3DMember[j];
        }

      vecTransformedCorner[INDICATOR_MINIMUM].*pfVertex3DMember[i] += rpfTransformationMatrix[i][ELEMENT_COUNT_3D];
      vecTransformedCorner[INDICATOR_MAXIMUM].*pfVertex3DMember[i] += rpfTransformationMatrix[i][ELEMENT_COUNT_3D];
    }
  }

  float fStep[ELEMENT_COUNT_3D] = {(vecTransformedCorner[INDICATOR_MAXIMUM].fX - vecTransformedCorner[INDICATOR_MINIMUM].fX) / static_cast<float>(riDivisionCountX), (vecTransformedCorner[INDICATOR_MAXIMUM].fY - vecTransformedCorner[INDICATOR_MINIMUM].fY) / static_cast<float>(riDivisionCountY), (vecTransformedCorner[INDICATOR_MAXIMUM].fZ - vecTransformedCorner[INDICATOR_MINIMUM].fZ) / static_cast<float>(riDivisionCountZ)};

  for (int i = 0; i <= riDivisionCountX; i++)
    for (int j = 0; j <= riDivisionCountY; j++)
      for (int k = 0; k <= riDivisionCountZ; k++) {
        Vertex3D Vertex3DSamplingPoint{vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * static_cast<float>(i), vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * static_cast<float>(j), vecTransformedCorner[INDICATOR_MINIMUM].fZ + fStep[INDICATOR_Z] * static_cast<float>(k)};

        pfGridValueMatrix[i][j][k] = rPiecewisePolynomial3DData.At(rpfTransformationMatrix == nullptr ? Vertex3DSamplingPoint : TransformCoordinate(fInversedTransformationMatrix, Vertex3DSamplingPoint));
      }

  unsigned char uchOffset[GRID_VERTEX_COUNT_3D][ELEMENT_COUNT_3D] = {
    {0U, 0U, 0U},
    {0U, 0U, 1U},
    {0U, 1U, 1U},
    {0U, 1U, 0U},
    {1U, 0U, 0U},
    {1U, 0U, 1U},
    {1U, 1U, 1U},
    {1U, 1U, 0U}
  };

  for (int i = 0; i < riDivisionCountX; i++)
    for (int j = 0; j < riDivisionCountY; j++)
      for (int k = 0; k < riDivisionCountZ; k++) {
        unsigned char uchFlag = 0U;

        do {
          unsigned char uchDigit = 1U;

          for (const unsigned char* const& rpuchPointer : uchOffset) {
            if (pfGridValueMatrix[rpuchPointer[INDICATOR_X] + i][rpuchPointer[INDICATOR_Y] + j][rpuchPointer[INDICATOR_Z] + k] >= 0.0F)
              uchFlag |= uchDigit;

            uchDigit <<= 1U;
          }
        } while (false);

        if (!ushEdgePattern[uchFlag])
          continue;

        vector<Vertex3D> vecInterpolationPoint(GRID_EDGE_COUNT_3D, Vertex3D{0.0F, 0.0F, 0.0F});

        do {
          unsigned char uchFactor = 0U;

          if (ushEdgePattern[uchFlag] & 1U << uchFactor++)
            vecInterpolationPoint[uchFactor - 1] = Vertex3D{vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * static_cast<float>(i), vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * static_cast<float>(j), (pfGridValueMatrix[i][j][k + 1] * (vecTransformedCorner[INDICATOR_MINIMUM].fZ + fStep[INDICATOR_Z] * static_cast<float>(k)) - pfGridValueMatrix[i][j][k] * (vecTransformedCorner[INDICATOR_MINIMUM].fZ + fStep[INDICATOR_Z] * (static_cast<float>(k) + 1.0F))) / (pfGridValueMatrix[i][j][k + 1] - pfGridValueMatrix[i][j][k])};

          if (ushEdgePattern[uchFlag] & 1U << uchFactor++)
            vecInterpolationPoint[uchFactor - 1] = Vertex3D{vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * static_cast<float>(i), (pfGridValueMatrix[i][j + 1][k + 1] * (vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * static_cast<float>(j)) - pfGridValueMatrix[i][j][k + 1] * (vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * (static_cast<float>(j) + 1.0F))) / (pfGridValueMatrix[i][j + 1][k + 1] - pfGridValueMatrix[i][j][k + 1]), vecTransformedCorner[INDICATOR_MINIMUM].fZ + fStep[INDICATOR_Z] * (static_cast<float>(k) + 1.0F)};

          if (ushEdgePattern[uchFlag] & 1U << uchFactor++)
            vecInterpolationPoint[uchFactor - 1] = Vertex3D{vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * static_cast<float>(i), vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * (static_cast<float>(j) + 1.0F), (pfGridValueMatrix[i][j + 1][k + 1] * (vecTransformedCorner[INDICATOR_MINIMUM].fZ + fStep[INDICATOR_Z] * static_cast<float>(k)) - pfGridValueMatrix[i][j + 1][k] * (vecTransformedCorner[INDICATOR_MINIMUM].fZ + fStep[INDICATOR_Z] * (static_cast<float>(k) + 1.0F))) / (pfGridValueMatrix[i][j + 1][k + 1] - pfGridValueMatrix[i][j + 1][k])};

          if (ushEdgePattern[uchFlag] & 1U << uchFactor++)
            vecInterpolationPoint[uchFactor - 1] = Vertex3D{vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * static_cast<float>(i), (pfGridValueMatrix[i][j + 1][k] * (vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * static_cast<float>(j)) - pfGridValueMatrix[i][j][k] * (vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * (static_cast<float>(j) + 1.0F))) / (pfGridValueMatrix[i][j + 1][k] - pfGridValueMatrix[i][j][k]), vecTransformedCorner[INDICATOR_MINIMUM].fZ + fStep[INDICATOR_Z] * static_cast<float>(k)};

          if (ushEdgePattern[uchFlag] & 1U << uchFactor++)
            vecInterpolationPoint[uchFactor - 1] = Vertex3D{vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * (static_cast<float>(i) + 1.0F), vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * static_cast<float>(j), (pfGridValueMatrix[i + 1][j][k + 1] * (vecTransformedCorner[INDICATOR_MINIMUM].fZ + fStep[INDICATOR_Z] * static_cast<float>(k)) - pfGridValueMatrix[i + 1][j][k] * (vecTransformedCorner[INDICATOR_MINIMUM].fZ + fStep[INDICATOR_Z] * (static_cast<float>(k) + 1.0F))) / (pfGridValueMatrix[i + 1][j][k + 1] - pfGridValueMatrix[i + 1][j][k])};

          if (ushEdgePattern[uchFlag] & 1U << uchFactor++)
            vecInterpolationPoint[uchFactor - 1] = Vertex3D{vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * (static_cast<float>(i) + 1.0F), (pfGridValueMatrix[i + 1][j + 1][k + 1] * (vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * static_cast<float>(j)) - pfGridValueMatrix[i + 1][j][k + 1] * (vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * (static_cast<float>(j) + 1.0F))) / (pfGridValueMatrix[i + 1][j + 1][k + 1] - pfGridValueMatrix[i + 1][j][k + 1]), vecTransformedCorner[INDICATOR_MINIMUM].fZ + fStep[INDICATOR_Z] * (static_cast<float>(k) + 1.0F)};

          if (ushEdgePattern[uchFlag] & 1U << uchFactor++)
            vecInterpolationPoint[uchFactor - 1] = Vertex3D{vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * (static_cast<float>(i) + 1.0F), vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * (static_cast<float>(j) + 1.0F), (pfGridValueMatrix[i + 1][j + 1][k + 1] * (vecTransformedCorner[INDICATOR_MINIMUM].fZ + fStep[INDICATOR_Z] * static_cast<float>(k)) - pfGridValueMatrix[i + 1][j + 1][k] * (vecTransformedCorner[INDICATOR_MINIMUM].fZ + fStep[INDICATOR_Z] * (static_cast<float>(k) + 1.0F))) / (pfGridValueMatrix[i + 1][j + 1][k + 1] - pfGridValueMatrix[i + 1][j + 1][k])};

          if (ushEdgePattern[uchFlag] & 1U << uchFactor++)
            vecInterpolationPoint[uchFactor - 1] = Vertex3D{vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * (static_cast<float>(i) + 1.0F), (pfGridValueMatrix[i + 1][j + 1][k] * (vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * static_cast<float>(j)) - pfGridValueMatrix[i + 1][j][k] * (vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * (static_cast<float>(j) + 1.0F))) / (pfGridValueMatrix[i + 1][j + 1][k] - pfGridValueMatrix[i + 1][j][k]), vecTransformedCorner[INDICATOR_MINIMUM].fZ + fStep[INDICATOR_Z] * static_cast<float>(k)};

          if (ushEdgePattern[uchFlag] & 1U << uchFactor++)
            vecInterpolationPoint[uchFactor - 1] = Vertex3D{(pfGridValueMatrix[i + 1][j][k] * (vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * static_cast<float>(i)) - pfGridValueMatrix[i][j][k] * (vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * (static_cast<float>(i) + 1.0F))) / (pfGridValueMatrix[i + 1][j][k] - pfGridValueMatrix[i][j][k]), vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * static_cast<float>(j), vecTransformedCorner[INDICATOR_MINIMUM].fZ + fStep[INDICATOR_Z] * static_cast<float>(k)};

          if (ushEdgePattern[uchFlag] & 1U << uchFactor++)
            vecInterpolationPoint[uchFactor - 1] = Vertex3D{(pfGridValueMatrix[i + 1][j][k + 1] * (vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * static_cast<float>(i)) - pfGridValueMatrix[i][j][k + 1] * (vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * (static_cast<float>(i) + 1.0F))) / (pfGridValueMatrix[i + 1][j][k + 1] - pfGridValueMatrix[i][j][k + 1]), vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * static_cast<float>(j), vecTransformedCorner[INDICATOR_MINIMUM].fZ + fStep[INDICATOR_Z] * (static_cast<float>(k) + 1.0F)};

          if (ushEdgePattern[uchFlag] & 1U << uchFactor++)
            vecInterpolationPoint[uchFactor - 1] = Vertex3D{(pfGridValueMatrix[i + 1][j + 1][k + 1] * (vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * static_cast<float>(i)) - pfGridValueMatrix[i][j + 1][k + 1] * (vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * (static_cast<float>(i) + 1.0F))) / (pfGridValueMatrix[i + 1][j + 1][k + 1] - pfGridValueMatrix[i][j + 1][k + 1]), vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * (static_cast<float>(j) + 1.0F), vecTransformedCorner[INDICATOR_MINIMUM].fZ + fStep[INDICATOR_Z] * (static_cast<float>(k) + 1.0F)};

          if (ushEdgePattern[uchFlag] & 1U << uchFactor++)
            vecInterpolationPoint[uchFactor - 1] = Vertex3D{(pfGridValueMatrix[i + 1][j + 1][k] * (vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * static_cast<float>(i)) - pfGridValueMatrix[i][j + 1][k] * (vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * (static_cast<float>(i) + 1.0F))) / (pfGridValueMatrix[i + 1][j + 1][k] - pfGridValueMatrix[i][j + 1][k]), vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * (static_cast<float>(j) + 1.0F), vecTransformedCorner[INDICATOR_MINIMUM].fZ + fStep[INDICATOR_Z] * static_cast<float>(k)};
        } while (false);

        for (int l = 0; l * TRIANGLE_ELEMENT_COUNT < vecTrianglePattern[uchFlag].size(); l++) {
          for (int m = 0; m < TRIANGLE_ELEMENT_COUNT; m++)
            rofsFile << VERTEX_CHARACTER << ' ' << vecInterpolationPoint[vecTrianglePattern[uchFlag][l * TRIANGLE_ELEMENT_COUNT + m]].fX << ' ' << vecInterpolationPoint[vecTrianglePattern[uchFlag][l * TRIANGLE_ELEMENT_COUNT + m]].fY << ' ' << vecInterpolationPoint[vecTrianglePattern[uchFlag][l * TRIANGLE_ELEMENT_COUNT + m]].fZ << endl;

          iTriangleCount++;
        }
      }

  for (int i = 0; i <= riDivisionCountX; i++) {
    for (int j = 0; j <= riDivisionCountY; j++)
      delete[] pfGridValueMatrix[i][j];

    delete[] pfGridValueMatrix[i];
  }

  delete[] pfGridValueMatrix;

  int iVertexCount = 0;

  for (int i = 0; i < iTriangleCount; i++) {
    rofsFile << FACE_ELEMENT_CHARACTER << ' ';

    for (int j = 0; j < TRIANGLE_ELEMENT_COUNT; j++) {
      rofsFile << iVertexCount + j + 1 << DELIMITER_CHARACTER << DELIMITER_CHARACTER;

      if (TRIANGLE_ELEMENT_COUNT - j - 1)
        rofsFile << ' ';
      else
        rofsFile << endl;
    }

    iVertexCount += TRIANGLE_ELEMENT_COUNT;
  }
}

void MarchingCube(ofstream& rofsFile, const vector<PiecewisePolynomial3D>& rvecImplicitSurface, const float (*const& rpfTransformationMatrix)[ELEMENT_COUNT_3D + 1][ELEMENT_COUNT_3D + 1], const int& riDivisionCountX, const int& riDivisionCountY, const int& riDivisionCountZ)
{
  static const vector<vector<unsigned char>> vecTrianglePattern{
    vector<unsigned char>{},
    vector<unsigned char>{0U, 8U, 3U},
    vector<unsigned char>{0U, 1U, 9U},
    vector<unsigned char>{1U, 8U, 3U, 9U, 8U, 1U},
    vector<unsigned char>{1U, 2U, 10U},
    vector<unsigned char>{0U, 8U, 3U, 1U, 2U, 10U},
    vector<unsigned char>{9U, 2U, 10U, 0U, 2U, 9U},
    vector<unsigned char>{2U, 8U, 3U, 2U, 10U, 8U, 10U, 9U, 8U},
    vector<unsigned char>{3U, 11U, 2U},
    vector<unsigned char>{0U, 11U, 2U, 8U, 11U, 0U},
    vector<unsigned char>{1U, 9U, 0U, 2U, 3U, 11U},
    vector<unsigned char>{1U, 11U, 2U, 1U, 9U, 11U, 9U, 8U, 11U},
    vector<unsigned char>{3U, 10U, 1U, 11U, 10U, 3U},
    vector<unsigned char>{0U, 10U, 1U, 0U, 8U, 10U, 8U, 11U, 10U},
    vector<unsigned char>{3U, 9U, 0U, 3U, 11U, 9U, 11U, 10U, 9U},
    vector<unsigned char>{9U, 8U, 10U, 10U, 8U, 11U},
    vector<unsigned char>{4U, 7U, 8U},
    vector<unsigned char>{4U, 3U, 0U, 7U, 3U, 4U},
    vector<unsigned char>{0U, 1U, 9U, 8U, 4U, 7U},
    vector<unsigned char>{4U, 1U, 9U, 4U, 7U, 1U, 7U, 3U, 1U},
    vector<unsigned char>{1U, 2U, 10U, 8U, 4U, 7U},
    vector<unsigned char>{3U, 4U, 7U, 3U, 0U, 4U, 1U, 2U, 10U},
    vector<unsigned char>{9U, 2U, 10U, 9U, 0U, 2U, 8U, 4U, 7U},
    vector<unsigned char>{2U, 10U, 9U, 2U, 9U, 7U, 2U, 7U, 3U, 7U, 9U, 4U},
    vector<unsigned char>{8U, 4U, 7U, 3U, 11U, 2U},
    vector<unsigned char>{11U, 4U, 7U, 11U, 2U, 4U, 2U, 0U, 4U},
    vector<unsigned char>{9U, 0U, 1U, 8U, 4U, 7U, 2U, 3U, 11U},
    vector<unsigned char>{4U, 7U, 11U, 9U, 4U, 11U, 9U, 11U, 2U, 9U, 2U, 1U},
    vector<unsigned char>{3U, 10U, 1U, 3U, 11U, 10U, 7U, 8U, 4U},
    vector<unsigned char>{1U, 11U, 10U, 1U, 4U, 11U, 1U, 0U, 4U, 7U, 11U, 4U},
    vector<unsigned char>{4U, 7U, 8U, 9U, 0U, 11U, 9U, 11U, 10U, 11U, 0U, 3U},
    vector<unsigned char>{4U, 7U, 11U, 4U, 11U, 9U, 9U, 11U, 10U},
    vector<unsigned char>{9U, 5U, 4U},
    vector<unsigned char>{9U, 5U, 4U, 0U, 8U, 3U},
    vector<unsigned char>{0U, 5U, 4U, 1U, 5U, 0U},
    vector<unsigned char>{8U, 5U, 4U, 8U, 3U, 5U, 3U, 1U, 5U},
    vector<unsigned char>{1U, 2U, 10U, 9U, 5U, 4U},
    vector<unsigned char>{3U, 0U, 8U, 1U, 2U, 10U, 4U, 9U, 5U},
    vector<unsigned char>{5U, 2U, 10U, 5U, 4U, 2U, 4U, 0U, 2U},
    vector<unsigned char>{2U, 10U, 5U, 3U, 2U, 5U, 3U, 5U, 4U, 3U, 4U, 8U},
    vector<unsigned char>{9U, 5U, 4U, 2U, 3U, 11U},
    vector<unsigned char>{0U, 11U, 2U, 0U, 8U, 11U, 4U, 9U, 5U},
    vector<unsigned char>{0U, 5U, 4U, 0U, 1U, 5U, 2U, 3U, 11U},
    vector<unsigned char>{2U, 1U, 5U, 2U, 5U, 8U, 2U, 8U, 11U, 4U, 8U, 5U},
    vector<unsigned char>{10U, 3U, 11U, 10U, 1U, 3U, 9U, 5U, 4U},
    vector<unsigned char>{4U, 9U, 5U, 0U, 8U, 1U, 8U, 10U, 1U, 8U, 11U, 10U},
    vector<unsigned char>{5U, 4U, 0U, 5U, 0U, 11U, 5U, 11U, 10U, 11U, 0U, 3U},
    vector<unsigned char>{5U, 4U, 8U, 5U, 8U, 10U, 10U, 8U, 11U},
    vector<unsigned char>{9U, 7U, 8U, 5U, 7U, 9U},
    vector<unsigned char>{9U, 3U, 0U, 9U, 5U, 3U, 5U, 7U, 3U},
    vector<unsigned char>{0U, 7U, 8U, 0U, 1U, 7U, 1U, 5U, 7U},
    vector<unsigned char>{1U, 5U, 3U, 3U, 5U, 7U},
    vector<unsigned char>{9U, 7U, 8U, 9U, 5U, 7U, 10U, 1U, 2U},
    vector<unsigned char>{10U, 1U, 2U, 9U, 5U, 0U, 5U, 3U, 0U, 5U, 7U, 3U},
    vector<unsigned char>{8U, 0U, 2U, 8U, 2U, 5U, 8U, 5U, 7U, 10U, 5U, 2U},
    vector<unsigned char>{2U, 10U, 5U, 2U, 5U, 3U, 3U, 5U, 7U},
    vector<unsigned char>{7U, 9U, 5U, 7U, 8U, 9U, 3U, 11U, 2U},
    vector<unsigned char>{9U, 5U, 7U, 9U, 7U, 2U, 9U, 2U, 0U, 2U, 7U, 11U},
    vector<unsigned char>{2U, 3U, 11U, 0U, 1U, 8U, 1U, 7U, 8U, 1U, 5U, 7U},
    vector<unsigned char>{11U, 2U, 1U, 11U, 1U, 7U, 7U, 1U, 5U},
    vector<unsigned char>{9U, 5U, 8U, 8U, 5U, 7U, 10U, 1U, 3U, 10U, 3U, 11U},
    vector<unsigned char>{5U, 7U, 0U, 5U, 0U, 9U, 7U, 11U, 0U, 1U, 0U, 10U, 11U, 10U, 0U},
    vector<unsigned char>{11U, 10U, 0U, 11U, 0U, 3U, 10U, 5U, 0U, 8U, 0U, 7U, 5U, 7U, 0U},
    vector<unsigned char>{11U, 10U, 5U, 7U, 11U, 5U},
    vector<unsigned char>{10U, 6U, 5U},
    vector<unsigned char>{0U, 8U, 3U, 5U, 10U, 6U},
    vector<unsigned char>{9U, 0U, 1U, 5U, 10U, 6U},
    vector<unsigned char>{1U, 8U, 3U, 1U, 9U, 8U, 5U, 10U, 6U},
    vector<unsigned char>{1U, 6U, 5U, 2U, 6U, 1U},
    vector<unsigned char>{1U, 6U, 5U, 1U, 2U, 6U, 3U, 0U, 8U},
    vector<unsigned char>{9U, 6U, 5U, 9U, 0U, 6U, 0U, 2U, 6U},
    vector<unsigned char>{5U, 9U, 8U, 5U, 8U, 2U, 5U, 2U, 6U, 3U, 2U, 8U},
    vector<unsigned char>{2U, 3U, 11U, 10U, 6U, 5U},
    vector<unsigned char>{11U, 0U, 8U, 11U, 2U, 0U, 10U, 6U, 5U},
    vector<unsigned char>{0U, 1U, 9U, 2U, 3U, 11U, 5U, 10U, 6U},
    vector<unsigned char>{5U, 10U, 6U, 1U, 9U, 2U, 9U, 11U, 2U, 9U, 8U, 11U},
    vector<unsigned char>{6U, 3U, 11U, 6U, 5U, 3U, 5U, 1U, 3U},
    vector<unsigned char>{0U, 8U, 11U, 0U, 11U, 5U, 0U, 5U, 1U, 5U, 11U, 6U},
    vector<unsigned char>{3U, 11U, 6U, 0U, 3U, 6U, 0U, 6U, 5U, 0U, 5U, 9U},
    vector<unsigned char>{6U, 5U, 9U, 6U, 9U, 11U, 11U, 9U, 8U},
    vector<unsigned char>{5U, 10U, 6U, 4U, 7U, 8U},
    vector<unsigned char>{4U, 3U, 0U, 4U, 7U, 3U, 6U, 5U, 10U},
    vector<unsigned char>{1U, 9U, 0U, 5U, 10U, 6U, 8U, 4U, 7U},
    vector<unsigned char>{10U, 6U, 5U, 1U, 9U, 7U, 1U, 7U, 3U, 7U, 9U, 4U},
    vector<unsigned char>{6U, 1U, 2U, 6U, 5U, 1U, 4U, 7U, 8U},
    vector<unsigned char>{1U, 2U, 5U, 5U, 2U, 6U, 3U, 0U, 4U, 3U, 4U, 7U},
    vector<unsigned char>{8U, 4U, 7U, 9U, 0U, 5U, 0U, 6U, 5U, 0U, 2U, 6U},
    vector<unsigned char>{7U, 3U, 9U, 7U, 9U, 4U, 3U, 2U, 9U, 5U, 9U, 6U, 2U, 6U, 9U},
    vector<unsigned char>{3U, 11U, 2U, 7U, 8U, 4U, 10U, 6U, 5U},
    vector<unsigned char>{5U, 10U, 6U, 4U, 7U, 2U, 4U, 2U, 0U, 2U, 7U, 11U},
    vector<unsigned char>{0U, 1U, 9U, 4U, 7U, 8U, 2U, 3U, 11U, 5U, 10U, 6U},
    vector<unsigned char>{9U, 2U, 1U, 9U, 11U, 2U, 9U, 4U, 11U, 7U, 11U, 4U, 5U, 10U, 6U},
    vector<unsigned char>{8U, 4U, 7U, 3U, 11U, 5U, 3U, 5U, 1U, 5U, 11U, 6U},
    vector<unsigned char>{5U, 1U, 11U, 5U, 11U, 6U, 1U, 0U, 11U, 7U, 11U, 4U, 0U, 4U, 11U},
    vector<unsigned char>{0U, 5U, 9U, 0U, 6U, 5U, 0U, 3U, 6U, 11U, 6U, 3U, 8U, 4U, 7U},
    vector<unsigned char>{6U, 5U, 9U, 6U, 9U, 11U, 4U, 7U, 9U, 7U, 11U, 9U},
    vector<unsigned char>{10U, 4U, 9U, 6U, 4U, 10U},
    vector<unsigned char>{4U, 10U, 6U, 4U, 9U, 10U, 0U, 8U, 3U},
    vector<unsigned char>{10U, 0U, 1U, 10U, 6U, 0U, 6U, 4U, 0U},
    vector<unsigned char>{8U, 3U, 1U, 8U, 1U, 6U, 8U, 6U, 4U, 6U, 1U, 10U},
    vector<unsigned char>{1U, 4U, 9U, 1U, 2U, 4U, 2U, 6U, 4U},
    vector<unsigned char>{3U, 0U, 8U, 1U, 2U, 9U, 2U, 4U, 9U, 2U, 6U, 4U},
    vector<unsigned char>{0U, 2U, 4U, 4U, 2U, 6U},
    vector<unsigned char>{8U, 3U, 2U, 8U, 2U, 4U, 4U, 2U, 6U},
    vector<unsigned char>{10U, 4U, 9U, 10U, 6U, 4U, 11U, 2U, 3U},
    vector<unsigned char>{0U, 8U, 2U, 2U, 8U, 11U, 4U, 9U, 10U, 4U, 10U, 6U},
    vector<unsigned char>{3U, 11U, 2U, 0U, 1U, 6U, 0U, 6U, 4U, 6U, 1U, 10U},
    vector<unsigned char>{6U, 4U, 1U, 6U, 1U, 10U, 4U, 8U, 1U, 2U, 1U, 11U, 8U, 11U, 1U},
    vector<unsigned char>{9U, 6U, 4U, 9U, 3U, 6U, 9U, 1U, 3U, 11U, 6U, 3U},
    vector<unsigned char>{8U, 11U, 1U, 8U, 1U, 0U, 11U, 6U, 1U, 9U, 1U, 4U, 6U, 4U, 1U},
    vector<unsigned char>{3U, 11U, 6U, 3U, 6U, 0U, 0U, 6U, 4U},
    vector<unsigned char>{6U, 4U, 8U, 11U, 6U, 8U},
    vector<unsigned char>{7U, 10U, 6U, 7U, 8U, 10U, 8U, 9U, 10U},
    vector<unsigned char>{0U, 7U, 3U, 0U, 10U, 7U, 0U, 9U, 10U, 6U, 7U, 10U},
    vector<unsigned char>{10U, 6U, 7U, 1U, 10U, 7U, 1U, 7U, 8U, 1U, 8U, 0U},
    vector<unsigned char>{10U, 6U, 7U, 10U, 7U, 1U, 1U, 7U, 3U},
    vector<unsigned char>{1U, 2U, 6U, 1U, 6U, 8U, 1U, 8U, 9U, 8U, 6U, 7U},
    vector<unsigned char>{2U, 6U, 9U, 2U, 9U, 1U, 6U, 7U, 9U, 0U, 9U, 3U, 7U, 3U, 9U},
    vector<unsigned char>{7U, 8U, 0U, 7U, 0U, 6U, 6U, 0U, 2U},
    vector<unsigned char>{7U, 3U, 2U, 6U, 7U, 2U},
    vector<unsigned char>{2U, 3U, 11U, 10U, 6U, 8U, 10U, 8U, 9U, 8U, 6U, 7U},
    vector<unsigned char>{2U, 0U, 7U, 2U, 7U, 11U, 0U, 9U, 7U, 6U, 7U, 10U, 9U, 10U, 7U},
    vector<unsigned char>{1U, 8U, 0U, 1U, 7U, 8U, 1U, 10U, 7U, 6U, 7U, 10U, 2U, 3U, 11U},
    vector<unsigned char>{11U, 2U, 1U, 11U, 1U, 7U, 10U, 6U, 1U, 6U, 7U, 1U},
    vector<unsigned char>{8U, 9U, 6U, 8U, 6U, 7U, 9U, 1U, 6U, 11U, 6U, 3U, 1U, 3U, 6U},
    vector<unsigned char>{0U, 9U, 1U, 11U, 6U, 7U},
    vector<unsigned char>{7U, 8U, 0U, 7U, 0U, 6U, 3U, 11U, 0U, 11U, 6U, 0U},
    vector<unsigned char>{7U, 11U, 6U},
    vector<unsigned char>{7U, 6U, 11U},
    vector<unsigned char>{3U, 0U, 8U, 11U, 7U, 6U},
    vector<unsigned char>{0U, 1U, 9U, 11U, 7U, 6U},
    vector<unsigned char>{8U, 1U, 9U, 8U, 3U, 1U, 11U, 7U, 6U},
    vector<unsigned char>{10U, 1U, 2U, 6U, 11U, 7U},
    vector<unsigned char>{1U, 2U, 10U, 3U, 0U, 8U, 6U, 11U, 7U},
    vector<unsigned char>{2U, 9U, 0U, 2U, 10U, 9U, 6U, 11U, 7U},
    vector<unsigned char>{6U, 11U, 7U, 2U, 10U, 3U, 10U, 8U, 3U, 10U, 9U, 8U},
    vector<unsigned char>{7U, 2U, 3U, 6U, 2U, 7U},
    vector<unsigned char>{7U, 0U, 8U, 7U, 6U, 0U, 6U, 2U, 0U},
    vector<unsigned char>{2U, 7U, 6U, 2U, 3U, 7U, 0U, 1U, 9U},
    vector<unsigned char>{1U, 6U, 2U, 1U, 8U, 6U, 1U, 9U, 8U, 8U, 7U, 6U},
    vector<unsigned char>{10U, 7U, 6U, 10U, 1U, 7U, 1U, 3U, 7U},
    vector<unsigned char>{10U, 7U, 6U, 1U, 7U, 10U, 1U, 8U, 7U, 1U, 0U, 8U},
    vector<unsigned char>{0U, 3U, 7U, 0U, 7U, 10U, 0U, 10U, 9U, 6U, 10U, 7U},
    vector<unsigned char>{7U, 6U, 10U, 7U, 10U, 8U, 8U, 10U, 9U},
    vector<unsigned char>{6U, 8U, 4U, 11U, 8U, 6U},
    vector<unsigned char>{3U, 6U, 11U, 3U, 0U, 6U, 0U, 4U, 6U},
    vector<unsigned char>{8U, 6U, 11U, 8U, 4U, 6U, 9U, 0U, 1U},
    vector<unsigned char>{9U, 4U, 6U, 9U, 6U, 3U, 9U, 3U, 1U, 11U, 3U, 6U},
    vector<unsigned char>{6U, 8U, 4U, 6U, 11U, 8U, 2U, 10U, 1U},
    vector<unsigned char>{1U, 2U, 10U, 3U, 0U, 11U, 0U, 6U, 11U, 0U, 4U, 6U},
    vector<unsigned char>{4U, 11U, 8U, 4U, 6U, 11U, 0U, 2U, 9U, 2U, 10U, 9U},
    vector<unsigned char>{10U, 9U, 3U, 10U, 3U, 2U, 9U, 4U, 3U, 11U, 3U, 6U, 4U, 6U, 3U},
    vector<unsigned char>{8U, 2U, 3U, 8U, 4U, 2U, 4U, 6U, 2U},
    vector<unsigned char>{0U, 4U, 2U, 4U, 6U, 2U},
    vector<unsigned char>{1U, 9U, 0U, 2U, 3U, 4U, 2U, 4U, 6U, 4U, 3U, 8U},
    vector<unsigned char>{1U, 9U, 4U, 1U, 4U, 2U, 2U, 4U, 6U},
    vector<unsigned char>{8U, 1U, 3U, 8U, 6U, 1U, 8U, 4U, 6U, 6U, 10U, 1U},
    vector<unsigned char>{10U, 1U, 0U, 10U, 0U, 6U, 6U, 0U, 4U},
    vector<unsigned char>{4U, 6U, 3U, 4U, 3U, 8U, 6U, 10U, 3U, 0U, 3U, 9U, 10U, 9U, 3U},
    vector<unsigned char>{10U, 9U, 4U, 6U, 10U, 4U},
    vector<unsigned char>{4U, 9U, 5U, 7U, 6U, 11U},
    vector<unsigned char>{0U, 8U, 3U, 4U, 9U, 5U, 11U, 7U, 6U},
    vector<unsigned char>{5U, 0U, 1U, 5U, 4U, 0U, 7U, 6U, 11U},
    vector<unsigned char>{11U, 7U, 6U, 8U, 3U, 4U, 3U, 5U, 4U, 3U, 1U, 5U},
    vector<unsigned char>{9U, 5U, 4U, 10U, 1U, 2U, 7U, 6U, 11U},
    vector<unsigned char>{6U, 11U, 7U, 1U, 2U, 10U, 0U, 8U, 3U, 4U, 9U, 5U},
    vector<unsigned char>{7U, 6U, 11U, 5U, 4U, 10U, 4U, 2U, 10U, 4U, 0U, 2U},
    vector<unsigned char>{3U, 4U, 8U, 3U, 5U, 4U, 3U, 2U, 5U, 10U, 5U, 2U, 11U, 7U, 6U},
    vector<unsigned char>{7U, 2U, 3U, 7U, 6U, 2U, 5U, 4U, 9U},
    vector<unsigned char>{9U, 5U, 4U, 0U, 8U, 6U, 0U, 6U, 2U, 6U, 8U, 7U},
    vector<unsigned char>{3U, 6U, 2U, 3U, 7U, 6U, 1U, 5U, 0U, 5U, 4U, 0U},
    vector<unsigned char>{6U, 2U, 8U, 6U, 8U, 7U, 2U, 1U, 8U, 4U, 8U, 5U, 1U, 5U, 8U},
    vector<unsigned char>{9U, 5U, 4U, 10U, 1U, 6U, 1U, 7U, 6U, 1U, 3U, 7U},
    vector<unsigned char>{1U, 6U, 10U, 1U, 7U, 6U, 1U, 0U, 7U, 8U, 7U, 0U, 9U, 5U, 4U},
    vector<unsigned char>{4U, 0U, 10U, 4U, 10U, 5U, 0U, 3U, 10U, 6U, 10U, 7U, 3U, 7U, 10U},
    vector<unsigned char>{7U, 6U, 10U, 7U, 10U, 8U, 5U, 4U, 10U, 4U, 8U, 10U},
    vector<unsigned char>{6U, 9U, 5U, 6U, 11U, 9U, 11U, 8U, 9U},
    vector<unsigned char>{3U, 6U, 11U, 0U, 6U, 3U, 0U, 5U, 6U, 0U, 9U, 5U},
    vector<unsigned char>{0U, 11U, 8U, 0U, 5U, 11U, 0U, 1U, 5U, 5U, 6U, 11U},
    vector<unsigned char>{6U, 11U, 3U, 6U, 3U, 5U, 5U, 3U, 1U},
    vector<unsigned char>{1U, 2U, 10U, 9U, 5U, 11U, 9U, 11U, 8U, 11U, 5U, 6U},
    vector<unsigned char>{0U, 11U, 3U, 0U, 6U, 11U, 0U, 9U, 6U, 5U, 6U, 9U, 1U, 2U, 10U},
    vector<unsigned char>{11U, 8U, 5U, 11U, 5U, 6U, 8U, 0U, 5U, 10U, 5U, 2U, 0U, 2U, 5U},
    vector<unsigned char>{6U, 11U, 3U, 6U, 3U, 5U, 2U, 10U, 3U, 10U, 5U, 3U},
    vector<unsigned char>{5U, 8U, 9U, 5U, 2U, 8U, 5U, 6U, 2U, 3U, 8U, 2U},
    vector<unsigned char>{9U, 5U, 6U, 9U, 6U, 0U, 0U, 6U, 2U},
    vector<unsigned char>{1U, 5U, 8U, 1U, 8U, 0U, 5U, 6U, 8U, 3U, 8U, 2U, 6U, 2U, 8U},
    vector<unsigned char>{1U, 5U, 6U, 2U, 1U, 6U},
    vector<unsigned char>{1U, 3U, 6U, 1U, 6U, 10U, 3U, 8U, 6U, 5U, 6U, 9U, 8U, 9U, 6U},
    vector<unsigned char>{10U, 1U, 0U, 10U, 0U, 6U, 9U, 5U, 0U, 5U, 6U, 0U},
    vector<unsigned char>{0U, 3U, 8U, 5U, 6U, 10U},
    vector<unsigned char>{10U, 5U, 6U},
    vector<unsigned char>{11U, 5U, 10U, 7U, 5U, 11U},
    vector<unsigned char>{11U, 5U, 10U, 11U, 7U, 5U, 8U, 3U, 0U},
    vector<unsigned char>{5U, 11U, 7U, 5U, 10U, 11U, 1U, 9U, 0U},
    vector<unsigned char>{10U, 7U, 5U, 10U, 11U, 7U, 9U, 8U, 1U, 8U, 3U, 1U},
    vector<unsigned char>{11U, 1U, 2U, 11U, 7U, 1U, 7U, 5U, 1U},
    vector<unsigned char>{0U, 8U, 3U, 1U, 2U, 7U, 1U, 7U, 5U, 7U, 2U, 11U},
    vector<unsigned char>{9U, 7U, 5U, 9U, 2U, 7U, 9U, 0U, 2U, 2U, 11U, 7U},
    vector<unsigned char>{7U, 5U, 2U, 7U, 2U, 11U, 5U, 9U, 2U, 3U, 2U, 8U, 9U, 8U, 2U},
    vector<unsigned char>{2U, 5U, 10U, 2U, 3U, 5U, 3U, 7U, 5U},
    vector<unsigned char>{8U, 2U, 0U, 8U, 5U, 2U, 8U, 7U, 5U, 10U, 2U, 5U},
    vector<unsigned char>{9U, 0U, 1U, 5U, 10U, 3U, 5U, 3U, 7U, 3U, 10U, 2U},
    vector<unsigned char>{9U, 8U, 2U, 9U, 2U, 1U, 8U, 7U, 2U, 10U, 2U, 5U, 7U, 5U, 2U},
    vector<unsigned char>{1U, 3U, 5U, 3U, 7U, 5U},
    vector<unsigned char>{0U, 8U, 7U, 0U, 7U, 1U, 1U, 7U, 5U},
    vector<unsigned char>{9U, 0U, 3U, 9U, 3U, 5U, 5U, 3U, 7U},
    vector<unsigned char>{9U, 8U, 7U, 5U, 9U, 7U},
    vector<unsigned char>{5U, 8U, 4U, 5U, 10U, 8U, 10U, 11U, 8U},
    vector<unsigned char>{5U, 0U, 4U, 5U, 11U, 0U, 5U, 10U, 11U, 11U, 3U, 0U},
    vector<unsigned char>{0U, 1U, 9U, 8U, 4U, 10U, 8U, 10U, 11U, 10U, 4U, 5U},
    vector<unsigned char>{10U, 11U, 4U, 10U, 4U, 5U, 11U, 3U, 4U, 9U, 4U, 1U, 3U, 1U, 4U},
    vector<unsigned char>{2U, 5U, 1U, 2U, 8U, 5U, 2U, 11U, 8U, 4U, 5U, 8U},
    vector<unsigned char>{0U, 4U, 11U, 0U, 11U, 3U, 4U, 5U, 11U, 2U, 11U, 1U, 5U, 1U, 11U},
    vector<unsigned char>{0U, 2U, 5U, 0U, 5U, 9U, 2U, 11U, 5U, 4U, 5U, 8U, 11U, 8U, 5U},
    vector<unsigned char>{9U, 4U, 5U, 2U, 11U, 3U},
    vector<unsigned char>{2U, 5U, 10U, 3U, 5U, 2U, 3U, 4U, 5U, 3U, 8U, 4U},
    vector<unsigned char>{5U, 10U, 2U, 5U, 2U, 4U, 4U, 2U, 0U},
    vector<unsigned char>{3U, 10U, 2U, 3U, 5U, 10U, 3U, 8U, 5U, 4U, 5U, 8U, 0U, 1U, 9U},
    vector<unsigned char>{5U, 10U, 2U, 5U, 2U, 4U, 1U, 9U, 2U, 9U, 4U, 2U},
    vector<unsigned char>{8U, 4U, 5U, 8U, 5U, 3U, 3U, 5U, 1U},
    vector<unsigned char>{0U, 4U, 5U, 1U, 0U, 5U},
    vector<unsigned char>{8U, 4U, 5U, 8U, 5U, 3U, 9U, 0U, 5U, 0U, 3U, 5U},
    vector<unsigned char>{9U, 4U, 5U},
    vector<unsigned char>{4U, 11U, 7U, 4U, 9U, 11U, 9U, 10U, 11U},
    vector<unsigned char>{0U, 8U, 3U, 4U, 9U, 7U, 9U, 11U, 7U, 9U, 10U, 11U},
    vector<unsigned char>{1U, 10U, 11U, 1U, 11U, 4U, 1U, 4U, 0U, 7U, 4U, 11U},
    vector<unsigned char>{3U, 1U, 4U, 3U, 4U, 8U, 1U, 10U, 4U, 7U, 4U, 11U, 10U, 11U, 4U},
    vector<unsigned char>{4U, 11U, 7U, 9U, 11U, 4U, 9U, 2U, 11U, 9U, 1U, 2U},
    vector<unsigned char>{9U, 7U, 4U, 9U, 11U, 7U, 9U, 1U, 11U, 2U, 11U, 1U, 0U, 8U, 3U},
    vector<unsigned char>{11U, 7U, 4U, 11U, 4U, 2U, 2U, 4U, 0U},
    vector<unsigned char>{11U, 7U, 4U, 11U, 4U, 2U, 8U, 3U, 4U, 3U, 2U, 4U},
    vector<unsigned char>{2U, 9U, 10U, 2U, 7U, 9U, 2U, 3U, 7U, 7U, 4U, 9U},
    vector<unsigned char>{9U, 10U, 7U, 9U, 7U, 4U, 10U, 2U, 7U, 8U, 7U, 0U, 2U, 0U, 7U},
    vector<unsigned char>{3U, 7U, 10U, 3U, 10U, 2U, 7U, 4U, 10U, 1U, 10U, 0U, 4U, 0U, 10U},
    vector<unsigned char>{1U, 10U, 2U, 8U, 7U, 4U},
    vector<unsigned char>{4U, 9U, 1U, 4U, 1U, 7U, 7U, 1U, 3U},
    vector<unsigned char>{4U, 9U, 1U, 4U, 1U, 7U, 0U, 8U, 1U, 8U, 7U, 1U},
    vector<unsigned char>{4U, 0U, 3U, 7U, 4U, 3U},
    vector<unsigned char>{4U, 8U, 7U},
    vector<unsigned char>{9U, 10U, 8U, 10U, 11U, 8U},
    vector<unsigned char>{3U, 0U, 9U, 3U, 9U, 11U, 11U, 9U, 10U},
    vector<unsigned char>{0U, 1U, 10U, 0U, 10U, 8U, 8U, 10U, 11U},
    vector<unsigned char>{3U, 1U, 10U, 11U, 3U, 10U},
    vector<unsigned char>{1U, 2U, 11U, 1U, 11U, 9U, 9U, 11U, 8U},
    vector<unsigned char>{3U, 0U, 9U, 3U, 9U, 11U, 1U, 2U, 9U, 2U, 11U, 9U},
    vector<unsigned char>{0U, 2U, 11U, 8U, 0U, 11U},
    vector<unsigned char>{3U, 2U, 11U},
    vector<unsigned char>{2U, 3U, 8U, 2U, 8U, 10U, 10U, 8U, 9U},
    vector<unsigned char>{9U, 10U, 2U, 0U, 9U, 2U},
    vector<unsigned char>{2U, 3U, 8U, 2U, 8U, 10U, 0U, 1U, 8U, 1U, 10U, 8U},
    vector<unsigned char>{1U, 10U, 2U},
    vector<unsigned char>{1U, 3U, 8U, 9U, 1U, 8U},
    vector<unsigned char>{0U, 9U, 1U},
    vector<unsigned char>{0U, 3U, 8U},
    vector<unsigned char>{}
  };
  static const unsigned short ushEdgePattern[1U << GRID_VERTEX_COUNT_3D] = {0x0U, 0x109U, 0x203U, 0x30AU, 0x406U, 0x50FU, 0x605U, 0x70CU, 0x80CU, 0x905U, 0xA0FU, 0xB06U, 0xC0AU, 0xD03U, 0xE09U, 0xF00U, 0x190U, 0x99U, 0x393U, 0x29AU, 0x596U, 0x49FU, 0x795U, 0x69CU, 0x99CU, 0x895U, 0xB9FU, 0xA96U, 0xD9AU, 0xC93U, 0xF99U, 0xE90U, 0x230U, 0x339U, 0x33U, 0x13AU, 0x636U, 0x73FU, 0x435U, 0x53CU, 0xA3CU, 0xB35U, 0x83FU, 0x936U, 0xE3AU, 0xF33U, 0xC39U, 0xD30U, 0x3A0U, 0x2A9U, 0x1A3U, 0xAAU, 0x7A6U, 0x6AFU, 0x5A5U, 0x4ACU, 0xBACU, 0xAA5U, 0x9AFU, 0x8A6U, 0xFAAU, 0xEA3U, 0xDA9U, 0xCA0U, 0x460U, 0x569U, 0x663U, 0x76AU, 0x66U, 0x16FU, 0x265U, 0x36CU, 0xC6CU, 0xD65U, 0xE6FU, 0xF66U, 0x86AU, 0x963U, 0xA69U, 0xB60U, 0x5F0U, 0x4F9U, 0x7F3U, 0x6FAU, 0x1F6U, 0xFFU, 0x3F5U, 0x2FCU, 0xDFCU, 0xCF5U, 0xFFFU, 0xEF6U, 0x9FAU, 0x8F3U, 0xBF9U, 0xAF0U, 0x650U, 0x759U, 0x453U, 0x55AU, 0x256U, 0x35FU, 0x55U, 0x15CU, 0xE5CU, 0xF55U, 0xC5FU, 0xD56U, 0xA5AU, 0xB53U, 0x859U, 0x950U, 0x7C0U, 0x6C9U, 0x5C3U, 0x4CAU, 0x3C6U, 0x2CFU, 0x1C5U, 0xCCU, 0xFCCU, 0xEC5U, 0xDCFU, 0xCC6U, 0xBCAU, 0xAC3U, 0x9C9U, 0x8C0U, 0x8C0U, 0x9C9U, 0xAC3U, 0xBCAU, 0xCC6U, 0xDCFU, 0xEC5U, 0xFCCU, 0xCCU, 0x1C5U, 0x2CFU, 0x3C6U, 0x4CAU, 0x5C3U, 0x6C9U, 0x7C0U, 0x950U, 0x859U, 0xB53U, 0xA5AU, 0xD56U, 0xC5FU, 0xF55U, 0xE5CU, 0x15CU, 0x55U, 0x35FU, 0x256U, 0x55AU, 0x453U, 0x759U, 0x650U, 0xAF0U, 0xBF9U, 0x8F3U, 0x9FAU, 0xEF6U, 0xFFFU, 0xCF5U, 0xDFCU, 0x2FCU, 0x3F5U, 0xFFU, 0x1F6U, 0x6FAU, 0x7F3U, 0x4F9U, 0x5F0U, 0xB60U, 0xA69U, 0x963U, 0x86AU, 0xF66U, 0xE6FU, 0xD65U, 0xC6CU, 0x36CU, 0x265U, 0x16FU, 0x66U, 0x76AU, 0x663U, 0x569U, 0x460U, 0xCA0U, 0xDA9U, 0xEA3U, 0xFAAU, 0x8A6U, 0x9AFU, 0xAA5U, 0xBACU, 0x4ACU, 0x5A5U, 0x6AFU, 0x7A6U, 0xAAU, 0x1A3U, 0x2A9U, 0x3A0U, 0xD30U, 0xC39U, 0xF33U, 0xE3AU, 0x936U, 0x83FU, 0xB35U, 0xA3CU, 0x53CU, 0x435U, 0x73FU, 0x636U, 0x13AU, 0x33U, 0x339U, 0x230U, 0xE90U, 0xF99U, 0xC93U, 0xD9AU, 0xA96U, 0xB9FU, 0x895U, 0x99CU, 0x69CU, 0x795U, 0x49FU, 0x596U, 0x29AU, 0x393U, 0x99U, 0x190U, 0xF00U, 0xE09U, 0xD03U, 0xC0AU, 0xB06U, 0xA0FU, 0x905U, 0x80CU, 0x70CU, 0x605U, 0x50FU, 0x406U, 0x30AU, 0x203U, 0x109U, 0x0U};

  float*** pfGridValueMatrix = new float**[riDivisionCountX + 1];
  int* piTriangleCountArray = new int[rvecImplicitSurface.size()]();

  for (int i = 0; i <= riDivisionCountX; i++) {
    pfGridValueMatrix[i] = new float*[riDivisionCountY + 1];

    for (int j = 0; j <= riDivisionCountY; j++)
      pfGridValueMatrix[i][j] = new float[riDivisionCountZ + 1];
  }

  for (int i = 0; i < rvecImplicitSurface.size(); i++) {
    vector<Vertex3D> vecTransformedCorner{rvecImplicitSurface[i].MinimumCorner(), rvecImplicitSurface[i].MaximumCorner()};
    float fInversedTransformationMatrix[ELEMENT_COUNT_3D + 1][ELEMENT_COUNT_3D + 1];

    for (int j = 0; j <= ELEMENT_COUNT_3D; j++)
      for (int k = 0; k <= ELEMENT_COUNT_3D; k++)
        fInversedTransformationMatrix[j][k] = rpfTransformationMatrix[i][j][k];

    InverseTransformation<ELEMENT_COUNT_3D + 1>(fInversedTransformationMatrix);

    float Vertex3D::*pfVertex3DMember[ELEMENT_COUNT_3D] = {&Vertex3D::fX, &Vertex3D::fY, &Vertex3D::fZ};
    vector<Vertex3D> vecOriginalCorner{rvecImplicitSurface[i].MinimumCorner(), rvecImplicitSurface[i].MaximumCorner()};

    for (int j = 0; j < ELEMENT_COUNT_3D; j++) {
      vecTransformedCorner[INDICATOR_MINIMUM].*pfVertex3DMember[j] = vecTransformedCorner[INDICATOR_MAXIMUM].*pfVertex3DMember[j] = 0.0F;

      for (int k = 0; k < ELEMENT_COUNT_3D; k++)
        if (rpfTransformationMatrix[i][j][k] * vecOriginalCorner[INDICATOR_MINIMUM].*pfVertex3DMember[k] < rpfTransformationMatrix[i][j][k] * vecOriginalCorner[INDICATOR_MAXIMUM].*pfVertex3DMember[k]) {
          vecTransformedCorner[INDICATOR_MINIMUM].*pfVertex3DMember[j] += rpfTransformationMatrix[i][j][k] * vecOriginalCorner[INDICATOR_MINIMUM].*pfVertex3DMember[k];
          vecTransformedCorner[INDICATOR_MAXIMUM].*pfVertex3DMember[j] += rpfTransformationMatrix[i][j][k] * vecOriginalCorner[INDICATOR_MAXIMUM].*pfVertex3DMember[k];
        } else {
          vecTransformedCorner[INDICATOR_MINIMUM].*pfVertex3DMember[j] += rpfTransformationMatrix[i][j][k] * vecOriginalCorner[INDICATOR_MAXIMUM].*pfVertex3DMember[k];
          vecTransformedCorner[INDICATOR_MAXIMUM].*pfVertex3DMember[j] += rpfTransformationMatrix[i][j][k] * vecOriginalCorner[INDICATOR_MINIMUM].*pfVertex3DMember[k];
        }

      vecTransformedCorner[INDICATOR_MINIMUM].*pfVertex3DMember[j] += rpfTransformationMatrix[i][j][ELEMENT_COUNT_3D];
      vecTransformedCorner[INDICATOR_MAXIMUM].*pfVertex3DMember[j] += rpfTransformationMatrix[i][j][ELEMENT_COUNT_3D];
    }

    float fStep[ELEMENT_COUNT_3D] = {(vecTransformedCorner[INDICATOR_MAXIMUM].fX - vecTransformedCorner[INDICATOR_MINIMUM].fX) / static_cast<float>(riDivisionCountX), (vecTransformedCorner[INDICATOR_MAXIMUM].fY - vecTransformedCorner[INDICATOR_MINIMUM].fY) / static_cast<float>(riDivisionCountY), (vecTransformedCorner[INDICATOR_MAXIMUM].fZ - vecTransformedCorner[INDICATOR_MINIMUM].fZ) / static_cast<float>(riDivisionCountZ)};

    for (int j = 0; j <= riDivisionCountX; j++)
      for (int k = 0; k <= riDivisionCountY; k++)
        for (int l = 0; l <= riDivisionCountZ; l++) {
          Vertex3D Vertex3DSamplingPoint{vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * static_cast<float>(j), vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * static_cast<float>(k), vecTransformedCorner[INDICATOR_MINIMUM].fZ + fStep[INDICATOR_Z] * static_cast<float>(l)};

          pfGridValueMatrix[j][k][l] = rvecImplicitSurface[i].At(TransformCoordinate(fInversedTransformationMatrix, Vertex3DSamplingPoint));
        }

    unsigned char uchOffset[GRID_VERTEX_COUNT_3D][ELEMENT_COUNT_3D] = {
      {0U, 0U, 0U},
      {0U, 0U, 1U},
      {0U, 1U, 1U},
      {0U, 1U, 0U},
      {1U, 0U, 0U},
      {1U, 0U, 1U},
      {1U, 1U, 1U},
      {1U, 1U, 0U}
    };

    for (int j = 0; j < riDivisionCountX; j++)
      for (int k = 0; k < riDivisionCountY; k++)
        for (int l = 0; l < riDivisionCountZ; l++) {
          unsigned char uchFlag = 0U;

          do {
            unsigned char uchDigit = 1U;

            for (const unsigned char* const& rpuchPointer : uchOffset) {
              if (pfGridValueMatrix[rpuchPointer[INDICATOR_X] + i][rpuchPointer[INDICATOR_Y] + j][rpuchPointer[INDICATOR_Z] + k] >= 0.0F)
                uchFlag |= uchDigit;

              uchDigit <<= 1U;
            }
          } while (false);

          if (!ushEdgePattern[uchFlag])
            continue;

          vector<Vertex3D> vecInterpolationPoint(GRID_EDGE_COUNT_3D, Vertex3D{0.0F, 0.0F, 0.0F});
          
          do {
            unsigned char uchFactor = 0U;

            if (ushEdgePattern[uchFlag] & 1U << uchFactor++)
              vecInterpolationPoint[uchFactor - 1] = Vertex3D{vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * static_cast<float>(j), vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * static_cast<float>(k), (pfGridValueMatrix[j][k][l + 1] * (vecTransformedCorner[INDICATOR_MINIMUM].fZ + fStep[INDICATOR_Z] * static_cast<float>(l)) - pfGridValueMatrix[j][k][l] * (vecTransformedCorner[INDICATOR_MINIMUM].fZ + fStep[INDICATOR_Z] * (static_cast<float>(l) + 1.0F))) / (pfGridValueMatrix[j][k][l + 1] - pfGridValueMatrix[j][k][l])};

            if (ushEdgePattern[uchFlag] & 1U << uchFactor++)
              vecInterpolationPoint[uchFactor - 1] = Vertex3D{vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * static_cast<float>(j), (pfGridValueMatrix[j][k + 1][l + 1] * (vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * static_cast<float>(k)) - pfGridValueMatrix[j][k][l + 1] * (vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * (static_cast<float>(k) + 1.0F))) / (pfGridValueMatrix[j][k + 1][l + 1] - pfGridValueMatrix[j][k][l + 1]), vecTransformedCorner[INDICATOR_MINIMUM].fZ + fStep[INDICATOR_Z] * (static_cast<float>(l) + 1.0F)};

            if (ushEdgePattern[uchFlag] & 1U << uchFactor++)
              vecInterpolationPoint[uchFactor - 1] = Vertex3D{vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * static_cast<float>(j), vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * (static_cast<float>(k) + 1.0F), (pfGridValueMatrix[j][k + 1][l + 1] * (vecTransformedCorner[INDICATOR_MINIMUM].fZ + fStep[INDICATOR_Z] * static_cast<float>(l)) - pfGridValueMatrix[j][k + 1][l] * (vecTransformedCorner[INDICATOR_MINIMUM].fZ + fStep[INDICATOR_Z] * (static_cast<float>(l) + 1.0F))) / (pfGridValueMatrix[j][k + 1][l + 1] - pfGridValueMatrix[j][k + 1][l])};

            if (ushEdgePattern[uchFlag] & 1U << uchFactor++)
              vecInterpolationPoint[uchFactor - 1] = Vertex3D{vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * static_cast<float>(j), (pfGridValueMatrix[j][k + 1][l] * (vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * static_cast<float>(k)) - pfGridValueMatrix[j][k][l] * (vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * (static_cast<float>(k) + 1.0F))) / (pfGridValueMatrix[j][k + 1][l] - pfGridValueMatrix[j][k][l]), vecTransformedCorner[INDICATOR_MINIMUM].fZ + fStep[INDICATOR_Z] * static_cast<float>(l)};

            if (ushEdgePattern[uchFlag] & 1U << uchFactor++)
              vecInterpolationPoint[uchFactor - 1] = Vertex3D{vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * (static_cast<float>(j) + 1.0F), vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * static_cast<float>(k), (pfGridValueMatrix[j + 1][k][l + 1] * (vecTransformedCorner[INDICATOR_MINIMUM].fZ + fStep[INDICATOR_Z] * static_cast<float>(l)) - pfGridValueMatrix[j + 1][k][l] * (vecTransformedCorner[INDICATOR_MINIMUM].fZ + fStep[INDICATOR_Z] * (static_cast<float>(l) + 1.0F))) / (pfGridValueMatrix[j + 1][k][l + 1] - pfGridValueMatrix[j + 1][k][l])};

            if (ushEdgePattern[uchFlag] & 1U << uchFactor++)
              vecInterpolationPoint[uchFactor - 1] = Vertex3D{vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * (static_cast<float>(j) + 1.0F), (pfGridValueMatrix[j + 1][k + 1][l + 1] * (vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * static_cast<float>(k)) - pfGridValueMatrix[j + 1][k][l + 1] * (vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * (static_cast<float>(k) + 1.0F))) / (pfGridValueMatrix[j + 1][k + 1][l + 1] - pfGridValueMatrix[j + 1][k][l + 1]), vecTransformedCorner[INDICATOR_MINIMUM].fZ + fStep[INDICATOR_Z] * (static_cast<float>(l) + 1.0F)};

            if (ushEdgePattern[uchFlag] & 1U << uchFactor++)
              vecInterpolationPoint[uchFactor - 1] = Vertex3D{vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * (static_cast<float>(j) + 1.0F), vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * (static_cast<float>(k) + 1.0F), (pfGridValueMatrix[j + 1][k + 1][l + 1] * (vecTransformedCorner[INDICATOR_MINIMUM].fZ + fStep[INDICATOR_Z] * static_cast<float>(l)) - pfGridValueMatrix[j + 1][k + 1][l] * (vecTransformedCorner[INDICATOR_MINIMUM].fZ + fStep[INDICATOR_Z] * (static_cast<float>(l) + 1.0F))) / (pfGridValueMatrix[j + 1][k + 1][l + 1] - pfGridValueMatrix[j + 1][k + 1][l])};

            if (ushEdgePattern[uchFlag] & 1U << uchFactor++)
              vecInterpolationPoint[uchFactor - 1] = Vertex3D{vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * (static_cast<float>(j) + 1.0F), (pfGridValueMatrix[j + 1][k + 1][l] * (vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * static_cast<float>(k)) - pfGridValueMatrix[j + 1][k][l] * (vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * (static_cast<float>(k) + 1.0F))) / (pfGridValueMatrix[j + 1][k + 1][l] - pfGridValueMatrix[j + 1][k][l]), vecTransformedCorner[INDICATOR_MINIMUM].fZ + fStep[INDICATOR_Z] * static_cast<float>(l)};

            if (ushEdgePattern[uchFlag] & 1U << uchFactor++)
              vecInterpolationPoint[uchFactor - 1] = Vertex3D{(pfGridValueMatrix[j + 1][k][l] * (vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * static_cast<float>(j)) - pfGridValueMatrix[j][k][l] * (vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * (static_cast<float>(j) + 1.0F))) / (pfGridValueMatrix[j + 1][k][l] - pfGridValueMatrix[j][k][l]), vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * static_cast<float>(k), vecTransformedCorner[INDICATOR_MINIMUM].fZ + fStep[INDICATOR_Z] * static_cast<float>(l)};

            if (ushEdgePattern[uchFlag] & 1U << uchFactor++)
              vecInterpolationPoint[uchFactor - 1] = Vertex3D{(pfGridValueMatrix[j + 1][k][l + 1] * (vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * static_cast<float>(j)) - pfGridValueMatrix[j][k][l + 1] * (vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * (static_cast<float>(j) + 1.0F))) / (pfGridValueMatrix[j + 1][k][l + 1] - pfGridValueMatrix[j][k][l + 1]), vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * static_cast<float>(k), vecTransformedCorner[INDICATOR_MINIMUM].fZ + fStep[INDICATOR_Z] * (static_cast<float>(l) + 1.0F)};

            if (ushEdgePattern[uchFlag] & 1U << uchFactor++)
              vecInterpolationPoint[uchFactor - 1] = Vertex3D{(pfGridValueMatrix[j + 1][k + 1][l + 1] * (vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * static_cast<float>(j)) - pfGridValueMatrix[j][k + 1][l + 1] * (vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * (static_cast<float>(j) + 1.0F))) / (pfGridValueMatrix[j + 1][k + 1][l + 1] - pfGridValueMatrix[j][k + 1][l + 1]), vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * (static_cast<float>(k) + 1.0F), vecTransformedCorner[INDICATOR_MINIMUM].fZ + fStep[INDICATOR_Z] * (static_cast<float>(l) + 1.0F)};

            if (ushEdgePattern[uchFlag] & 1U << uchFactor++)
              vecInterpolationPoint[uchFactor - 1] = Vertex3D{(pfGridValueMatrix[j + 1][k + 1][l] * (vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * static_cast<float>(j)) - pfGridValueMatrix[j][k + 1][l] * (vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * (static_cast<float>(j) + 1.0F))) / (pfGridValueMatrix[j + 1][k + 1][l] - pfGridValueMatrix[j][k + 1][l]), vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * (static_cast<float>(k) + 1.0F), vecTransformedCorner[INDICATOR_MINIMUM].fZ + fStep[INDICATOR_Z] * static_cast<float>(l)};
          } while (false);

          for (int m = 0; m * TRIANGLE_ELEMENT_COUNT < vecTrianglePattern[uchFlag].size(); m++) {
            for (int n = 0; n < TRIANGLE_ELEMENT_COUNT; n++)
              rofsFile << VERTEX_CHARACTER << ' ' << vecInterpolationPoint[vecTrianglePattern[uchFlag][m * TRIANGLE_ELEMENT_COUNT + n]].fX << ' ' << vecInterpolationPoint[vecTrianglePattern[uchFlag][m * TRIANGLE_ELEMENT_COUNT + n]].fY << ' ' << vecInterpolationPoint[vecTrianglePattern[uchFlag][m * TRIANGLE_ELEMENT_COUNT + n]].fZ << endl;

            piTriangleCountArray[i]++;
          }
        }

    if (rvecImplicitSurface.size() - i - 1)
      rofsFile << COMMENT_CHARACTER << endl;
  }

  for (int i = 0; i <= riDivisionCountX; i++) {
    for (int j = 0; j <= riDivisionCountY; j++)
      delete[] pfGridValueMatrix[i][j];

    delete[] pfGridValueMatrix[i];
  }

  delete[] pfGridValueMatrix;

  int iVertexCount = 0;

  for (int i = 0; i < rvecImplicitSurface.size(); i++) {
    for (int j = 0; j < piTriangleCountArray[i]; j++) {
      rofsFile << FACE_ELEMENT_CHARACTER << ' ';

      for (int k = 0; k < TRIANGLE_ELEMENT_COUNT; k++) {
        rofsFile << iVertexCount + k + 1 << DELIMITER_CHARACTER << DELIMITER_CHARACTER;

        if (TRIANGLE_ELEMENT_COUNT - k - 1)
          rofsFile << ' ';
        else
          rofsFile << endl;
      }

      iVertexCount += TRIANGLE_ELEMENT_COUNT;
    }

    if (rvecImplicitSurface.size() - i - 1)
      rofsFile << COMMENT_CHARACTER << endl;
  }

  delete[] piTriangleCountArray;
}