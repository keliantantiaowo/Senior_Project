#include <iostream>
#include <iomanip>
#include <chrono>
#include "Import File.h"
#include "Export File.h"
#include "Generate Implicit Surface.h"
#include "Collision Detection.h"

using namespace std;

string strFileName;
const char* pchFileNameList[] = {"../データ/三角形", "../データ/正方形", "../データ/円", "../データ/星 2D", "../データ/T型 2D", "../データ/十字型 2D", "../データ/ぎざぎざ円"};
// const char* pchFileNameList[] = {"../データ/三角錐", "../データ/キューブ", "../データ/球", "../データ/星 3D", "../データ/T型 3D", "../データ/十字型 3D", "../データ/大星型十二面体"};
const int iPolygonCount = GET_ARRAY_SIZE(pchFileNameList), iTestCount = 10E4, iInterval = 100;

int main(void)
{
  vector<vector<vector<Vertex2D>>> vecVertex(iPolygonCount, vector<vector<Vertex2D>>(0));
  vector<vector<Vertex2D>> vecTransformedVertex(iPolygonCount, vector<Vertex2D>(0)), vecVertexNormal(iPolygonCount, vector<Vertex2D>(0));
  vector<Mesh2D> vecPolygonSet(iPolygonCount, Mesh2D(vector<vector<Vertex2D>>(0))), vecTransformedPolygonSet(iPolygonCount, Mesh2D(vector<vector<Vertex2D>>(0)));
  vector<PiecewisePolynomial2D> vecImplicitSurface{};

  for (int i = 0; i < iPolygonCount; i++) {
    ImportObjectFile2D(pchFileNameList[i], vecVertex[i]);

    CalculateVertexNormal2D(*vecVertex[i].begin(), vecVertexNormal[i], Normalized);

    vecPolygonSet[i] = Mesh2D(vecVertex[i]);

    vecImplicitSurface.push_back(PiecewisePolynomial2D(*vecVertex[i].begin(), vecVertexNormal[i], 1000, 1000, 0.5F));
  }

  chrono::microseconds othConventionalMethodTime(0), othProposedMethodTime(0);
  float fTransformationMatrix[iPolygonCount][ELEMENT_COUNT_2D + 1][ELEMENT_COUNT_2D + 1];
  int iExecutionCount = 0, iPositiveCount = 0, iNegativeCount = 0, iFalsePositiveCount = 0, iFalseNegativeCount = 0;

  for (int i = 0; i < iTestCount; i++) {
    for (int j = 0; j < iPolygonCount; j++) {
      RandomizeTransformation(fTransformationMatrix[j], Vertex2D{-2.0F, -2.0F}, Vertex2D{2.0F, 2.0F});

      TransformCoordinate(fTransformationMatrix[j], *vecVertex[j].begin(), vecTransformedVertex[j]);

      (vecTransformedPolygonSet[j] = vecPolygonSet[j]).ApplyTransformation(fTransformationMatrix[j]);
    }

    for (int j = 0; j < iPolygonCount - 1; j++)
      for (int k = j + 1; k < iPolygonCount; k++) {
        chrono::high_resolution_clock::time_point othStart, othEnd;
        bool bDetectionResult[DETECTION_METHOD_COUNT];

        othStart = chrono::high_resolution_clock::now();
        bDetectionResult[INDICATOR_CONVENTIONAL_METHOD] = DetectCollision2D(vecTransformedPolygonSet[j], vecTransformedPolygonSet[k]);
        othEnd = chrono::high_resolution_clock::now();
        othConventionalMethodTime += chrono::duration_cast<chrono::microseconds>(othEnd - othStart);
        othStart = chrono::high_resolution_clock::now();
        bDetectionResult[INDICATOR_PROPOSED_METHOD] = DetectCollision2D(vecTransformedVertex[j], vecTransformedVertex[k], vecImplicitSurface[j], vecImplicitSurface[k], fTransformationMatrix[j], fTransformationMatrix[k]);
        othEnd = chrono::high_resolution_clock::now();
        othProposedMethodTime += chrono::duration_cast<chrono::microseconds>(othEnd - othStart);
        iExecutionCount++;

        if (bDetectionResult[INDICATOR_CONVENTIONAL_METHOD]) {
          iPositiveCount++;

          if (bDetectionResult[INDICATOR_CONVENTIONAL_METHOD] == bDetectionResult[INDICATOR_PROPOSED_METHOD])
            ;
          else {
            if (!(iFalseNegativeCount % iInterval) && j <= 2 && k > 2) {
              ostringstream ossData("../フォルスネガティブ");

              ossData.seekp(0, ios::end);

              ossData << ' ' << setw(DIGIT) << setfill(FILL_CHARACTER) << iExecutionCount << ' ';

              ExportObjectFile2D(ossData.str() + "メッシュ " + to_string(j), vecTransformedPolygonSet[j]);
              ExportDataFile2D(ossData.str() + "陰関数曲面 " + to_string(j), vecImplicitSurface[j], fTransformationMatrix[j], 1000, 1000);
              ExportObjectFile2D(ossData.str() + "メッシュ " + to_string(k), vecTransformedPolygonSet[k]);
              ExportDataFile2D(ossData.str() + "陰関数曲面 " + to_string(k), vecImplicitSurface[k], fTransformationMatrix[k], 1000, 1000);
            }

            iFalseNegativeCount++;
          }
        } else {
          iNegativeCount++;

          if (bDetectionResult[INDICATOR_CONVENTIONAL_METHOD] == bDetectionResult[INDICATOR_PROPOSED_METHOD])
            ;
          else {
            ostringstream ossData("../フォルスポジティブ");

            ossData.seekp(0, ios::end);

            ossData << ' ' << setw(DIGIT) << setfill(FILL_CHARACTER) << iExecutionCount << ' ';

            ExportObjectFile2D(ossData.str() + "メッシュ " + to_string(j), vecTransformedPolygonSet[j]);
            ExportDataFile2D(ossData.str() + "陰関数曲面 " + to_string(j), vecImplicitSurface[j], fTransformationMatrix[j], 1000, 1000);
            ExportObjectFile2D(ossData.str() + "メッシュ " + to_string(k), vecTransformedPolygonSet[k]);
            ExportDataFile2D(ossData.str() + "陰関数曲面 " + to_string(k), vecImplicitSurface[k], fTransformationMatrix[k], 1000, 1000);

            iFalsePositiveCount++;
          }
        }
      }
  }

  cout << "iExecutionCount = " << iExecutionCount << ", iPositiveCount = " << iPositiveCount << ", iNegativeCount = " << iNegativeCount << ", iFalsePositiveCount = " << iFalsePositiveCount << ", iFalseNegativeCount = " << iFalseNegativeCount << endl << "othConventionalMethodTime = " << othConventionalMethodTime.count() << "µs, othProposedMethodTime = " << othProposedMethodTime.count() << "µs" << endl;

/*
  vector<vector<vector<Vertex3D>>> vecVertex(iPolygonCount, vector<vector<Vertex3D>>(0)), vecSurfaceNormal(iPolygonCount, vector<vector<Vertex3D>>(0));
  vector<vector<Vertex3D>> vecTransformedVertex(iPolygonCount, vector<Vertex3D>(0)), vecVertexNormal(iPolygonCount, vector<Vertex3D>(0));
  vector<vector<vector<vector<unsigned>>>> vecTriangleVertexIndex(iPolygonCount, vector<vector<vector<unsigned>>>(0)), vecTriangleNormalIndex(iPolygonCount, vector<vector<vector<unsigned>>>(0));
  vector<Mesh3D> vecPolygonSet(iPolygonCount, Mesh3D(vector<vector<Vertex3D>>(0))), vecTransformedPolygonSet(iPolygonCount, Mesh3D(vector<vector<Vertex3D>>(0)));
  vector<PiecewisePolynomial3D> vecImplicitSurface{};

  for (int i = 0; i < iPolygonCount; i++) {
    ImportObjectFile3D(pchFileNameList[i], vecVertex[i], vecSurfaceNormal[i], vecTriangleVertexIndex[i], vecTriangleNormalIndex[i]);

    CalculateVertexNormal3D(*vecVertex[i].begin(), vecVertexNormal[i], *vecSurfaceNormal[i].begin(), *vecTriangleVertexIndex[i].begin(), *vecTriangleNormalIndex[i].begin(), Normalized);

    vecPolygonSet[i] = Mesh3D(vecVertex[i]);

    vecImplicitSurface.push_back(PiecewisePolynomial3D(*vecVertex[i].begin(), vecVertexNormal[i], *vecSurfaceNormal[i].begin(), *vecTriangleVertexIndex[i].begin(), *vecTriangleNormalIndex[i].begin(), 200, 200, 200, 0.5F));
  }
  
  float fTransformationMatrix[iPolygonCount][ELEMENT_COUNT_3D + 1][ELEMENT_COUNT_3D + 1];
  chrono::microseconds othConventionalMethodTime(0), othProposedMethodTime(0);
  int iExecutionCount = 0, iPositiveCount = 0, iNegativeCount = 0, iFalsePositiveCount = 0, iFalseNegativeCount = 0;

  for (int i = 0; i < iTestCount; i++) {
    for (int j = 0; j < iPolygonCount; j++) {
      RandomizeTransformation(fTransformationMatrix[j], Vertex3D{-2.0F, -2.0F, -2.0F}, Vertex3D{2.0F, 2.0F, 2.0F});

      TransformCoordinate(fTransformationMatrix[j], *vecVertex[j].begin(), vecTransformedVertex[j]);

      (vecTransformedPolygonSet[j] = vecPolygonSet[j]).ApplyTransformation(fTransformationMatrix[j]);
    }

    for (int j = 0; j < iPolygonCount - 1; j++)
      for (int k = j + 1; k < iPolygonCount; k++) {
        chrono::high_resolution_clock::time_point othStart, othEnd;
        bool bDetectionResult[DETECTION_METHOD_COUNT];

        othStart = chrono::high_resolution_clock::now();
        bDetectionResult[INDICATOR_CONVENTIONAL_METHOD] = DetectCollision3D(vecTransformedPolygonSet[j], vecTransformedPolygonSet[k]);
        othEnd = chrono::high_resolution_clock::now();
        othConventionalMethodTime += chrono::duration_cast<chrono::microseconds>(othEnd - othStart);
        othStart = chrono::high_resolution_clock::now();
        bDetectionResult[INDICATOR_PROPOSED_METHOD] = DetectCollision3D(vecTransformedVertex[j], vecTransformedVertex[k], vecImplicitSurface[j], vecImplicitSurface[k], fTransformationMatrix[j], fTransformationMatrix[k]);
        othEnd = chrono::high_resolution_clock::now();
        othProposedMethodTime += chrono::duration_cast<chrono::microseconds>(othEnd - othStart);
        iExecutionCount++;

        if (bDetectionResult[INDICATOR_CONVENTIONAL_METHOD]) {
          iPositiveCount++;

          if (bDetectionResult[INDICATOR_CONVENTIONAL_METHOD] == bDetectionResult[INDICATOR_PROPOSED_METHOD])
            ;
          else {
            if (!(iFalseNegativeCount % iInterval)) {
              ostringstream ossData("../フォルスネガティブ");

              ossData.seekp(0, ios::end);

              ossData << ' ' << setw(DIGIT) << setfill(FILL_CHARACTER) << iFalseNegativeCount + 1 << ' ';

              ExportObjectFile3D(ossData.str() + "メッシュ " + to_string(j), vecTransformedPolygonSet[j]);
              ExportObjectFile3D(ossData.str() + "陰関数曲面 " + to_string(j), vecImplicitSurface[j], fTransformationMatrix[j], 400, 400, 400);
              ExportObjectFile3D(ossData.str() + "メッシュ " + to_string(k), vecTransformedPolygonSet[k]);
              ExportObjectFile3D(ossData.str() + "陰関数曲面 " + to_string(k), vecImplicitSurface[k], fTransformationMatrix[k], 400, 400, 400);
            }

            iFalseNegativeCount++;
          }
        } else {
          iNegativeCount++;

          if (bDetectionResult[INDICATOR_CONVENTIONAL_METHOD] == bDetectionResult[INDICATOR_PROPOSED_METHOD])
            ;
          else {
            ostringstream ossData("../フォルスポジティブ");

            ossData.seekp(0, ios::end);

            ossData << ' ' << setw(DIGIT) << setfill(FILL_CHARACTER) << iFalsePositiveCount + 1 << ' ';

            ExportObjectFile3D(ossData.str() + "メッシュ " + to_string(j), vecTransformedPolygonSet[j]);
            ExportObjectFile3D(ossData.str() + "陰関数曲面 " + to_string(j), vecImplicitSurface[j], fTransformationMatrix[j], 400, 400, 400);
            ExportObjectFile3D(ossData.str() + "メッシュ " + to_string(k), vecTransformedPolygonSet[k]);
            ExportObjectFile3D(ossData.str() + "陰関数曲面 " + to_string(k), vecImplicitSurface[k], fTransformationMatrix[k], 400, 400, 400);

            iFalsePositiveCount++;
          }
        }
      }
  }

  cout << "iExecutionCount = " << iExecutionCount << ", iPositiveCount = " << iPositiveCount << ", iNegativeCount = " << iNegativeCount << ", iFalsePositiveCount = " << iFalsePositiveCount << ", iFalseNegativeCount = " << iFalseNegativeCount << endl << "othConventionalMethodTime = " << othConventionalMethodTime.count() << "µs, othProposedMethodTime = " << othProposedMethodTime.count() << "µs" << endl;
*/

/*
  vector<Vertex2D> vecVertex(0), vecVertexNormal(0);

  for (int i = 0; i < iPolygonCount; i++) {
    ImportObjectFile2D(strFileName = pchFileNameList[i], vecVertex);

    CalculateVertexNormal2D(vecVertex, vecVertexNormal, Normalized);

    PiecewisePolynomial2D PiecewisePolynomial2DData(vecVertex, vecVertexNormal, 1000, 1000, 0.5F);

    ExportObjectFile2D(strFileName + " 陰関数曲面", PiecewisePolynomial2DData, nullptr, 1000, 1000);
    ExportDataFile2D(strFileName + " 陰関数曲面", PiecewisePolynomial2DData, nullptr, 1000, 1000);
  }
*/

/*
  vector<Vertex3D> vecVertex(0), vecSurfaceNormal(0), vecVertexNormal(0);
  vector<vector<unsigned>> vecTriangleVertexIndex(0), vecTriangleNormalIndex(0);

  for (int i = 0; i < iPolygonCount; i++) {
    ImportObjectFile3D(strFileName = pchFileNameList[i], vecVertex, vecSurfaceNormal, vecTriangleVertexIndex, vecTriangleNormalIndex);

    CalculateVertexNormal3D(vecVertex, vecVertexNormal, vecSurfaceNormal, vecTriangleVertexIndex, vecTriangleNormalIndex, Normalized);

    PiecewisePolynomial3D PiecewisePolynomial3DData(vecVertex, vecVertexNormal, vecSurfaceNormal, vecTriangleVertexIndex, vecTriangleNormalIndex, 200, 200, 200, 0.5F);

    ExportObjectFile3D(strFileName + " 陰関数曲面", PiecewisePolynomial3DData, nullptr, 400, 400, 400);
  }
*/
}