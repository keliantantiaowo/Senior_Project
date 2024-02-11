#include <sstream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#if defined(_WIN32) || defined(_WIN64)
#include <direct.h>
#elif defined(__linux) || defined(__MACH__)
#include <unistd.h>
#endif
#include <sys/stat.h>
#include "Definition.h"
#include "Display Result.h"
#include "Export File.h"
#include "Geometric Operation.h"
#include "Generate Implicit Surface.h"
#include "Marching Square & Cube.h"

using namespace std;

void ExportObjectFile2D(string strFileName, const ConvexHull2D& rConvexHull2DData)
{
  ofstream ofsFile;

  if (strFileName.empty()) {
    cout << "書き出すファイル名を入力してください：" << flush;

    getline(cin, strFileName);

    cout << endl;
  }

  ofsFile.open(strFileName + OBJECT_FILE_EXTENSION);

  if (!ofsFile.is_open()) {
    cerr << "ファイル「" << strFileName << OBJECT_FILE_EXTENSION << "」を開けませんでした。" << endl;

    exit(EXIT_FAILURE);
  }

  for (list<Vertex2D>::const_iterator i = rConvexHull2DData.Begin(); i != rConvexHull2DData.End(); i++)
    ofsFile << VERTEX_CHARACTER << ' ' << i->fX << ' ' << i->fY << ' ' << 0.0F << endl;

  ofsFile << FACE_ELEMENT_CHARACTER << ' ';

  for (int i = 0; i < rConvexHull2DData.VertexCount(); i++) {
    ofsFile << i + 1 << DELIMITER_CHARACTER << DELIMITER_CHARACTER;

    if (rConvexHull2DData.VertexCount() - i - 1)
      ofsFile << ' ';
    else
      ofsFile << endl;
  }

  ofsFile.close();
}

void ExportObjectFile2D(string strFileName, const Mesh2D& rMesh2DData)
{
  ofstream ofsFile;

  if (strFileName.empty()) {
    cout << "書き出すファイル名を入力してください；" << flush;

    getline(cin, strFileName);

    cout << endl;
  }

  ofsFile.open(strFileName + OBJECT_FILE_EXTENSION);

  if (!ofsFile.is_open()) {
    cerr << "ファイル「" << strFileName << OBJECT_FILE_EXTENSION << "」を開けませんでした。" << endl;

    exit(EXIT_FAILURE);
  }

  for (int i = 0; i < rMesh2DData.Size(); i++)
    for (list<Vertex2D>::const_iterator j = rMesh2DData[i].Begin(); j != rMesh2DData[i].End(); j++)
      ofsFile << VERTEX_CHARACTER << ' ' << j->fX << ' ' << j->fY << ' ' << 0.0F << endl;

  int iCount = 0;

  for (int i = 0; i < rMesh2DData.Size(); i++) {
    ofsFile << FACE_ELEMENT_CHARACTER << ' ';

    for (int j = 0; j < rMesh2DData[i].VertexCount(); j++) {
      ofsFile << iCount + 1 << DELIMITER_CHARACTER << DELIMITER_CHARACTER;

      iCount++;

      if (rMesh2DData[i].VertexCount() - j - 1)
        ofsFile << ' ';
      else
        ofsFile << endl;
    }
  }
}

void ExportObjectFile2D(string strFileName, const vector<Mesh2D>& rvecPolygonSet)
{
  ofstream ofsFile;

  if (strFileName.empty()) {
    cout << "書き出すファイル名を入力してください：" << flush;

    getline(cin, strFileName);

    cout << endl;
  }

  ofsFile.open(strFileName + OBJECT_FILE_EXTENSION);

  if (!ofsFile.is_open()) {
    cerr << "ファイル「" << strFileName << OBJECT_FILE_EXTENSION << "」を開けませんでした。" << endl;

    exit(EXIT_FAILURE);
  }

  for (int i = 0; i < rvecPolygonSet.size(); i++)
    for (int j = 0; j < rvecPolygonSet[i].Size(); j++)
      for (list<Vertex2D>::const_iterator k = rvecPolygonSet[i][j].Begin(); k != rvecPolygonSet[i][j].End(); k++)
        ofsFile << VERTEX_CHARACTER << ' ' << k->fX << ' ' << k->fY << ' ' << 0.0F << endl;

  int iCount = 0;

  for (int i = 0; i < rvecPolygonSet.size(); i++)
    for (int j = 0; j < rvecPolygonSet[i].Size(); j++) {
      ofsFile << FACE_ELEMENT_CHARACTER << ' ';

      for (int k = 0; k < rvecPolygonSet[i][j].VertexCount(); k++) {
        ofsFile << iCount + 1 << DELIMITER_CHARACTER << DELIMITER_CHARACTER;

        iCount++;

        if (rvecPolygonSet[i][j].VertexCount() - k - 1)
          ofsFile << ' ';
        else
          ofsFile << endl;
      }
    }

  ofsFile.close();
}

void ExportObjectFile2D(string strFileName, const PiecewisePolynomial2D& rPiecewisePolynomial2DData, const float (*const& rpfTransformationMatrix)[ELEMENT_COUNT_2D + 1], const int& riDivisionCountX, const int& riDivisionCountY)
{
  ofstream ofsFile;

  if (strFileName.empty()) {
    cout << "書き出すファイル名を入力してください：" << flush;

    getline(cin, strFileName);
  }

  ofsFile.open(strFileName + OBJECT_FILE_EXTENSION);

  if (!ofsFile.is_open()) {
    cerr << "ファイル「" << strFileName << OBJECT_FILE_EXTENSION << "」を開けませんでした。" << endl;

    exit(EXIT_FAILURE);
  }

  MarchingSquare(ofsFile, rPiecewisePolynomial2DData, rpfTransformationMatrix, riDivisionCountX, riDivisionCountY);

  ofsFile.close();
}

void ExportObjectFile2D(string strFileName, const vector<PiecewisePolynomial2D>& rvecImplicitSurface, const float (*const& rpfTransformationMatrix)[ELEMENT_COUNT_2D + 1][ELEMENT_COUNT_2D + 1], const int& riDivisionCountX, const int& riDivisionCountY)
{
  ofstream ofsFile;

  if (strFileName.empty()) {
    cout << "書き出すファイル名を入力してください：" << flush;

    getline(cin, strFileName);
  }

  ofsFile.open(strFileName + OBJECT_FILE_EXTENSION);

  if (!ofsFile.is_open()) {
    cerr << "ファイル「" << strFileName << OBJECT_FILE_EXTENSION << "」を開けませんでした。" << endl;

    exit(EXIT_FAILURE);
  }

  MarchingSquare(ofsFile, rvecImplicitSurface, rpfTransformationMatrix, riDivisionCountX, riDivisionCountY);

  ofsFile.close();
}

void ExportObjectFile3D(string strFileName, const ConvexHull3D& rConvexHull3DData)
{
  ofstream ofsFile;

  if (strFileName.empty()) {
    cout << "書き出すファイル名を入力してください：" << flush;

    getline(cin, strFileName);

    cout << endl;
  }

  ofsFile.open(strFileName + OBJECT_FILE_EXTENSION);

  if (!ofsFile.is_open()) {
    cerr << "ファイル「" << strFileName << OBJECT_FILE_EXTENSION << "」を開けませんでした。" << endl;

    exit(EXIT_FAILURE);
  }

  for (list<Vertex3D>::const_iterator i = rConvexHull3DData.VertexBegin(); i != rConvexHull3DData.VertexEnd(); i++)
    ofsFile << VERTEX_CHARACTER << ' ' << i->fX << ' ' << i->fY << ' ' << i->fZ << endl;

  switch (rConvexHull3DData.ShapeState()) {
  case Plane :
    ofsFile << FACE_ELEMENT_CHARACTER << ' ';

    for (int i = 0; i < rConvexHull3DData.VertexCount(); i++) {
      ofsFile << i + 1 << DELIMITER_CHARACTER << DELIMITER_CHARACTER;

      if (rConvexHull3DData.VertexCount() - i - 1)
        ofsFile << ' ';
      else
        ofsFile << endl;
    }

    break;
  case Solid :
    do {
      for (list<Surface3D>::const_iterator i = rConvexHull3DData.SurfaceBegin(); i != rConvexHull3DData.SurfaceEnd(); i++)
        ofsFile << VERTEX_CHARACTER << NORMAL_CHARACTER << ' ' << i->Vertex3DSurfaceNormal.fX << ' ' << i->Vertex3DSurfaceNormal.fY << ' ' << i->Vertex3DSurfaceNormal.fZ << endl;

      int iVertexCount = 0, iNormalCount = 0;

      for (list<Surface3D>::const_iterator i = rConvexHull3DData.SurfaceBegin(); i != rConvexHull3DData.SurfaceEnd(); i++) {
        ofsFile << FACE_ELEMENT_CHARACTER << ' ';

        const HalfEdge3D* pHalfEdge3DAnchor = i->pHalfEdge3DPointer;

        do {
          list<Vertex3D>::const_iterator iterPosition = rConvexHull3DData.VertexBegin();
          const Vertex3D* pVertex3DAnchor = pHalfEdge3DAnchor->pVertex3DPointer;

          iVertexCount = 0;

          do {
            if (&*iterPosition == pVertex3DAnchor)
              break;

            iterPosition++;
            iVertexCount++;
          } while (iterPosition != rConvexHull3DData.VertexEnd());

          ofsFile << iVertexCount + 1 << DELIMITER_CHARACTER << DELIMITER_CHARACTER << iNormalCount + 1;

          if ((pHalfEdge3DAnchor = pHalfEdge3DAnchor->pHalfEdge3DNext) == i->pHalfEdge3DPointer) {
            ofsFile << endl;

            break;
          } else
            ofsFile << ' ';
        } while (true);

        iNormalCount++;
      }
    } while (false);

    break;
  case Undefined :
    cerr << "凸包が正しく生成されていません。" << endl;

    exit(EXIT_FAILURE);
  default : break;
  }

  ofsFile.close();
}

void ExportObjectFile3D(string strFileName, const Mesh3D& rMesh3DData)
{
  ofstream ofsFile;

  if (strFileName.empty()) {
    cout << "書き出すファイル名を入力してください：" << flush;

    getline(cin, strFileName);

    cout << endl;
  }

  ofsFile.open(strFileName + OBJECT_FILE_EXTENSION);

  if (!ofsFile.is_open()) {
    cerr << "ファイル「" << strFileName << OBJECT_FILE_EXTENSION << "」を開けませんでした。" << endl;

    exit(EXIT_FAILURE);
  }

  for (int i = 0; i < rMesh3DData.Size(); i++)
    for (list<Vertex3D>::const_iterator j = rMesh3DData[i].VertexBegin(); j != rMesh3DData[i].VertexEnd(); j++)
      ofsFile << VERTEX_CHARACTER << ' ' << j->fX << ' ' << j->fY << ' ' << j->fZ << endl;

  for (int i = 0; i < rMesh3DData.Size(); i++)
    for (list<Surface3D>::const_iterator j = rMesh3DData[i].SurfaceBegin(); j != rMesh3DData[i].SurfaceEnd(); j++)
      ofsFile << VERTEX_CHARACTER << NORMAL_CHARACTER << ' ' << j->Vertex3DSurfaceNormal.fX << ' ' << j->Vertex3DSurfaceNormal.fY << ' ' << j->Vertex3DSurfaceNormal.fZ << endl;

  int iGlobalVertexCount = 0, iGlobalNormalCount = 0;

  for (int i = 0; i < rMesh3DData.Size(); i++)
    switch (rMesh3DData[i].ShapeState()) {
    case Plane :
      ofsFile << FACE_ELEMENT_CHARACTER << ' ';

      for (int j = 0; j < rMesh3DData[i].VertexCount(); j++) {
        ofsFile << iGlobalVertexCount + 1 << DELIMITER_CHARACTER << DELIMITER_CHARACTER << iGlobalNormalCount + 1;

        iGlobalVertexCount++;
      }

      iGlobalNormalCount++;

      break;
    case Solid :
      do {
        int iLocalVertexCount = 0, iLocalNormalCount = 0;

        for (list<Surface3D>::const_iterator j = rMesh3DData[i].SurfaceBegin(); j != rMesh3DData[i].SurfaceEnd(); j++) {
          ofsFile << FACE_ELEMENT_CHARACTER << ' ';

          const HalfEdge3D* pHalfEdge3DAnchor = j->pHalfEdge3DPointer;

          do {
            list<Vertex3D>::const_iterator iterPosition = rMesh3DData[i].VertexBegin();
            const Vertex3D* pVertex3DAnchor = pHalfEdge3DAnchor->pVertex3DPointer;

            iLocalVertexCount = 0;

            do {
              if (&*iterPosition == pVertex3DAnchor)
                break;

              iterPosition++;
              iLocalVertexCount++;
            } while (iterPosition != rMesh3DData[i].VertexEnd());

            ofsFile << iGlobalVertexCount + iLocalVertexCount + 1 << DELIMITER_CHARACTER << DELIMITER_CHARACTER << iGlobalNormalCount + iLocalNormalCount + 1;

            if ((pHalfEdge3DAnchor = pHalfEdge3DAnchor->pHalfEdge3DNext) == j->pHalfEdge3DPointer) {
              ofsFile << endl;

              break;
            } else
              ofsFile << ' ';
          } while (true);

          iLocalNormalCount++;
        }

        iGlobalVertexCount += rMesh3DData[i].VertexCount();
        iGlobalNormalCount += rMesh3DData[i].SurfaceCount();
      } while (false);

      break;
    case Undefined :
      cerr << "凸包が正しく生成されていません。" << endl;

      exit(EXIT_FAILURE);
    default : break;
    }

  ofsFile.close();
}

void ExportObjectFile3D(string strFileName, const vector<Mesh3D>& rvecPolygonSet)
{
  ofstream ofsFile;

  if (strFileName.empty()) {
    cout << "書き出すファイル名を入力してください：" << flush;

    getline(cin, strFileName);

    cout << endl;
  }

  ofsFile.open(strFileName + OBJECT_FILE_EXTENSION);

  if (!ofsFile.is_open()) {
    cerr << "ファイル「" << strFileName << OBJECT_FILE_EXTENSION << "」を開けませんでした。" << endl;

    exit(EXIT_FAILURE);
  }

  for (int i = 0; i < rvecPolygonSet.size(); i++)
    for (int j = 0; j < rvecPolygonSet[i].Size(); j++)
      for (list<Vertex3D>::const_iterator k = rvecPolygonSet[i][j].VertexBegin(); k != rvecPolygonSet[i][j].VertexEnd(); k++)
        ofsFile << VERTEX_CHARACTER << ' ' << k->fX << ' ' << k->fY << ' ' << k->fZ << endl;

  for (int i = 0; i < rvecPolygonSet.size(); i++)
    for (int j = 0; j < rvecPolygonSet[i].Size(); j++)
      for (list<Surface3D>::const_iterator k = rvecPolygonSet[i][j].SurfaceBegin(); k != rvecPolygonSet[i][j].SurfaceEnd(); k++)
        ofsFile << VERTEX_CHARACTER << NORMAL_CHARACTER << ' ' << k->Vertex3DSurfaceNormal.fX << ' ' << k->Vertex3DSurfaceNormal.fY << ' ' << k->Vertex3DSurfaceNormal.fZ << endl;

  int iGlobalVertexCount = 0, iGlobalNormalCount = 0;

  for (int i = 0; i < rvecPolygonSet.size(); i++)
    for (int j = 0; j < rvecPolygonSet[i].Size(); j++)
      switch (rvecPolygonSet[i][j].ShapeState()) {
      case Plane :
        ofsFile << FACE_ELEMENT_CHARACTER << ' ';

        for (int k = 0; k < rvecPolygonSet[i][j].VertexCount(); k++) {
          ofsFile << iGlobalVertexCount + 1 << DELIMITER_CHARACTER << DELIMITER_CHARACTER << iGlobalVertexCount + 1;

          iGlobalVertexCount++;

          if (rvecPolygonSet[i][j].VertexCount() - k - 1)
            ofsFile << ' ';
          else
            ofsFile << endl;
        }

        iGlobalNormalCount++;

        break;
      case Solid :
        do {
          int iLocalVertexCount = 0, iLocalNormalCount = 0;

          for (list<Surface3D>::const_iterator k = rvecPolygonSet[i][j].SurfaceBegin(); k != rvecPolygonSet[i][j].SurfaceEnd(); k++) {
            ofsFile << FACE_ELEMENT_CHARACTER << ' ';

            const HalfEdge3D* pHalfEdge3DAnchor = k->pHalfEdge3DPointer;

            do {
              list<Vertex3D>::const_iterator iterPosition = rvecPolygonSet[i][j].VertexBegin();
              const Vertex3D* pVertex3DAnchor = pHalfEdge3DAnchor->pVertex3DPointer;

              iLocalVertexCount = 0;

              do {
                if (&*iterPosition == pVertex3DAnchor)
                  break;

                iterPosition++;
                iLocalVertexCount++;
              } while (iterPosition != rvecPolygonSet[i][j].VertexEnd());

              ofsFile << iGlobalVertexCount + iLocalVertexCount + 1 << DELIMITER_CHARACTER << DELIMITER_CHARACTER << iGlobalNormalCount + iLocalNormalCount + 1;

              if ((pHalfEdge3DAnchor = pHalfEdge3DAnchor->pHalfEdge3DNext) == k->pHalfEdge3DPointer) {
                ofsFile << endl;

                break;
              } else
                ofsFile << ' ';
            } while (true);

            iLocalNormalCount++;
          }

          iGlobalVertexCount += rvecPolygonSet[i][j].VertexCount();
          iGlobalNormalCount += rvecPolygonSet[i][j].SurfaceCount();
        } while (false);

        break;
      case Undefined :
        cerr << "凸包が正しく生成されていません。" << endl;

        exit(EXIT_FAILURE);
      default : break;
      }

  ofsFile.close();
}

void ExportObjectFile3D(string strFileName, const PiecewisePolynomial3D& rPiecewisePolynomial3DData, const float (*const& rpfTransformationMatrix)[ELEMENT_COUNT_3D + 1], const int& riDivisionCountX, const int& riDivisionCountY, const int& riDivisionCountZ)
{
  ofstream ofsFile;

  if (strFileName.empty()) {
    cout << "書き出すファイル名を入力してください：" << flush;

    getline(cin, strFileName);
  }

  ofsFile.open(strFileName + OBJECT_FILE_EXTENSION);

  if (!ofsFile.is_open()) {
    cerr << "ファイル「" << strFileName << OBJECT_FILE_EXTENSION << "」を開けませんでした。" << endl;

    exit(EXIT_FAILURE);
  }

  MarchingCube(ofsFile, rPiecewisePolynomial3DData, rpfTransformationMatrix, riDivisionCountX, riDivisionCountY, riDivisionCountZ);

  ofsFile.close();
}

void ExportObjectFile3D(string strFileName, const vector<PiecewisePolynomial3D>& rvecImplicitSurface, const float (*const& rpfTransformationMatrix)[ELEMENT_COUNT_3D + 1][ELEMENT_COUNT_3D + 1], const int& riDivisionCountX, const int& riDivisionCountY, const int& riDivisionCountZ)
{
  ofstream ofsFile;

  if (strFileName.empty()) {
    cout << "書き出すファイル名を入力してください：" << flush;

    getline(cin, strFileName);
  }

  ofsFile.open(strFileName + OBJECT_FILE_EXTENSION);

  if (!ofsFile.is_open()) {
    cerr << "ファイル「" << strFileName << OBJECT_FILE_EXTENSION << "」を開けませんでした。" << endl;

    exit(EXIT_FAILURE);
  }

  MarchingCube(ofsFile, rvecImplicitSurface, rpfTransformationMatrix, riDivisionCountX, riDivisionCountY, riDivisionCountZ);

  ofsFile.close();
}

void ExportDataFile2D(string strFileName, const vector<Vertex2D>& rvecVertex, const vector<float>& rvecCoefficient, const float (*const& rpfTransformationMatrix)[ELEMENT_COUNT_2D + 1], const int& riDivisionCountX, const int& riDivisionCountY, const Vertex2D& rVertex2DMinimumCorner, const Vertex2D& rVertex2DMaximumCorner, const float& rfPadding)
{
  ofstream ofsFile;

  if (strFileName.empty()) {
    cout << "書き出すファイル名を入力してください：" << flush;

    getline(cin, strFileName);

    cout << endl;
  }

  ofsFile.open(strFileName + DATA_FILE_EXTENSION);

  if (!ofsFile.is_open()) {
    cerr << "ファイル「" << strFileName << DATA_FILE_EXTENSION << "」を開けませんでした。" << endl;

    exit(EXIT_FAILURE);
  }

  if (rpfTransformationMatrix == nullptr) {
    if (rVertex2DMinimumCorner.fX < rVertex2DMaximumCorner.fX && rVertex2DMinimumCorner.fY < rVertex2DMaximumCorner.fY) {
      float fStep[ELEMENT_COUNT_2D] = {(rVertex2DMaximumCorner.fX - rVertex2DMinimumCorner.fX) / static_cast<float>(riDivisionCountX), (rVertex2DMaximumCorner.fY - rVertex2DMinimumCorner.fY) / static_cast<float>(riDivisionCountY)};

      for (int i = 0; i < riDivisionCountX; i++) {
        for (int j = 0; j < riDivisionCountY; j++) {
          float fSamplingPoint[ELEMENT_COUNT_2D] = {rVertex2DMinimumCorner.fX + fStep[INDICATOR_X] * (static_cast<float>(i) + 0.5F), rVertex2DMinimumCorner.fY + fStep[INDICATOR_Y] * (static_cast<float>(j) + 0.5F)};

          ofsFile << fSamplingPoint[INDICATOR_X] << ' ' << fSamplingPoint[INDICATOR_Y] << ' ' << FunctionValue2D(fSamplingPoint, rvecVertex, rvecCoefficient) << endl;
        }

        if (riDivisionCountX - i - 1)
          ofsFile << endl;
      }
    } else {
      float Vertex2D::*pfVertex2DMember[ELEMENT_COUNT_2D] = {&Vertex2D::fX, &Vertex2D::fY};
      ExtremePoint ExtremePointMatrix[ELEMENT_COUNT_2D][EXTREME_ELEMENT_COUNT];

      for (int i = 0; i < ELEMENT_COUNT_2D; i++)
        for (int j = 0; j < EXTREME_ELEMENT_COUNT; j++) {
          ExtremePointMatrix[i][j].uIndex = 0U;
          ExtremePointMatrix[i][j].fDistance = *rvecVertex.begin().*pfVertex2DMember[i];
        }

      for (int i = 0; i < rvecVertex.size(); i++)
        for (int j = 0; j < ELEMENT_COUNT_2D; j++)
          if (rvecVertex[i].*pfVertex2DMember[j] < ExtremePointMatrix[j][INDICATOR_MINIMUM].fDistance) {
            ExtremePointMatrix[j][INDICATOR_MINIMUM].uIndex = static_cast<unsigned>(i);
            ExtremePointMatrix[j][INDICATOR_MINIMUM].fDistance = rvecVertex[i].*pfVertex2DMember[j];
          } else if (rvecVertex[i].*pfVertex2DMember[j] > ExtremePointMatrix[j][INDICATOR_MAXIMUM].fDistance) {
            ExtremePointMatrix[j][INDICATOR_MAXIMUM].uIndex = static_cast<unsigned>(i);
            ExtremePointMatrix[j][INDICATOR_MAXIMUM].fDistance = rvecVertex[i].*pfVertex2DMember[j];
          }

      float fDifference = 0.0F;

      if (ExtremePointMatrix[INDICATOR_X][INDICATOR_MAXIMUM].fDistance - ExtremePointMatrix[INDICATOR_X][INDICATOR_MINIMUM].fDistance >= ExtremePointMatrix[INDICATOR_Y][INDICATOR_MAXIMUM].fDistance - ExtremePointMatrix[INDICATOR_Y][INDICATOR_MINIMUM].fDistance)
        fDifference = ExtremePointMatrix[INDICATOR_X][INDICATOR_MAXIMUM].fDistance - ExtremePointMatrix[INDICATOR_X][INDICATOR_MINIMUM].fDistance;
      else
        fDifference = ExtremePointMatrix[INDICATOR_Y][INDICATOR_MAXIMUM].fDistance - ExtremePointMatrix[INDICATOR_Y][INDICATOR_MINIMUM].fDistance;

      float fPaddedCorner[ELEMENT_COUNT_2D][EXTREME_ELEMENT_COUNT];

      for (int i = 0; i < ELEMENT_COUNT_2D; i++) {
        fPaddedCorner[i][INDICATOR_MINIMUM] = (ExtremePointMatrix[i][INDICATOR_MINIMUM].fDistance + ExtremePointMatrix[i][INDICATOR_MAXIMUM].fDistance) / 2.0F - (rfPadding + 0.5F) * fDifference;
        fPaddedCorner[i][INDICATOR_MAXIMUM] = (ExtremePointMatrix[i][INDICATOR_MINIMUM].fDistance + ExtremePointMatrix[i][INDICATOR_MAXIMUM].fDistance) / 2.0F + (rfPadding + 0.5F) * fDifference;
      }

      float fStep[ELEMENT_COUNT_2D] = {(fPaddedCorner[INDICATOR_X][INDICATOR_MAXIMUM] - fPaddedCorner[INDICATOR_X][INDICATOR_MINIMUM]) / static_cast<float>(riDivisionCountX), (fPaddedCorner[INDICATOR_Y][INDICATOR_MAXIMUM] - fPaddedCorner[INDICATOR_Y][INDICATOR_MINIMUM]) / static_cast<float>(riDivisionCountY)};

      for (int i = 0; i < riDivisionCountX; i++) {
        for (int j = 0; j < riDivisionCountY; j++) {
          float fSamplingPoint[ELEMENT_COUNT_2D] = {fPaddedCorner[INDICATOR_X][INDICATOR_MINIMUM] + fStep[INDICATOR_X] * (static_cast<float>(i) + 0.5F), fPaddedCorner[INDICATOR_Y][INDICATOR_MINIMUM] + fStep[INDICATOR_Y] * (static_cast<float>(j) + 0.5F)};

          ofsFile << fSamplingPoint[INDICATOR_X] << ' ' << fSamplingPoint[INDICATOR_Y] << ' ' << FunctionValue2D(fSamplingPoint, rvecVertex, rvecCoefficient) << endl;
        }

        if (riDivisionCountX - i - 1)
          ofsFile << endl;
      }
    }
  } else {
    float fInversedTransformationMatrix[ELEMENT_COUNT_2D + 1][ELEMENT_COUNT_2D + 1];

    for (int i = 0; i <= ELEMENT_COUNT_2D; i++)
      for (int j = 0; j <= ELEMENT_COUNT_2D; j++)
        fInversedTransformationMatrix[i][j] = rpfTransformationMatrix[i][j];

    InverseTransformation<ELEMENT_COUNT_2D + 1>(fInversedTransformationMatrix);

    if (rVertex2DMinimumCorner.fX < rVertex2DMaximumCorner.fX && rVertex2DMinimumCorner.fY < rVertex2DMaximumCorner.fY) {
      float fStep[ELEMENT_COUNT_2D] = {(rVertex2DMaximumCorner.fX - rVertex2DMinimumCorner.fX) / static_cast<float>(riDivisionCountX), (rVertex2DMaximumCorner.fY - rVertex2DMinimumCorner.fY) / static_cast<float>(riDivisionCountY)};

      for (int i = 0; i < riDivisionCountX; i++) {
        for (int j = 0; j < riDivisionCountY; j++) {
          float fSamplingPoint[ELEMENT_COUNT_2D] = {rVertex2DMinimumCorner.fX + fStep[INDICATOR_X] * (static_cast<float>(i) + 0.5F), rVertex2DMinimumCorner.fY + fStep[INDICATOR_Y] * (static_cast<float>(j) + 0.5F)}, fTransformedSamplingPoint[ELEMENT_COUNT_2D] = {0.0F, 0.0F};

          TransformCoordinate<float, ELEMENT_COUNT_2D + 1>(fInversedTransformationMatrix, fSamplingPoint, fTransformedSamplingPoint);

          ofsFile << fSamplingPoint[INDICATOR_X] << ' ' << fSamplingPoint[INDICATOR_Y] << ' ' << FunctionValue2D(fTransformedSamplingPoint, rvecVertex, rvecCoefficient) << endl;
        }

        if (riDivisionCountX - i - 1)
          ofsFile << endl;
      }
    } else {
      float Vertex2D::*pfVertex2DMember[ELEMENT_COUNT_2D] = {&Vertex2D::fX, &Vertex2D::fY};
      ExtremePoint ExtremePointMatrix[ELEMENT_COUNT_2D][EXTREME_ELEMENT_COUNT];

      for (int i = 0; i < ELEMENT_COUNT_2D; i++)
        for (int j = 0; j < EXTREME_ELEMENT_COUNT; j++) {
          ExtremePointMatrix[i][j].uIndex = 0U;
          ExtremePointMatrix[i][j].fDistance = *rvecVertex.begin().*pfVertex2DMember[i];
        }

      for (int i = 0; i < rvecVertex.size(); i++)
        for (int j = 0; j < ELEMENT_COUNT_2D; j++)
          if (rvecVertex[i].*pfVertex2DMember[j] < ExtremePointMatrix[j][INDICATOR_MINIMUM].fDistance) {
            ExtremePointMatrix[j][INDICATOR_MINIMUM].uIndex = static_cast<unsigned>(i);
            ExtremePointMatrix[j][INDICATOR_MINIMUM].fDistance = rvecVertex[i].*pfVertex2DMember[j];
          } else if (rvecVertex[i].*pfVertex2DMember[j] > ExtremePointMatrix[j][INDICATOR_MAXIMUM].fDistance) {
            ExtremePointMatrix[j][INDICATOR_MAXIMUM].uIndex = static_cast<unsigned>(i);
            ExtremePointMatrix[j][INDICATOR_MAXIMUM].fDistance = rvecVertex[i].*pfVertex2DMember[j];
          }

      float fDifference = 0.0F;

      if (ExtremePointMatrix[INDICATOR_X][INDICATOR_MAXIMUM].fDistance - ExtremePointMatrix[INDICATOR_X][INDICATOR_MINIMUM].fDistance >= ExtremePointMatrix[INDICATOR_Y][INDICATOR_MAXIMUM].fDistance - ExtremePointMatrix[INDICATOR_Y][INDICATOR_MINIMUM].fDistance)
        fDifference = ExtremePointMatrix[INDICATOR_X][INDICATOR_MAXIMUM].fDistance - ExtremePointMatrix[INDICATOR_X][INDICATOR_MINIMUM].fDistance;
      else
        fDifference = ExtremePointMatrix[INDICATOR_Y][INDICATOR_MAXIMUM].fDistance - ExtremePointMatrix[INDICATOR_Y][INDICATOR_MINIMUM].fDistance;

      float fPaddedCorner[ELEMENT_COUNT_2D][EXTREME_ELEMENT_COUNT];

      for (int i = 0; i < ELEMENT_COUNT_2D; i++) {
        fPaddedCorner[i][INDICATOR_MINIMUM] = (ExtremePointMatrix[i][INDICATOR_MINIMUM].fDistance + ExtremePointMatrix[i][INDICATOR_MAXIMUM].fDistance) / 2.0F - (rfPadding + 0.5F) * fDifference;
        fPaddedCorner[i][INDICATOR_MAXIMUM] = (ExtremePointMatrix[i][INDICATOR_MINIMUM].fDistance + ExtremePointMatrix[i][INDICATOR_MAXIMUM].fDistance) / 2.0F + (rfPadding + 0.5F) * fDifference;
      }

      float fTransformedCorner[ELEMENT_COUNT_2D][EXTREME_ELEMENT_COUNT];

      for (int i = 0; i < ELEMENT_COUNT_2D; i++) {
        fTransformedCorner[i][INDICATOR_MINIMUM] = fTransformedCorner[i][INDICATOR_MAXIMUM] = 0.0F;

        for (int j = 0; j < ELEMENT_COUNT_2D; j++)
          if (rpfTransformationMatrix[i][j] * fPaddedCorner[j][INDICATOR_MINIMUM] < rpfTransformationMatrix[i][j] * fPaddedCorner[j][INDICATOR_MAXIMUM]) {
            fTransformedCorner[i][INDICATOR_MINIMUM] += rpfTransformationMatrix[i][j] * fPaddedCorner[j][INDICATOR_MINIMUM];
            fTransformedCorner[i][INDICATOR_MAXIMUM] += rpfTransformationMatrix[i][j] * fPaddedCorner[j][INDICATOR_MAXIMUM];
          } else {
            fTransformedCorner[i][INDICATOR_MINIMUM] += rpfTransformationMatrix[i][j] * fPaddedCorner[j][INDICATOR_MAXIMUM];
            fTransformedCorner[i][INDICATOR_MAXIMUM] += rpfTransformationMatrix[i][j] * fPaddedCorner[j][INDICATOR_MINIMUM];
          }

        fTransformedCorner[i][INDICATOR_MINIMUM] += rpfTransformationMatrix[i][ELEMENT_COUNT_2D];
        fTransformedCorner[i][INDICATOR_MAXIMUM] += rpfTransformationMatrix[i][ELEMENT_COUNT_2D];
      }

      float fStep[ELEMENT_COUNT_2D] = {(fTransformedCorner[INDICATOR_X][INDICATOR_MAXIMUM] - fTransformedCorner[INDICATOR_X][INDICATOR_MINIMUM]) / static_cast<float>(riDivisionCountX), (fTransformedCorner[INDICATOR_Y][INDICATOR_MAXIMUM] - fTransformedCorner[INDICATOR_Y][INDICATOR_MINIMUM]) / static_cast<float>(riDivisionCountY)};

      for (int i = 0; i < riDivisionCountX; i++) {
        for (int j = 0; j < riDivisionCountY; j++) {
          float fSamplingPoint[ELEMENT_COUNT_2D] = {fTransformedCorner[INDICATOR_X][INDICATOR_MINIMUM] + fStep[INDICATOR_X] * (static_cast<float>(i) + 0.5F), fTransformedCorner[INDICATOR_Y][INDICATOR_MINIMUM] + fStep[INDICATOR_Y] * (static_cast<float>(j) + 0.5F)}, fTransformedSamplingPoint[ELEMENT_COUNT_2D] = {0.0F, 0.0F};

          TransformCoordinate<float, ELEMENT_COUNT_2D + 1>(fInversedTransformationMatrix, fSamplingPoint, fTransformedSamplingPoint);

          ofsFile << fSamplingPoint[INDICATOR_X] << ' ' << fSamplingPoint[INDICATOR_Y] << ' ' << FunctionValue2D(fTransformedSamplingPoint, rvecVertex, rvecCoefficient) << endl;
        }

        if (riDivisionCountX - i - 1)
          ofsFile << endl;
      }
    }
  }

  ofsFile.close();
}

void ExportDataFile2D(string strFileName, const PiecewisePolynomial2D& rPiecewisePolynomial2DData, const float (*const& rpfTransformationMatrix)[ELEMENT_COUNT_2D + 1], const int& riDivisionCountX, const int& riDivisionCountY, const Vertex2D& rVertex2DMinimumCorner, const Vertex2D& rVertex2DMaximumCorner)
{
  ofstream ofsFile;

  if (strFileName.empty()) {
    cout << "書き出すファイル名を入力してください：" << flush;

    getline(cin, strFileName);

    cout << endl;
  }

  ofsFile.open(strFileName + DATA_FILE_EXTENSION);

  if (!ofsFile.is_open()) {
    cerr << "ファイル「" << strFileName << DATA_FILE_EXTENSION << "」を開けませんでした。" << endl;

    exit(EXIT_FAILURE);
  }

  if (rpfTransformationMatrix == nullptr) {
    if (rVertex2DMinimumCorner.fX < rVertex2DMaximumCorner.fX && rVertex2DMinimumCorner.fY < rVertex2DMaximumCorner.fY) {
      float fStep[ELEMENT_COUNT_2D] = {(rVertex2DMaximumCorner.fX - rVertex2DMinimumCorner.fX) / static_cast<float>(riDivisionCountX), (rVertex2DMaximumCorner.fY - rVertex2DMinimumCorner.fY) / static_cast<float>(riDivisionCountY)};

      for (int i = 0; i < riDivisionCountX; i++) {
        for (int j = 0; j < riDivisionCountY; j++) {
          Vertex2D Vertex2DSamplingPoint{rVertex2DMinimumCorner.fX + fStep[INDICATOR_X] * (static_cast<float>(i) + 0.5F), rVertex2DMinimumCorner.fY + fStep[INDICATOR_Y] * (static_cast<float>(j) + 0.5F)};

          ofsFile << Vertex2DSamplingPoint.fX << ' ' << Vertex2DSamplingPoint.fY << ' ' << rPiecewisePolynomial2DData.At(Vertex2DSamplingPoint) << endl;
        }

        if (riDivisionCountX - i - 1)
          ofsFile << endl;
      }
    } else {
      float fStep[ELEMENT_COUNT_2D] = {(rPiecewisePolynomial2DData.MaximumCorner().fX - rPiecewisePolynomial2DData.MinimumCorner().fX) / static_cast<float>(riDivisionCountX), (rPiecewisePolynomial2DData.MaximumCorner().fY - rPiecewisePolynomial2DData.MinimumCorner().fY) / static_cast<float>(riDivisionCountY)};

      for (int i = 0; i < riDivisionCountX; i++) {
        for (int j = 0; j < riDivisionCountY; j++) {
          Vertex2D Vertex2DSamplingPoint{rPiecewisePolynomial2DData.MinimumCorner().fX + fStep[INDICATOR_X] * (static_cast<float>(i) + 0.5F), rPiecewisePolynomial2DData.MinimumCorner().fY + fStep[INDICATOR_Y] * (static_cast<float>(j) + 0.5F)};

          ofsFile << Vertex2DSamplingPoint.fX << ' ' << Vertex2DSamplingPoint.fY << ' ' << rPiecewisePolynomial2DData.At(Vertex2DSamplingPoint) << endl;
        }

        if (riDivisionCountX - i - 1)
          ofsFile << endl;
      }
    }
  } else {
    float fInversedTransformationMatrix[ELEMENT_COUNT_2D + 1][ELEMENT_COUNT_2D + 1];

    for (int i = 0; i <= ELEMENT_COUNT_2D; i++)
      for (int j = 0; j <= ELEMENT_COUNT_2D; j++)
        fInversedTransformationMatrix[i][j] = rpfTransformationMatrix[i][j];

    InverseTransformation<ELEMENT_COUNT_2D + 1>(fInversedTransformationMatrix);

    if (rVertex2DMinimumCorner.fX < rVertex2DMaximumCorner.fX && rVertex2DMinimumCorner.fY < rVertex2DMaximumCorner.fY) {
      float fStep[ELEMENT_COUNT_2D] = {(rVertex2DMaximumCorner.fX - rVertex2DMinimumCorner.fX) / static_cast<float>(riDivisionCountX), (rVertex2DMaximumCorner.fY - rVertex2DMinimumCorner.fY) / static_cast<float>(riDivisionCountY)};

      for (int i = 0; i < riDivisionCountX; i++) {
        for (int j = 0; j < riDivisionCountY; j++) {
          Vertex2D Vertex2DSamplingPoint{rVertex2DMinimumCorner.fX + fStep[INDICATOR_X] * (static_cast<float>(i) + 0.5F), rVertex2DMinimumCorner.fY + fStep[INDICATOR_Y] * (static_cast<float>(j) + 0.5F)};

          ofsFile << Vertex2DSamplingPoint.fX << ' ' << Vertex2DSamplingPoint.fY << ' ' << rPiecewisePolynomial2DData.At(TransformCoordinate(fInversedTransformationMatrix, Vertex2DSamplingPoint)) << endl;
        }

        if (riDivisionCountX - i - 1)
          ofsFile << endl;
      }
    } else {
      float Vertex2D::*pfVertex2DMember[ELEMENT_COUNT_2D] = {&Vertex2D::fX, &Vertex2D::fY};
      vector<Vertex2D> vecOriginalCorner{rPiecewisePolynomial2DData.MinimumCorner(), rPiecewisePolynomial2DData.MaximumCorner()}, vecTransformedCorner(ELEMENT_COUNT_2D);

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

      float fStep[ELEMENT_COUNT_2D] = {(vecTransformedCorner[INDICATOR_MAXIMUM].fX - vecTransformedCorner[INDICATOR_MINIMUM].fX) / static_cast<float>(riDivisionCountX), (vecTransformedCorner[INDICATOR_MAXIMUM].fY - vecTransformedCorner[INDICATOR_MINIMUM].fY) / static_cast<float>(riDivisionCountY)};

      for (int i = 0; i < riDivisionCountX; i++) {
        for (int j = 0; j < riDivisionCountY; j++) {
          Vertex2D Vertex2DSamplingPoint{vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * (static_cast<float>(i) + 0.5F), vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * (static_cast<float>(j) + 0.5F)};

          ofsFile << Vertex2DSamplingPoint.fX << ' ' << Vertex2DSamplingPoint.fY << ' ' << rPiecewisePolynomial2DData.At(TransformCoordinate(fInversedTransformationMatrix, Vertex2DSamplingPoint)) << endl;
        }

        if (riDivisionCountX - i - 1)
          ofsFile << endl;
      }
    }
  }

  ofsFile.close();
}

void ExportDataFile3D(string strFileName, const vector<Vertex3D>& rvecVertex, const vector<float>& rvecCoefficient, const float (*const& rpfTransformationMatrix)[ELEMENT_COUNT_3D + 1], const int& riDivisionCountX, const int& riDivisionCountY, const int& riDivisionCountZ, const Vertex3D& rVertex3DMinimumCorner, const Vertex3D& rVertex3DMaximumCorner, const float& rfPadding)
{
  ofstream ofsFile;

  if (strFileName.empty()) {
    cout << "書き出すファイル名を入力してください：" << flush;

    getline(cin, strFileName);

    cout << endl;
  }

  ofsFile.open(strFileName + DATA_FILE_EXTENSION);

  if (!ofsFile.is_open()) {
    cerr << "ファイル「" << strFileName << DATA_FILE_EXTENSION << "」を開けませんでした。" << endl;

    exit(EXIT_FAILURE);
  }

  if (rpfTransformationMatrix == nullptr) {
    if (rVertex3DMinimumCorner.fX < rVertex3DMaximumCorner.fX && rVertex3DMinimumCorner.fY < rVertex3DMaximumCorner.fY && rVertex3DMinimumCorner.fZ < rVertex3DMaximumCorner.fZ) {
      float fStep[ELEMENT_COUNT_3D] = {(rVertex3DMaximumCorner.fX - rVertex3DMinimumCorner.fX) / static_cast<float>(riDivisionCountX), (rVertex3DMaximumCorner.fY - rVertex3DMinimumCorner.fY) / static_cast<float>(riDivisionCountY), (rVertex3DMaximumCorner.fZ - rVertex3DMinimumCorner.fZ) / static_cast<float>(riDivisionCountZ)};

      for (int i = 0; i < riDivisionCountX; i++) {
        for (int j = 0; j < riDivisionCountY; j++) {
          for (int k = 0; k < riDivisionCountZ; k++) {
            float fSamplingPoint[ELEMENT_COUNT_3D] = {rVertex3DMinimumCorner.fX + fStep[INDICATOR_X] * (static_cast<float>(i) + 0.5F), rVertex3DMinimumCorner.fY + fStep[INDICATOR_Y] * (static_cast<float>(j) + 0.5F), rVertex3DMinimumCorner.fZ + fStep[INDICATOR_Z] * (static_cast<float>(k) + 0.5F)};

            ofsFile << fSamplingPoint[INDICATOR_X] << ' ' << fSamplingPoint[INDICATOR_Y] << ' ' << fSamplingPoint[INDICATOR_Z] << ' ' << FunctionValue3D(fSamplingPoint, rvecVertex, rvecCoefficient) << endl;
          }

          if (riDivisionCountY - j - 1)
            ofsFile << endl;
        }

        if (riDivisionCountX - i - 1)
          ofsFile << endl;
      }
    } else {
      float Vertex3D::*pfVertex3DMember[ELEMENT_COUNT_3D] = {&Vertex3D::fX, &Vertex3D::fY, &Vertex3D::fZ};
      ExtremePoint ExtremePointMatrix[ELEMENT_COUNT_3D][EXTREME_ELEMENT_COUNT];

      for (int i = 0; i < ELEMENT_COUNT_3D; i++)
        for (int j = 0; j < EXTREME_ELEMENT_COUNT; j++) {
          ExtremePointMatrix[i][j].uIndex = 0U;
          ExtremePointMatrix[i][j].fDistance = *rvecVertex.begin().*pfVertex3DMember[i];
        }

      for (int i = 0; i < rvecVertex.size(); i++)
        for (int j = 0; j < ELEMENT_COUNT_3D; j++) {
          if (rvecVertex[i].*pfVertex3DMember[j] < ExtremePointMatrix[j][INDICATOR_MINIMUM].fDistance) {
            ExtremePointMatrix[j][INDICATOR_MINIMUM].uIndex = static_cast<unsigned>(i);
            ExtremePointMatrix[j][INDICATOR_MINIMUM].fDistance = rvecVertex[i].*pfVertex3DMember[j];
          } else if (rvecVertex[i].*pfVertex3DMember[j] > ExtremePointMatrix[j][INDICATOR_MAXIMUM].fDistance) {
            ExtremePointMatrix[j][INDICATOR_MAXIMUM].uIndex = static_cast<unsigned>(i);
            ExtremePointMatrix[j][INDICATOR_MAXIMUM].fDistance = rvecVertex[i].*pfVertex3DMember[j];
          }
        }

      float fDifference = 0.0F;

      if (ExtremePointMatrix[INDICATOR_X][INDICATOR_MAXIMUM].fDistance - ExtremePointMatrix[INDICATOR_X][INDICATOR_MINIMUM].fDistance >= ExtremePointMatrix[INDICATOR_Y][INDICATOR_MAXIMUM].fDistance - ExtremePointMatrix[INDICATOR_Y][INDICATOR_MINIMUM].fDistance)
        if (ExtremePointMatrix[INDICATOR_X][INDICATOR_MAXIMUM].fDistance - ExtremePointMatrix[INDICATOR_X][INDICATOR_MINIMUM].fDistance >= ExtremePointMatrix[INDICATOR_Z][INDICATOR_MAXIMUM].fDistance - ExtremePointMatrix[INDICATOR_Z][INDICATOR_MINIMUM].fDistance)
          fDifference = ExtremePointMatrix[INDICATOR_X][INDICATOR_MAXIMUM].fDistance - ExtremePointMatrix[INDICATOR_X][INDICATOR_MINIMUM].fDistance;
        else
          fDifference = ExtremePointMatrix[INDICATOR_Z][INDICATOR_MAXIMUM].fDistance - ExtremePointMatrix[INDICATOR_Z][INDICATOR_MINIMUM].fDistance;
      else if (ExtremePointMatrix[INDICATOR_Y][INDICATOR_MAXIMUM].fDistance - ExtremePointMatrix[INDICATOR_Y][INDICATOR_MINIMUM].fDistance >= ExtremePointMatrix[INDICATOR_Z][INDICATOR_MAXIMUM].fDistance - ExtremePointMatrix[INDICATOR_Z][INDICATOR_MINIMUM].fDistance)
        fDifference = ExtremePointMatrix[INDICATOR_Y][INDICATOR_MAXIMUM].fDistance - ExtremePointMatrix[INDICATOR_Y][INDICATOR_MINIMUM].fDistance;
      else
        fDifference = ExtremePointMatrix[INDICATOR_Z][INDICATOR_MAXIMUM].fDistance - ExtremePointMatrix[INDICATOR_Z][INDICATOR_MINIMUM].fDistance;

      float fPaddedCorner[ELEMENT_COUNT_3D][EXTREME_ELEMENT_COUNT];

      for (int i = 0; i < ELEMENT_COUNT_3D; i++) {
        fPaddedCorner[i][INDICATOR_MINIMUM] = (ExtremePointMatrix[i][INDICATOR_MINIMUM].fDistance + ExtremePointMatrix[i][INDICATOR_MAXIMUM].fDistance) / 2.0F - (rfPadding + 0.5F) * fDifference;
        fPaddedCorner[i][INDICATOR_MAXIMUM] = (ExtremePointMatrix[i][INDICATOR_MINIMUM].fDistance + ExtremePointMatrix[i][INDICATOR_MAXIMUM].fDistance) / 2.0F + (rfPadding + 0.5F) * fDifference;
      }

      float fStep[ELEMENT_COUNT_3D] = {(fPaddedCorner[INDICATOR_X][INDICATOR_MAXIMUM] - fPaddedCorner[INDICATOR_X][INDICATOR_MINIMUM]) / static_cast<float>(riDivisionCountX), (fPaddedCorner[INDICATOR_Y][INDICATOR_MAXIMUM] - fPaddedCorner[INDICATOR_Y][INDICATOR_MINIMUM]) / static_cast<float>(riDivisionCountY), (fPaddedCorner[INDICATOR_Z][INDICATOR_MAXIMUM] - fPaddedCorner[INDICATOR_Z][INDICATOR_MINIMUM]) / static_cast<float>(riDivisionCountZ)};

      for (int i = 0; i < riDivisionCountX; i++) {
        for (int j = 0; j < riDivisionCountY; j++) {
          for (int k = 0; k < riDivisionCountZ; k++) {
            float fSamplingPoint[ELEMENT_COUNT_3D] = {fPaddedCorner[INDICATOR_X][INDICATOR_MINIMUM] + fStep[INDICATOR_X] * (static_cast<float>(i) + 0.5F), fPaddedCorner[INDICATOR_Y][INDICATOR_MINIMUM] + fStep[INDICATOR_Y] * (static_cast<float>(j) + 0.5F), fPaddedCorner[INDICATOR_Z][INDICATOR_MINIMUM] + fStep[INDICATOR_Z] * (static_cast<float>(k) + 0.5F)};

            ofsFile << fSamplingPoint[INDICATOR_X] << ' ' << fSamplingPoint[INDICATOR_Y] << ' ' << fSamplingPoint[INDICATOR_Z] << ' ' << FunctionValue3D(fSamplingPoint, rvecVertex, rvecCoefficient) << endl;
          }

          if (riDivisionCountY - j - 1)
            ofsFile << endl;
        }

        if (riDivisionCountX - i - 1)
          ofsFile << endl;
      }
    }
  } else {
    float fInversedTransformationMatrix[ELEMENT_COUNT_3D + 1][ELEMENT_COUNT_3D + 1];

    for (int i = 0; i <= ELEMENT_COUNT_3D; i++)
      for (int j = 0; j <= ELEMENT_COUNT_3D; j++)
        fInversedTransformationMatrix[i][j] = rpfTransformationMatrix[i][j];

    InverseTransformation<ELEMENT_COUNT_3D + 1>(fInversedTransformationMatrix);

    if (rVertex3DMinimumCorner.fX < rVertex3DMaximumCorner.fX && rVertex3DMinimumCorner.fY < rVertex3DMaximumCorner.fY && rVertex3DMinimumCorner.fZ < rVertex3DMaximumCorner.fZ) {
      float fStep[ELEMENT_COUNT_3D] = {(rVertex3DMaximumCorner.fX - rVertex3DMinimumCorner.fX) / static_cast<float>(riDivisionCountX), (rVertex3DMaximumCorner.fY - rVertex3DMinimumCorner.fY) / static_cast<float>(riDivisionCountY), (rVertex3DMaximumCorner.fZ - rVertex3DMinimumCorner.fZ) / static_cast<float>(riDivisionCountZ)};

      for (int i = 0; i < riDivisionCountX; i++) {
        for (int j = 0; j < riDivisionCountY; j++) {
          for (int k = 0; k < riDivisionCountZ; k++) {
            float fSamplingPoint[ELEMENT_COUNT_3D] = {rVertex3DMinimumCorner.fX + fStep[INDICATOR_X] * (static_cast<float>(i) + 0.5F), rVertex3DMinimumCorner.fY + fStep[INDICATOR_Y] * (static_cast<float>(j) + 0.5F), rVertex3DMinimumCorner.fZ + fStep[INDICATOR_Z] * (static_cast<float>(k) + 0.5F)}, fTransformedSamplingPoint[ELEMENT_COUNT_3D] = {0.0F, 0.0F, 0.0F};

            TransformCoordinate<float, ELEMENT_COUNT_3D + 1>(fInversedTransformationMatrix, fSamplingPoint, fTransformedSamplingPoint);

            ofsFile << fSamplingPoint[INDICATOR_X] << ' ' << fSamplingPoint[INDICATOR_Y] << ' ' << fSamplingPoint[INDICATOR_Z] << ' ' << FunctionValue3D(fTransformedSamplingPoint, rvecVertex, rvecCoefficient) << endl;
          }

          if (riDivisionCountY - j - 1)
            ofsFile << endl;
        }

        if (riDivisionCountX - i - 1)
          ofsFile << endl;
      }
    } else {
      float Vertex3D::*pfVertex3DMember[ELEMENT_COUNT_3D] = {&Vertex3D::fX, &Vertex3D::fY, &Vertex3D::fZ};
      ExtremePoint ExtremePointMatrix[ELEMENT_COUNT_3D][EXTREME_ELEMENT_COUNT];

      for (int i = 0; i < ELEMENT_COUNT_3D; i++)
        for (int j = 0; j < EXTREME_ELEMENT_COUNT; j++) {
          ExtremePointMatrix[i][j].uIndex = 0U;
          ExtremePointMatrix[i][j].fDistance = *rvecVertex.begin().*pfVertex3DMember[i];
        }

      for (int i = 0; i < rvecVertex.size(); i++)
        for (int j = 0; j < ELEMENT_COUNT_3D; j++) {
          if (rvecVertex[i].*pfVertex3DMember[j] < ExtremePointMatrix[j][INDICATOR_MINIMUM].fDistance) {
            ExtremePointMatrix[j][INDICATOR_MINIMUM].uIndex = static_cast<unsigned>(i);
            ExtremePointMatrix[j][INDICATOR_MINIMUM].fDistance = rvecVertex[i].*pfVertex3DMember[j];
          } else if (rvecVertex[i].*pfVertex3DMember[j] > ExtremePointMatrix[j][INDICATOR_MAXIMUM].fDistance) {
            ExtremePointMatrix[j][INDICATOR_MAXIMUM].uIndex = static_cast<unsigned>(i);
            ExtremePointMatrix[j][INDICATOR_MAXIMUM].fDistance = rvecVertex[i].*pfVertex3DMember[j];
          }
        }

      float fDifference = 0.0F;

      if (ExtremePointMatrix[INDICATOR_X][INDICATOR_MAXIMUM].fDistance - ExtremePointMatrix[INDICATOR_X][INDICATOR_MINIMUM].fDistance >= ExtremePointMatrix[INDICATOR_Y][INDICATOR_MAXIMUM].fDistance - ExtremePointMatrix[INDICATOR_Y][INDICATOR_MINIMUM].fDistance)
        if (ExtremePointMatrix[INDICATOR_X][INDICATOR_MAXIMUM].fDistance - ExtremePointMatrix[INDICATOR_X][INDICATOR_MINIMUM].fDistance >= ExtremePointMatrix[INDICATOR_Z][INDICATOR_MAXIMUM].fDistance - ExtremePointMatrix[INDICATOR_Z][INDICATOR_MINIMUM].fDistance)
          fDifference = ExtremePointMatrix[INDICATOR_X][INDICATOR_MAXIMUM].fDistance - ExtremePointMatrix[INDICATOR_X][INDICATOR_MINIMUM].fDistance;
        else
          fDifference = ExtremePointMatrix[INDICATOR_Z][INDICATOR_MAXIMUM].fDistance - ExtremePointMatrix[INDICATOR_Z][INDICATOR_MINIMUM].fDistance;
      else if (ExtremePointMatrix[INDICATOR_Y][INDICATOR_MAXIMUM].fDistance - ExtremePointMatrix[INDICATOR_Y][INDICATOR_MINIMUM].fDistance >= ExtremePointMatrix[INDICATOR_Z][INDICATOR_MAXIMUM].fDistance - ExtremePointMatrix[INDICATOR_Z][INDICATOR_MINIMUM].fDistance)
        fDifference = ExtremePointMatrix[INDICATOR_Y][INDICATOR_MAXIMUM].fDistance - ExtremePointMatrix[INDICATOR_Y][INDICATOR_MINIMUM].fDistance;
      else
        fDifference = ExtremePointMatrix[INDICATOR_Z][INDICATOR_MAXIMUM].fDistance - ExtremePointMatrix[INDICATOR_Z][INDICATOR_MINIMUM].fDistance;

      float fPaddedCorner[ELEMENT_COUNT_3D][EXTREME_ELEMENT_COUNT];

      for (int i = 0; i < ELEMENT_COUNT_3D; i++) {
        fPaddedCorner[i][INDICATOR_MINIMUM] = (ExtremePointMatrix[i][INDICATOR_MINIMUM].fDistance + ExtremePointMatrix[i][INDICATOR_MAXIMUM].fDistance) / 2.0F - (rfPadding + 0.5F) * fDifference;
        fPaddedCorner[i][INDICATOR_MAXIMUM] = (ExtremePointMatrix[i][INDICATOR_MINIMUM].fDistance + ExtremePointMatrix[i][INDICATOR_MAXIMUM].fDistance) / 2.0F + (rfPadding + 0.5F) * fDifference;
      }

      float fTransformedCorner[ELEMENT_COUNT_3D][EXTREME_ELEMENT_COUNT];

      for (int i = 0; i < ELEMENT_COUNT_3D; i++) {
        fTransformedCorner[i][INDICATOR_MINIMUM] = fTransformedCorner[i][INDICATOR_MAXIMUM] = 0.0F;

        for (int j = 0; j < ELEMENT_COUNT_3D; j++)
          if (rpfTransformationMatrix[i][j] * fPaddedCorner[j][INDICATOR_MINIMUM] < rpfTransformationMatrix[i][j] * fPaddedCorner[j][INDICATOR_MAXIMUM]) {
            fTransformedCorner[i][INDICATOR_MINIMUM] += rpfTransformationMatrix[i][j] * fPaddedCorner[j][INDICATOR_MINIMUM];
            fTransformedCorner[i][INDICATOR_MAXIMUM] += rpfTransformationMatrix[i][j] * fPaddedCorner[j][INDICATOR_MAXIMUM];
          } else {
            fTransformedCorner[i][INDICATOR_MINIMUM] += rpfTransformationMatrix[i][j] * fPaddedCorner[j][INDICATOR_MAXIMUM];
            fTransformedCorner[i][INDICATOR_MAXIMUM] += rpfTransformationMatrix[i][j] * fPaddedCorner[j][INDICATOR_MINIMUM];
          }

        fTransformedCorner[i][INDICATOR_MINIMUM] += rpfTransformationMatrix[i][ELEMENT_COUNT_3D];
        fTransformedCorner[i][INDICATOR_MAXIMUM] += rpfTransformationMatrix[i][ELEMENT_COUNT_3D];
      }

      float fStep[ELEMENT_COUNT_3D] = {(fTransformedCorner[INDICATOR_X][INDICATOR_MAXIMUM] - fTransformedCorner[INDICATOR_X][INDICATOR_MINIMUM]) / static_cast<float>(riDivisionCountX), (fTransformedCorner[INDICATOR_Y][INDICATOR_MAXIMUM] - fTransformedCorner[INDICATOR_Y][INDICATOR_MINIMUM]) / static_cast<float>(riDivisionCountY), (fTransformedCorner[INDICATOR_Z][INDICATOR_MAXIMUM] - fTransformedCorner[INDICATOR_Z][INDICATOR_MINIMUM]) / static_cast<float>(riDivisionCountZ)};

      for (int i = 0; i < riDivisionCountX; i++) {
        for (int j = 0; j < riDivisionCountY; j++) {
          for (int k = 0; k < riDivisionCountZ; k++) {
            float fSamplingPoint[ELEMENT_COUNT_3D] = {fTransformedCorner[INDICATOR_X][INDICATOR_MINIMUM] + fStep[INDICATOR_X] * (static_cast<float>(i) + 0.5F), fTransformedCorner[INDICATOR_Y][INDICATOR_MINIMUM] + fStep[INDICATOR_Y] * (static_cast<float>(j) + 0.5F), fTransformedCorner[INDICATOR_Z][INDICATOR_MINIMUM] + fStep[INDICATOR_Z] * (static_cast<float>(k) + 0.5F)}, fTransformedSamplingPoint[ELEMENT_COUNT_3D] = {0.0F, 0.0F, 0.0F};

            TransformCoordinate<float, ELEMENT_COUNT_3D + 1>(fInversedTransformationMatrix, fSamplingPoint, fTransformedSamplingPoint);

            ofsFile << fSamplingPoint[INDICATOR_X] << ' ' << fSamplingPoint[INDICATOR_Y] << ' ' << fSamplingPoint[INDICATOR_Z] << ' ' << FunctionValue3D(fTransformedSamplingPoint, rvecVertex, rvecCoefficient) << endl;
          }

          if (riDivisionCountY - j - 1)
            ofsFile << endl;
        }

        if (riDivisionCountX - i - 1)
          ofsFile << endl;
      }
    }
  }

  ofsFile.close();
}

void ExportDataFile3D(string strFileName, const PiecewisePolynomial3D& rPiecewisePolynomial3DData, const float (*const& rpfTransformationMatrix)[ELEMENT_COUNT_3D + 1], const int& riDivisionCountX, const int& riDivisionCountY, const int& riDivisionCountZ, const Vertex3D& rVertex3DMinimumCorner, const Vertex3D& rVertex3DMaximumCorner)
{
  ofstream ofsFile;

  if (strFileName.empty()) {
    cout << "書き出すファイル名を入力してください：" << flush;

    getline(cin, strFileName);

    cout << endl;
  }

  ofsFile.open(strFileName + DATA_FILE_EXTENSION);

  if (!ofsFile.is_open()) {
    cerr << "ファイル「" << strFileName << DATA_FILE_EXTENSION << "」を開けませんでした。" << endl;

    exit(EXIT_FAILURE);
  }

  if (rpfTransformationMatrix == nullptr) {
    if (rVertex3DMinimumCorner.fX < rVertex3DMaximumCorner.fX && rVertex3DMinimumCorner.fY < rVertex3DMaximumCorner.fY && rVertex3DMinimumCorner.fZ < rVertex3DMaximumCorner.fZ) {
      float fStep[ELEMENT_COUNT_3D] = {(rVertex3DMaximumCorner.fX - rVertex3DMinimumCorner.fX) / static_cast<float>(riDivisionCountX), (rVertex3DMaximumCorner.fY - rVertex3DMinimumCorner.fY) / static_cast<float>(riDivisionCountY), (rVertex3DMaximumCorner.fZ - rVertex3DMinimumCorner.fZ) / static_cast<float>(riDivisionCountZ)};

      for (int i = 0; i < riDivisionCountX; i++) {
        for (int j = 0; j < riDivisionCountY; j++) {
          for (int k = 0; k < riDivisionCountZ; k++) {
            Vertex3D Vertex3DSamplingPoint{rVertex3DMinimumCorner.fX + fStep[INDICATOR_X] * (static_cast<float>(i) + 0.5F), rVertex3DMinimumCorner.fY + fStep[INDICATOR_Y] * (static_cast<float>(j) + 0.5F), rVertex3DMinimumCorner.fZ + fStep[INDICATOR_Z] * (static_cast<float>(k) + 0.5F)};

            ofsFile << Vertex3DSamplingPoint.fX << ' ' << Vertex3DSamplingPoint.fY << ' ' << Vertex3DSamplingPoint.fZ << ' ' << rPiecewisePolynomial3DData.At(Vertex3DSamplingPoint) << endl;
          }

          if (riDivisionCountY - j - 1)
            ofsFile << endl;
        }

        if (riDivisionCountX - i - 1)
          ofsFile << endl;
      }
    } else {
      float fStep[ELEMENT_COUNT_3D] = {(rPiecewisePolynomial3DData.MaximumCorner().fX - rPiecewisePolynomial3DData.MinimumCorner().fX) / static_cast<float>(riDivisionCountX), (rPiecewisePolynomial3DData.MaximumCorner().fY - rPiecewisePolynomial3DData.MinimumCorner().fY) / static_cast<float>(riDivisionCountY), (rPiecewisePolynomial3DData.MaximumCorner().fZ - rPiecewisePolynomial3DData.MinimumCorner().fZ) / static_cast<float>(riDivisionCountZ)};

      for (int i = 0; i < riDivisionCountX; i++) {
        for (int j = 0; j < riDivisionCountY; j++) {
          for (int k = 0; k < riDivisionCountZ; k++) {
            Vertex3D Vertex3DSamplingPoint{rPiecewisePolynomial3DData.MinimumCorner().fX + fStep[INDICATOR_X] * (static_cast<float>(i) + 0.5F), rPiecewisePolynomial3DData.MinimumCorner().fY + fStep[INDICATOR_Y] * (static_cast<float>(j) + 0.5F), rPiecewisePolynomial3DData.MinimumCorner().fZ + fStep[INDICATOR_Z] * (static_cast<float>(k) + 0.5F)};

            ofsFile << Vertex3DSamplingPoint.fX << ' ' << Vertex3DSamplingPoint.fY << ' ' << Vertex3DSamplingPoint.fZ << ' ' << rPiecewisePolynomial3DData.At(Vertex3DSamplingPoint) << endl;
          }

          if (riDivisionCountY - j - 1)
            ofsFile << endl;
        }

        if (riDivisionCountX - i - 1)
        ofsFile << endl;
      }
    }
  } else {
    float fInversedTransformationMatrix[ELEMENT_COUNT_3D + 1][ELEMENT_COUNT_3D + 1];

    for (int i = 0; i <= ELEMENT_COUNT_3D; i++)
      for (int j = 0; j <= ELEMENT_COUNT_3D; j++)
        fInversedTransformationMatrix[i][j] = rpfTransformationMatrix[i][j];

    InverseTransformation<ELEMENT_COUNT_3D + 1>(fInversedTransformationMatrix);

    if (rVertex3DMinimumCorner.fX < rVertex3DMaximumCorner.fX && rVertex3DMinimumCorner.fY < rVertex3DMaximumCorner.fY && rVertex3DMinimumCorner.fZ < rVertex3DMaximumCorner.fZ) {
      float fStep[ELEMENT_COUNT_3D] = {(rVertex3DMaximumCorner.fX - rVertex3DMinimumCorner.fX) / static_cast<float>(riDivisionCountX), (rVertex3DMaximumCorner.fY - rVertex3DMinimumCorner.fY) / static_cast<float>(riDivisionCountY), (rVertex3DMaximumCorner.fZ - rVertex3DMinimumCorner.fZ) / static_cast<float>(riDivisionCountZ)};

      for (int i = 0; i < riDivisionCountX; i++) {
        for (int j = 0; j < riDivisionCountY; j++) {
          for (int k = 0; k < riDivisionCountZ; k++) {
            Vertex3D Vertex3DSamplingPoint{rVertex3DMinimumCorner.fX + fStep[INDICATOR_X] * (static_cast<float>(i) + 0.5F), rVertex3DMinimumCorner.fY + fStep[INDICATOR_Y] * (static_cast<float>(j) + 0.5F), rVertex3DMinimumCorner.fZ + fStep[INDICATOR_Z] * (static_cast<float>(k) + 0.5F)};

            ofsFile << Vertex3DSamplingPoint.fX << ' ' << Vertex3DSamplingPoint.fY << ' ' << Vertex3DSamplingPoint.fZ << ' ' << rPiecewisePolynomial3DData.At(TransformCoordinate(fInversedTransformationMatrix, Vertex3DSamplingPoint)) << endl;
          }

          if (riDivisionCountY - j - 1)
            ofsFile << endl;
        }

        if (riDivisionCountX - i - 1)
          ofsFile << endl;
      }
    } else {
      float Vertex3D::*pfVertex3DMember[ELEMENT_COUNT_3D] = {&Vertex3D::fX, &Vertex3D::fY, &Vertex3D::fZ};
      vector<Vertex3D> vecOriginalCorner{rPiecewisePolynomial3DData.MinimumCorner(), rPiecewisePolynomial3DData.MaximumCorner()}, vecTransformedCorner(ELEMENT_COUNT_3D);

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

      float fStep[ELEMENT_COUNT_3D] = {(vecTransformedCorner[INDICATOR_MAXIMUM].fX - vecTransformedCorner[INDICATOR_MINIMUM].fX) / static_cast<float>(riDivisionCountX), (vecTransformedCorner[INDICATOR_MAXIMUM].fY - vecTransformedCorner[INDICATOR_MINIMUM].fY) / static_cast<float>(riDivisionCountY), (vecTransformedCorner[INDICATOR_MAXIMUM].fZ - vecTransformedCorner[INDICATOR_MINIMUM].fZ) / static_cast<float>(riDivisionCountZ)};

      for (int i = 0; i < riDivisionCountX; i++) {
        for (int j = 0; j < riDivisionCountY; j++) {
          for (int k = 0; k < riDivisionCountZ; k++) {
            Vertex3D Vertex3DSamplingPoint{vecTransformedCorner[INDICATOR_MINIMUM].fX + fStep[INDICATOR_X] * (static_cast<float>(i) + 0.5F), vecTransformedCorner[INDICATOR_MINIMUM].fY + fStep[INDICATOR_Y] * (static_cast<float>(j) + 0.5F), vecTransformedCorner[INDICATOR_MINIMUM].fZ + fStep[INDICATOR_Z] * (static_cast<float>(k) + 0.5F)};

            ofsFile << Vertex3DSamplingPoint.fX << ' ' << Vertex3DSamplingPoint.fY << ' ' << Vertex3DSamplingPoint.fZ << ' ' << rPiecewisePolynomial3DData.At(TransformCoordinate(fInversedTransformationMatrix, Vertex3DSamplingPoint)) << endl;
          }

          if (riDivisionCountY - j - 1)
            ofsFile << endl;
        }

        if (riDivisionCountX - i - 1)
        ofsFile << endl;
      }
    }
  }

  ofsFile.close();
}

void ExportPartialConvexHull2D(const string& rstrFileName, const ConvexHull2D& rConvexHull2DData)
{
  ostringstream ossData(rstrFileName);
  ios::fmtflags othFlag = cout.flags();
  char chCharacter = cout.fill();
  string::size_type othPosition = rstrFileName.rfind(DELIMITER_CHARACTER);
  struct stat othBuffer;

  ossData.seekp(0, ios::end);

  if (stat(ossData.str().c_str(), &othBuffer))
#if defined(_WIN32) || defined(_WIN64)
    _mkdir(ossData.str().c_str());
#elif defined (__linux) || defined(__MACH__)
    mkdir(ossData.str().c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
#else
    ;
#endif

  if (othPosition == string::npos)
    ossData << DELIMITER_CHARACTER << rstrFileName;
  else
    ossData << rstrFileName.substr(othPosition);

  ossData << " 凸包 " << setw(DIGIT) << setfill(FILL_CHARACTER) << rConvexHull2DData.VertexCount() - TRIANGLE_ELEMENT_COUNT + 1;
  cout << setiosflags(othFlag) << setfill(chCharacter);

  ExportObjectFile2D(ossData.str(), rConvexHull2DData);
}

void ExportPartialConvexHull3D(const string& rstrFileName, const ConvexHull3D& rConvexHull3DData)
{
  ostringstream ossData(rstrFileName);
  ios::fmtflags othFlag = cout.flags();
  char chCharacter = cout.fill();
  string::size_type othPosition = rstrFileName.rfind(DELIMITER_CHARACTER);
  struct stat othBuffer;

  ossData.seekp(0, ios::end);

  if (stat(ossData.str().c_str(), &othBuffer))
#if defined(_WIN32) || defined(_WIN64)
    _mkdir(ossData.str().c_str());
#elif defined(__linux) || defined(__MACH__)
    mkdir(ossData.str().c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
#else
    ;
#endif

  if (othPosition == string::npos)
    ossData << DELIMITER_CHARACTER << rstrFileName;
  else
    ossData << rstrFileName.substr(othPosition);

  ossData << " 凸包 ";

  switch (rConvexHull3DData.ShapeState()) {
  case Plane :
    ossData << setw(DIGIT) << setfill(FILL_CHARACTER) << rConvexHull3DData.VertexCount() - TRIANGLE_ELEMENT_COUNT + 1;

    break;
  case Solid :
    ossData << setw(DIGIT) << setfill(FILL_CHARACTER) << rConvexHull3DData.VertexCount() - TETRAHEDRON_ELEMENT_COUNT + 1;

    break;
  default : break;
  }

  cout << setiosflags(othFlag) << setfill(chCharacter);

  ExportObjectFile3D(ossData.str(), rConvexHull3DData);
}