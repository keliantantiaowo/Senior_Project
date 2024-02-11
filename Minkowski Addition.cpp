#include "Geometric Operation.h"
#include "Minkowski Addition.h"

using namespace std;

ConvexHull2D MinkowskiSum2D(const ConvexHull2D& rConvexHull2DDataA, const ConvexHull2D& rConvexHull2DDataB, const float (*const& rpfTransformationMatrixA)[ELEMENT_COUNT_2D + 1], const float (*const& rpfTransformationMatrixB)[ELEMENT_COUNT_2D + 1])
{
  vector<Vertex2D> vecTransformedVertexA(rConvexHull2DDataA.VertexCount(), Vertex2D{0.0F, 0.0F}), vecTransformedVertexB(rConvexHull2DDataB.VertexCount(), Vertex2D{0.0F, 0.0F});
  list<Vertex2D>::const_iterator iterPositionA = rConvexHull2DDataA.Begin(), iterPositionB = rConvexHull2DDataB.Begin();

  for (int i = 0; i < vecTransformedVertexA.size(); i++)
    vecTransformedVertexA[i] = rpfTransformationMatrixA == nullptr ? *iterPositionA++ : TransformCoordinate(rpfTransformationMatrixA, *iterPositionA++);

  for (int i = 0; i < vecTransformedVertexB.size(); i++)
    vecTransformedVertexB[i] = rpfTransformationMatrixB == nullptr ? *iterPositionB++ : TransformCoordinate(rpfTransformationMatrixB, *iterPositionB++);

  vector<Vertex2D> vecComposedVertex(vecTransformedVertexA.size() * vecTransformedVertexB.size(), Vertex2D{0.0F, 0.0F});

  for (int i = 0; i < vecTransformedVertexA.size(); i++)
    for (int j = 0; j < vecTransformedVertexB.size(); j++)
      vecComposedVertex[vecTransformedVertexB.size() * i + j] = Vertex2D{vecTransformedVertexA[i].fX + vecTransformedVertexB[j].fX, vecTransformedVertexA[i].fY + vecTransformedVertexB[j].fY};

  return ConvexHull2D(vecComposedVertex);
}

ConvexHull3D MinkowskiSum3D(const ConvexHull3D& rConvexHull3DDataA, const ConvexHull3D& rConvexHull3DDataB, const float (*const& rpfTransformationMatrixA)[ELEMENT_COUNT_3D + 1], const float (*const& rpfTransformationMatrixB)[ELEMENT_COUNT_3D + 1])
{
  vector<Vertex3D> vecTransformedVertexA(rConvexHull3DDataA.VertexCount(), Vertex3D{0.0F, 0.0F, 0.0F}), vecTransformedVertexB(rConvexHull3DDataB.VertexCount(), Vertex3D{0.0F, 0.0F, 0.0F});
  list<Vertex3D>::const_iterator iterPositionA = rConvexHull3DDataA.VertexBegin(), iterPositionB = rConvexHull3DDataB.VertexBegin();

  for (int i = 0; i < vecTransformedVertexA.size(); i++)
    vecTransformedVertexA[i] = rpfTransformationMatrixA == nullptr ? *iterPositionA++ : TransformCoordinate(rpfTransformationMatrixA, *iterPositionA++);

  for (int i = 0; i < vecTransformedVertexB.size(); i++)
    vecTransformedVertexB[i] = rpfTransformationMatrixB == nullptr ? *iterPositionB++ : TransformCoordinate(rpfTransformationMatrixB, *iterPositionB++);

  vector<Vertex3D> vecComposedVertex(vecTransformedVertexA.size() * vecTransformedVertexB.size(), Vertex3D{0.0F, 0.0F, 0.0F});

  for (int i = 0; i < vecTransformedVertexA.size(); i++)
    for (int j = 0; j < vecTransformedVertexB.size(); j++)
      vecComposedVertex[vecTransformedVertexB.size() * i + j] = Vertex3D{vecTransformedVertexA[i].fX + vecTransformedVertexB[j].fX, vecTransformedVertexA[i].fY + vecTransformedVertexB[j].fY, vecTransformedVertexA[i].fZ + vecTransformedVertexB[j].fZ};

  return ConvexHull3D(vecComposedVertex);
}

ConvexHull2D MinkowskiDifference2D(const ConvexHull2D& rConvexHull2DDataA, const ConvexHull2D& rConvexHull2DDataB, const float (*const& rpfTransformationMatrixA)[ELEMENT_COUNT_2D + 1], const float (*const& rpfTransformationMatrixB)[ELEMENT_COUNT_2D + 1])
{
  vector<Vertex2D> vecTransformedVertexA(rConvexHull2DDataA.VertexCount(), Vertex2D{0.0F, 0.0F}), vecTransformedVertexB(rConvexHull2DDataB.VertexCount(), Vertex2D{0.0F, 0.0F});
  list<Vertex2D>::const_iterator iterPositionA = rConvexHull2DDataA.Begin(), iterPositionB = rConvexHull2DDataB.Begin();

  for (int i = 0; i < vecTransformedVertexA.size(); i++)
    vecTransformedVertexA[i] = rpfTransformationMatrixA == nullptr ? *iterPositionA++ : TransformCoordinate(rpfTransformationMatrixA, *iterPositionA++);

  for (int i = 0; i < vecTransformedVertexB.size(); i++)
    vecTransformedVertexB[i] = rpfTransformationMatrixB == nullptr ? *iterPositionB++ : TransformCoordinate(rpfTransformationMatrixB, *iterPositionB++);

  vector<Vertex2D> vecComposedVertex(vecTransformedVertexA.size() * vecTransformedVertexB.size(), Vertex2D{0.0F, 0.0F});

  for (int i = 0; i < vecTransformedVertexA.size(); i++)
    for (int j = 0; j < vecTransformedVertexB.size(); j++)
      vecComposedVertex[vecTransformedVertexB.size() * i + j] = Vertex2D{vecTransformedVertexA[i].fX - vecTransformedVertexB[j].fX, vecTransformedVertexA[i].fY - vecTransformedVertexB[j].fY};

  return ConvexHull2D(vecComposedVertex);
}

ConvexHull3D MinkowskiDifference3D(const ConvexHull3D& rConvexHull3DDataA, const ConvexHull3D& rConvexHull3DDataB, const float (*const& rpfTransformationMatrixA)[ELEMENT_COUNT_3D + 1], const float (*const& rpfTransformationMatrixB)[ELEMENT_COUNT_3D + 1])
{
  vector<Vertex3D> vecTransformedVertexA(rConvexHull3DDataA.VertexCount(), Vertex3D{0.0F, 0.0F, 0.0F}), vecTransformedVertexB(rConvexHull3DDataB.VertexCount(), Vertex3D{0.0F, 0.0F, 0.0F});
  list<Vertex3D>::const_iterator iterPositionA = rConvexHull3DDataA.VertexBegin(), iterPositionB = rConvexHull3DDataB.VertexBegin();

  for (int i = 0; i < vecTransformedVertexA.size(); i++)
    vecTransformedVertexA[i] = rpfTransformationMatrixA == nullptr ? *iterPositionA++ : TransformCoordinate(rpfTransformationMatrixA, *iterPositionA++);

  for (int i = 0; i < vecTransformedVertexB.size(); i++)
    vecTransformedVertexB[i] = rpfTransformationMatrixB == nullptr ? *iterPositionB++ : TransformCoordinate(rpfTransformationMatrixB, *iterPositionB++);

  vector<Vertex3D> vecComposedVertex(vecTransformedVertexA.size() * vecTransformedVertexB.size(), Vertex3D{0.0F, 0.0F, 0.0F});

  for (int i = 0; i < vecTransformedVertexA.size(); i++)
    for (int j = 0; j < vecTransformedVertexB.size(); j++)
      vecComposedVertex[vecTransformedVertexB.size() * i + j] = Vertex3D{vecTransformedVertexA[i].fX - vecTransformedVertexB[j].fX, vecTransformedVertexA[i].fY - vecTransformedVertexB[j].fY, vecTransformedVertexA[i].fZ - vecTransformedVertexB[j].fZ};

  return ConvexHull3D(vecComposedVertex);
}