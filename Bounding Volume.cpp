#include <algorithm>
#include <iterator>
#include <limits>
#include "Geometric Operation.h"
#include "Bounding Volume.h"

using namespace std;

AABB2D::AABB2D(const vector<Vertex2D>& rvecVertex) : Vertex2DCenter(), RangeInfo2DData()
{
  float fMinimumX, fMaximumX, fMinimumY, fMaximumY;

  fMinimumX = fMinimumY = numeric_limits<float>::max();
  fMaximumX = fMaximumY = numeric_limits<float>::lowest();

  for (int i = 0; i < rvecVertex.size(); i++) {
    if (rvecVertex[i].fX < fMinimumX)
      fMinimumX = rvecVertex[i].fX;

    if (rvecVertex[i].fX > fMaximumX)
      fMaximumX = rvecVertex[i].fX;

    if (rvecVertex[i].fY < fMinimumY)
      fMinimumY = rvecVertex[i].fY;

    if (rvecVertex[i].fY > fMaximumY)
      fMaximumY = rvecVertex[i].fY;
  }

  Vertex2DCenter = Vertex2D{(fMinimumX + fMaximumX) / 2.0F, (fMinimumY + fMaximumY) / 2.0F};
  RangeInfo2DData = RangeInfo2D{(fMaximumX - fMinimumX) / 2.0F, (fMaximumY - fMinimumY) / 2.0F};
}

AABB2D::AABB2D(const AABB2D& rAABB2DData) : Vertex2DCenter(rAABB2DData.Vertex2DCenter), RangeInfo2DData(rAABB2DData.RangeInfo2DData)
{
}

AABB2D::AABB2D(AABB2D&& rAABB2DData) : Vertex2DCenter(move(rAABB2DData.Vertex2DCenter)), RangeInfo2DData(move(rAABB2DData.RangeInfo2DData))
{
}

AABB2D::~AABB2D(void)
{
}

AABB2D& AABB2D::operator=(const AABB2D& rAABB2DData)
{
  Vertex2DCenter = rAABB2DData.Vertex2DCenter;
  RangeInfo2DData = rAABB2DData.RangeInfo2DData;

  return *this;
}

AABB2D& AABB2D::operator=(AABB2D&& rAABB2DData)
{
  Vertex2DCenter = move(rAABB2DData.Vertex2DCenter);
  RangeInfo2DData = move(rAABB2DData.RangeInfo2DData);

  return *this;
}

Vertex2D AABB2D::Center(void) const
{
  return Vertex2DCenter;
}

AABB2D::RangeInfo2D AABB2D::Size(void) const
{
  return RangeInfo2DData;
}

AABB3D::AABB3D(const vector<Vertex3D>& rvecVertex) : Vertex3DCenter(), RangeInfo3DData()
{
  float fMinimumX, fMaximumX, fMinimumY, fMaximumY, fMinimumZ, fMaximumZ;

  fMinimumX = fMinimumY = fMinimumZ = numeric_limits<float>::max();
  fMaximumX = fMaximumY = fMaximumZ = numeric_limits<float>::lowest();

  for (int i = 0; i < rvecVertex.size(); i++) {
    if (rvecVertex[i].fX < fMinimumX)
      fMinimumX = rvecVertex[i].fX;

    if (rvecVertex[i].fX > fMaximumX)
      fMaximumX = rvecVertex[i].fX;

    if (rvecVertex[i].fY < fMinimumY)
      fMinimumY = rvecVertex[i].fY;

    if (rvecVertex[i].fY > fMaximumY)
      fMaximumY = rvecVertex[i].fY;

    if (rvecVertex[i].fZ < fMinimumZ)
      fMinimumZ = rvecVertex[i].fZ;

    if (rvecVertex[i].fZ > fMaximumZ)
      fMaximumZ = rvecVertex[i].fZ;
  }

  Vertex3DCenter = Vertex3D{(fMinimumX + fMaximumX) / 2.0F, (fMinimumY + fMaximumY) / 2.0F, (fMinimumZ + fMaximumZ) / 2.0F};
  RangeInfo3DData = RangeInfo3D{(fMaximumX - fMinimumX) / 2.0F, (fMaximumY - fMinimumY) / 2.0F, (fMaximumZ - fMinimumZ) / 2.0F};
}

AABB3D::AABB3D(const AABB3D& rAABB3DData) : Vertex3DCenter(rAABB3DData.Vertex3DCenter), RangeInfo3DData(rAABB3DData.RangeInfo3DData)
{
}

AABB3D::AABB3D(AABB3D&& rAABB3DData) : Vertex3DCenter(move(rAABB3DData.Vertex3DCenter)), RangeInfo3DData(move(rAABB3DData.RangeInfo3DData))
{
}

AABB3D::~AABB3D(void)
{
}

AABB3D& AABB3D::operator=(const AABB3D& rAABB3DData)
{
  Vertex3DCenter = rAABB3DData.Vertex3DCenter;
  RangeInfo3DData = rAABB3DData.RangeInfo3DData;

  return *this;
}

AABB3D& AABB3D::operator=(AABB3D&& rAABB3DData)
{
  Vertex3DCenter = move(rAABB3DData.Vertex3DCenter);
  RangeInfo3DData = move(rAABB3DData.RangeInfo3DData);

  return *this;
}

Vertex3D AABB3D::Center(void) const
{
  return Vertex3DCenter;
}

AABB3D::RangeInfo3D AABB3D::Size(void) const
{
  return RangeInfo3DData;
}

ConvexHull2D::ConvexHull2D(const vector<Vertex2D>& rvecVertex) : listHullVertex(0), fEpsilon(0.0F)
{
  if (rvecVertex.size() < TRIANGLE_ELEMENT_COUNT)
    throw invalid_argument("頂点数が不足するため、凸包を生成できません。");

  ExteriorSet ExteriorSetDataA{list<unsigned>(0), ExtremePoint{numeric_limits<unsigned>::max(), 0.0F}}, ExteriorSetDataB{list<unsigned>(0), ExtremePoint{numeric_limits<unsigned>::max(), 0.0F}};

  GenerateInitialHull(rvecVertex, ExteriorSetDataA, ExteriorSetDataB);

  if (ExteriorSetDataA.listIndexSet.empty() && ExteriorSetDataB.listIndexSet.empty())
    throw invalid_argument("頂点数が不足するため、凸包を生成できません。");

  if (ExteriorSetDataA.listIndexSet.size())
    ExpandHull(rvecVertex, listHullVertex.begin(), prev(listHullVertex.end()), ExteriorSetDataA);

  if (ExteriorSetDataB.listIndexSet.size())
    ExpandHull(rvecVertex, prev(listHullVertex.end()), listHullVertex.begin(), ExteriorSetDataB);
}

ConvexHull2D::ConvexHull2D(const ConvexHull2D& rConvexHull2DData) : listHullVertex(rConvexHull2DData.listHullVertex), fEpsilon(rConvexHull2DData.fEpsilon)
{
}

ConvexHull2D::ConvexHull2D(ConvexHull2D&& rConvexHull2DData) : listHullVertex(move(rConvexHull2DData.listHullVertex)), fEpsilon(move(rConvexHull2DData.fEpsilon))
{
}

ConvexHull2D::~ConvexHull2D(void)
{
}

ConvexHull2D& ConvexHull2D::operator=(const ConvexHull2D& rConvexHull2DData)
{
  listHullVertex = rConvexHull2DData.listHullVertex;
  fEpsilon = rConvexHull2DData.fEpsilon;

  return *this;
}

ConvexHull2D& ConvexHull2D::operator=(ConvexHull2D&& rConvexHull2DData)
{
  listHullVertex = move(rConvexHull2DData.listHullVertex);
  fEpsilon = move(rConvexHull2DData.fEpsilon);

  return *this;
}

Vertex2D& ConvexHull2D::operator[](const int& riIndex)
{
  int iCount = 0;

  if (riIndex < 0 || riIndex > listHullVertex.size() - 1)
    throw range_error("インデックスが不正です。");

  for (list<Vertex2D>::iterator i = listHullVertex.begin(); i != listHullVertex.end(); i++)
    if (riIndex == iCount++)
      return *i;
}

const Vertex2D& ConvexHull2D::operator[](const int& riIndex) const
{
  int iCount = 0;

  if (riIndex < 0 || riIndex > listHullVertex.size() - 1)
    throw range_error("インデックスが不正です。");

  for (list<Vertex2D>::const_iterator i = listHullVertex.cbegin(); i != listHullVertex.cend(); i++)
    if (riIndex == iCount++)
      return *i;
}

list<Vertex2D>::iterator ConvexHull2D::Begin(void)
{
  return listHullVertex.begin();
}

list<Vertex2D>::const_iterator ConvexHull2D::Begin(void) const
{
  return listHullVertex.cbegin();
}

list<Vertex2D>::iterator ConvexHull2D::End(void)
{
  return listHullVertex.end();
}

list<Vertex2D>::const_iterator ConvexHull2D::End(void) const
{
  return listHullVertex.cend();
}

int ConvexHull2D::VertexCount(void) const
{
  return listHullVertex.size();
}

float ConvexHull2D::Epsilon(void) const
{
  return fEpsilon;
}

void ConvexHull2D::ApplyTransformation(const float (*const& rpfTransformationMatrix)[ELEMENT_COUNT_2D + 1])
{
  for (Vertex2D& rVertex2DData : listHullVertex)
    rVertex2DData = TransformCoordinate(rpfTransformationMatrix, rVertex2DData);
}

void ConvexHull2D::GenerateInitialHull(const vector<Vertex2D>& rvecVertex, ExteriorSet& rExteriorSetDataA, ExteriorSet& rExteriorSetDataB)
{
  unsigned uIndexArray[EXTREME_ELEMENT_COUNT];

  FindExtremePoint(rvecVertex, uIndexArray);

  if (Distance(rvecVertex[uIndexArray[INDICATOR_MINIMUM]], rvecVertex[uIndexArray[INDICATOR_MAXIMUM]]) < numeric_limits<float>::epsilon())
    throw invalid_argument("ポリゴンが点に縮退しています。");

  fEpsilon = ((fabs(rvecVertex[uIndexArray[INDICATOR_MINIMUM]].fX) + fabs(rvecVertex[uIndexArray[INDICATOR_MAXIMUM]].fX)) + (fabs(rvecVertex[uIndexArray[INDICATOR_MINIMUM]].fY) + fabs(rvecVertex[uIndexArray[INDICATOR_MAXIMUM]].fY))) * numeric_limits<float>::epsilon();

  for (int i = 0; i < EXTREME_ELEMENT_COUNT; i++)
    listHullVertex.push_back(rvecVertex[uIndexArray[i]]);

  for (int i = 0; i < rvecVertex.size(); i++) {
    float fDistance = Distance(rvecVertex[i], *listHullVertex.begin(), *prev(listHullVertex.end()));

    if (fDistance > fEpsilon) {
      rExteriorSetDataA.listIndexSet.push_back(static_cast<unsigned>(i));

      if (rExteriorSetDataA.ExtremePointData.fDistance < fDistance) {
        rExteriorSetDataA.ExtremePointData.uIndex = static_cast<unsigned>(i);
        rExteriorSetDataA.ExtremePointData.fDistance = fDistance;
      }
    } else if (fDistance < -fEpsilon) {
      rExteriorSetDataB.listIndexSet.push_back(static_cast<unsigned>(i));

      if (rExteriorSetDataB.ExtremePointData.fDistance < -fDistance) {
        rExteriorSetDataB.ExtremePointData.uIndex = static_cast<unsigned>(i);
        rExteriorSetDataB.ExtremePointData.fDistance = -fDistance;
      }
    } else
      ;
  }

  rExteriorSetDataA.ExtremePointData.fDistance = rExteriorSetDataB.ExtremePointData.fDistance = 0.0F;
}

void ConvexHull2D::FindExtremePoint(const vector<Vertex2D>& rvecVertex, unsigned* const& rpuIndexArray) const
{
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

  if (ExtremePointMatrix[INDICATOR_X][INDICATOR_MAXIMUM].fDistance - ExtremePointMatrix[INDICATOR_X][INDICATOR_MINIMUM].fDistance >= ExtremePointMatrix[INDICATOR_Y][INDICATOR_MAXIMUM].fDistance - ExtremePointMatrix[INDICATOR_Y][INDICATOR_MINIMUM].fDistance) {
    rpuIndexArray[INDICATOR_MINIMUM] = ExtremePointMatrix[INDICATOR_X][INDICATOR_MINIMUM].uIndex;
    rpuIndexArray[INDICATOR_MAXIMUM] = ExtremePointMatrix[INDICATOR_X][INDICATOR_MAXIMUM].uIndex;
  } else {
    rpuIndexArray[INDICATOR_MINIMUM] = ExtremePointMatrix[INDICATOR_Y][INDICATOR_MINIMUM].uIndex;
    rpuIndexArray[INDICATOR_MAXIMUM] = ExtremePointMatrix[INDICATOR_Y][INDICATOR_MAXIMUM].uIndex;
  }
}

void ConvexHull2D::ExpandHull(const vector<Vertex2D>& rvecVertex, const list<Vertex2D>::iterator& riterStart, const list<Vertex2D>::iterator& riterEnd, ExteriorSet& rExteriorSetData)
{
  list<Vertex2D>::iterator iterIntermediate = listHullVertex.insert(next(riterStart), rvecVertex[rExteriorSetData.ExtremePointData.uIndex]);
  ExteriorSet ExteriorSetAddition{list<unsigned>(0), ExtremePoint{numeric_limits<unsigned>::max(), 0.0F}};
  Vertex2D Vertex2DEdgeNormalA, Vertex2DEdgeNormalB;

  Vertex2DEdgeNormalA = EdgeNormal(Vertex2D{rvecVertex[rExteriorSetData.ExtremePointData.uIndex].fX - riterStart->fX, rvecVertex[rExteriorSetData.ExtremePointData.uIndex].fY - riterStart->fY});
  Vertex2DEdgeNormalB = EdgeNormal(Vertex2D{riterEnd->fX - rvecVertex[rExteriorSetData.ExtremePointData.uIndex].fX, riterEnd->fY - rvecVertex[rExteriorSetData.ExtremePointData.uIndex].fY});

  for (list<unsigned>::iterator i = rExteriorSetData.listIndexSet.begin(); i != rExteriorSetData.listIndexSet.end(); true) {
    float fDistance = 0.0F;

    if (fDistance = DotProduct(Vertex2D{rvecVertex[*i].fX - riterStart->fX, rvecVertex[*i].fY - riterStart->fY}, Vertex2DEdgeNormalA), fDistance > fEpsilon) {
      if (rExteriorSetData.ExtremePointData.fDistance < fDistance) {
        rExteriorSetData.ExtremePointData.uIndex = *i;
        rExteriorSetData.ExtremePointData.fDistance = fDistance;
      }

      i++;
    } else if (fDistance = DotProduct(Vertex2D{rvecVertex[*i].fX - riterEnd->fX, rvecVertex[*i].fY - riterEnd->fY}, Vertex2DEdgeNormalB), fDistance > fEpsilon) {
      ExteriorSetAddition.listIndexSet.push_back(*i);

      if (ExteriorSetAddition.ExtremePointData.fDistance < fDistance) {
        ExteriorSetAddition.ExtremePointData.uIndex = *i;
        ExteriorSetAddition.ExtremePointData.fDistance = fDistance;
      }

      i = rExteriorSetData.listIndexSet.erase(i);
    } else
      i = rExteriorSetData.listIndexSet.erase(i);
  }

  rExteriorSetData.ExtremePointData.fDistance = ExteriorSetAddition.ExtremePointData.fDistance = 0.0F;

  if (rExteriorSetData.listIndexSet.size())
    ExpandHull(rvecVertex, riterStart, iterIntermediate, rExteriorSetData);

  if (ExteriorSetAddition.listIndexSet.size())
    ExpandHull(rvecVertex, iterIntermediate, riterEnd, ExteriorSetAddition);
}

ConvexHull3D::ConvexHull3D(const vector<Vertex3D>& rvecVertex) : listHullVertex(0), listHullSurface(0), ShapeTypeData(Undefined), fEpsilon(0.0F)
{
  if (rvecVertex.size() < TRIANGLE_ELEMENT_COUNT)
    throw invalid_argument("頂点数が不足するため、凸包を生成できません。");

  GenerateInitialHull(rvecVertex);

  if (listHullVertex.size() < TETRAHEDRON_ELEMENT_COUNT) {
    ShapeTypeData = Plane;

    ExteriorSet ExteriorSetDataA{list<unsigned>(0), ExtremePoint{numeric_limits<unsigned>::max(), 0.0F}}, ExteriorSetDataB{list<unsigned>(0), ExtremePoint{numeric_limits<unsigned>::max(), 0.0F}}, ExteriorSetDataC{list<unsigned>(0), ExtremePoint{numeric_limits<unsigned>::max(), 0.0F}};

    listHullSurface.pop_back();

    InitializeExteriorSet(rvecVertex, ExteriorSetDataA, ExteriorSetDataB, ExteriorSetDataC);

    list<Vertex3D>::iterator iterPosition = next(listHullVertex.begin());

    if (ExteriorSetDataA.listIndexSet.size())
      ExpandHull(rvecVertex, listHullVertex.begin(), iterPosition, ExteriorSetDataA);

    if (ExteriorSetDataB.listIndexSet.size())
      ExpandHull(rvecVertex, iterPosition, prev(listHullVertex.end()), ExteriorSetDataB);

    if (ExteriorSetDataC.listIndexSet.size())
      ExpandHull(rvecVertex, prev(listHullVertex.end()), listHullVertex.begin(), ExteriorSetDataC);
  } else {
    ShapeTypeData = Solid;

    list<Surface3D>::iterator iterPosition = listHullSurface.begin();

    do
      if (iterPosition->ExteriorSetData.listIndexSet.size())
        ExpandHull(rvecVertex, iterPosition);
      else
        iterPosition++;
    while (iterPosition != listHullSurface.end());
  }
}

ConvexHull3D::ConvexHull3D(const ConvexHull3D& rConvexHull3DData) : listHullVertex(rConvexHull3DData.listHullVertex), listHullSurface(rConvexHull3DData.listHullSurface), ShapeTypeData(rConvexHull3DData.ShapeTypeData), fEpsilon(rConvexHull3DData.fEpsilon)
{
  CopySurface(rConvexHull3DData);
}

ConvexHull3D::ConvexHull3D(ConvexHull3D&& rConvexHull3DData) : listHullVertex(move(rConvexHull3DData.listHullVertex)), listHullSurface(move(rConvexHull3DData.listHullSurface)), ShapeTypeData(move(rConvexHull3DData.ShapeTypeData)), fEpsilon(move(rConvexHull3DData.fEpsilon))
{
  MoveSurface(move(rConvexHull3DData));
}

ConvexHull3D::~ConvexHull3D(void)
{
  DeleteSurface();
}

ConvexHull3D& ConvexHull3D::operator=(const ConvexHull3D& rConvexHull3DData)
{
  DeleteSurface();

  listHullVertex = rConvexHull3DData.listHullVertex;
  listHullSurface = rConvexHull3DData.listHullSurface;
  fEpsilon = rConvexHull3DData.fEpsilon;

  CopySurface(rConvexHull3DData);

  return *this;
}

ConvexHull3D& ConvexHull3D::operator=(ConvexHull3D&& rConvexHull3DData)
{
  DeleteSurface();

  listHullVertex = move(rConvexHull3DData.listHullVertex);
  listHullSurface = move(rConvexHull3DData.listHullSurface);
  fEpsilon = move(rConvexHull3DData.fEpsilon);

  MoveSurface(move(rConvexHull3DData));

  return *this;
}

list<Vertex3D>::iterator ConvexHull3D::VertexBegin(void)
{
  return listHullVertex.begin();
}

list<Vertex3D>::const_iterator ConvexHull3D::VertexBegin(void) const
{
  return listHullVertex.cbegin();
}

list<Vertex3D>::iterator ConvexHull3D::VertexEnd(void)
{
  return listHullVertex.end();
}

list<Vertex3D>::const_iterator ConvexHull3D::VertexEnd(void) const
{
  return listHullVertex.cend();
}

list<Surface3D>::iterator ConvexHull3D::SurfaceBegin(void)
{
  return listHullSurface.begin();
}

list<Surface3D>::const_iterator ConvexHull3D::SurfaceBegin(void) const
{
  return listHullSurface.cbegin();
}

list<Surface3D>::iterator ConvexHull3D::SurfaceEnd(void)
{
  return listHullSurface.end();
}

list<Surface3D>::const_iterator ConvexHull3D::SurfaceEnd(void) const
{
  return listHullSurface.cend();
}

int ConvexHull3D::VertexCount(void) const
{
  return listHullVertex.size();
}

int ConvexHull3D::SurfaceCount(void) const
{
  return listHullSurface.size();
}

const ShapeType ConvexHull3D::ShapeState(void) const
{
  return ShapeTypeData;
}

float ConvexHull3D::Epsilon(void) const
{
  return fEpsilon;
}

void ConvexHull3D::ApplyTransformation(const float (*const& rpfTransformationMatrix)[ELEMENT_COUNT_3D + 1])
{
  for (Vertex3D& rVertex3DData : listHullVertex)
    rVertex3DData = TransformCoordinate(rpfTransformationMatrix, rVertex3DData);
}

void ConvexHull3D::GenerateInitialHull(const vector<Vertex3D>& rvecVertex)
{
  unsigned uIndexArray[EXTREME_ELEMENT_COUNT];

  FindExtremePoint(rvecVertex, uIndexArray);

  if (Distance(rvecVertex[uIndexArray[INDICATOR_MINIMUM]], rvecVertex[uIndexArray[INDICATOR_MAXIMUM]]) < numeric_limits<float>::epsilon())
    throw invalid_argument("ポリゴンが点に縮退しています。");

  fEpsilon = ((fabs(rvecVertex[uIndexArray[INDICATOR_MINIMUM]].fX) + fabs(rvecVertex[uIndexArray[INDICATOR_MAXIMUM]].fX)) + (fabs(rvecVertex[uIndexArray[INDICATOR_MINIMUM]].fY) + fabs(rvecVertex[uIndexArray[INDICATOR_MAXIMUM]].fY)) + (fabs(rvecVertex[uIndexArray[INDICATOR_MINIMUM]].fZ) + fabs(rvecVertex[uIndexArray[INDICATOR_MAXIMUM]].fZ))) * numeric_limits<float>::epsilon();

  for (int i = 0; i < EXTREME_ELEMENT_COUNT; i++)
    listHullVertex.push_back(rvecVertex[uIndexArray[i]]);

  GenerateInitialTriangle(rvecVertex);
  GenerateInitialTetrahedron(rvecVertex);

  if (listHullVertex.size() == TRIANGLE_ELEMENT_COUNT)
    ShapeTypeData = Plane;
  else if (listHullVertex.size() == TETRAHEDRON_ELEMENT_COUNT)
    ShapeTypeData = Solid;
  else
    ;
}

void ConvexHull3D::GenerateInitialTriangle(const vector<Vertex3D>& rvecVertex)
{
  ExtremePoint ExtremePointData{numeric_limits<unsigned>::max(), 0.0F};

  for (int i = 0; i < rvecVertex.size(); i++) {
    float fDistance = Distance(rvecVertex[i], *listHullVertex.begin(), *prev(listHullVertex.end()));

    if (ExtremePointData.fDistance < fDistance) {
      ExtremePointData.uIndex = static_cast<unsigned>(i);
      ExtremePointData.fDistance = fDistance;
    }
  }

  if (ExtremePointData.fDistance < fEpsilon)
    throw invalid_argument("頂点数が不足するため、凸包を生成できません。");

  listHullVertex.push_back(rvecVertex[ExtremePointData.uIndex]);

  for (int i = 0; i < INITIAL_SURFACE_COUNT; i++)
    listHullSurface.push_back(Surface3D{new HalfEdge3D, Vertex3D{0.0F, 0.0F, 0.0F}, ExteriorSet{list<unsigned>(0), ExtremePoint{numeric_limits<unsigned>::max(), 0.0F}}});

  HalfEdge3D (*pHalfEdge3DFront) = listHullSurface.begin()->pHalfEdge3DPointer, (*pHalfEdge3DBack) = prev(listHullSurface.end())->pHalfEdge3DPointer;

  pHalfEdge3DFront->pVertex3DPointer = &*listHullVertex.begin();
  pHalfEdge3DFront->pHalfEdge3DPrevious = new HalfEdge3D;
  pHalfEdge3DFront->pHalfEdge3DNext = new HalfEdge3D;
  pHalfEdge3DFront->pSurface3DPointer = &*listHullSurface.begin();
  pHalfEdge3DFront->pHalfEdge3DPrevious->pVertex3DPointer = &*prev(listHullVertex.end());
  pHalfEdge3DFront->pHalfEdge3DPrevious->pHalfEdge3DPrevious = pHalfEdge3DFront->pHalfEdge3DNext;
  pHalfEdge3DFront->pHalfEdge3DPrevious->pHalfEdge3DNext = pHalfEdge3DFront;
  pHalfEdge3DFront->pHalfEdge3DPrevious->pSurface3DPointer = &*listHullSurface.begin();
  pHalfEdge3DFront->pHalfEdge3DNext->pVertex3DPointer = &*next(listHullVertex.begin());
  pHalfEdge3DFront->pHalfEdge3DNext->pHalfEdge3DPrevious = pHalfEdge3DFront;
  pHalfEdge3DFront->pHalfEdge3DNext->pHalfEdge3DNext = pHalfEdge3DFront->pHalfEdge3DPrevious;
  pHalfEdge3DFront->pHalfEdge3DNext->pSurface3DPointer = &*listHullSurface.begin();
  pHalfEdge3DBack->pVertex3DPointer = &*next(listHullVertex.begin());
  pHalfEdge3DBack->pHalfEdge3DPrevious = new HalfEdge3D;
  pHalfEdge3DBack->pHalfEdge3DNext = new HalfEdge3D;
  pHalfEdge3DBack->pSurface3DPointer = &*prev(listHullSurface.end());
  pHalfEdge3DBack->pHalfEdge3DPrevious->pVertex3DPointer = &*prev(listHullVertex.end());
  pHalfEdge3DBack->pHalfEdge3DPrevious->pHalfEdge3DPrevious = pHalfEdge3DBack->pHalfEdge3DNext;
  pHalfEdge3DBack->pHalfEdge3DPrevious->pHalfEdge3DNext = pHalfEdge3DBack;
  pHalfEdge3DBack->pHalfEdge3DPrevious->pSurface3DPointer = &*prev(listHullSurface.end());
  pHalfEdge3DBack->pHalfEdge3DNext->pVertex3DPointer = &*listHullVertex.begin();
  pHalfEdge3DBack->pHalfEdge3DNext->pHalfEdge3DPrevious = pHalfEdge3DBack;
  pHalfEdge3DBack->pHalfEdge3DNext->pHalfEdge3DNext = pHalfEdge3DBack->pHalfEdge3DPrevious;
  pHalfEdge3DBack->pHalfEdge3DNext->pSurface3DPointer = &*prev(listHullSurface.end());
  pHalfEdge3DFront->pHalfEdge3DPair = pHalfEdge3DBack;
  pHalfEdge3DFront->pHalfEdge3DPrevious->pHalfEdge3DPair = pHalfEdge3DBack->pHalfEdge3DNext;
  pHalfEdge3DFront->pHalfEdge3DNext->pHalfEdge3DPair = pHalfEdge3DBack->pHalfEdge3DPrevious;
  pHalfEdge3DBack->pHalfEdge3DPair = pHalfEdge3DFront;
  pHalfEdge3DBack->pHalfEdge3DPrevious->pHalfEdge3DPair = pHalfEdge3DFront->pHalfEdge3DNext;
  pHalfEdge3DBack->pHalfEdge3DNext->pHalfEdge3DPair = pHalfEdge3DFront->pHalfEdge3DPrevious;
}

void ConvexHull3D::GenerateInitialTetrahedron(const vector<Vertex3D>& rvecVertex)
{
  const Vertex3D* pVertex3DAnchor = listHullSurface.begin()->pHalfEdge3DPointer->pVertex3DPointer;
  ExtremePoint ExtremePointData{numeric_limits<unsigned>::max(), 0.0F};

  listHullSurface.begin()->Vertex3DSurfaceNormal = SurfaceNormal(*listHullSurface.begin(), fEpsilon);
  prev(listHullSurface.end())->Vertex3DSurfaceNormal = Vertex3D{-listHullSurface.begin()->Vertex3DSurfaceNormal.fX, -listHullSurface.begin()->Vertex3DSurfaceNormal.fY, -listHullSurface.begin()->Vertex3DSurfaceNormal.fZ};

  for (int i = 0; i < rvecVertex.size(); i++) {
    float fDistance = DotProduct(Vertex3D{rvecVertex[i].fX - pVertex3DAnchor->fX, rvecVertex[i].fY - pVertex3DAnchor->fY, rvecVertex[i].fZ - pVertex3DAnchor->fZ}, listHullSurface.begin()->Vertex3DSurfaceNormal);

    if (fabs(ExtremePointData.fDistance) < fabs(fDistance)) {
      ExtremePointData.uIndex = static_cast<unsigned>(i);
      ExtremePointData.fDistance = fDistance;
    }
  }

  if (fabs(ExtremePointData.fDistance) < fEpsilon)
    return ;

  listHullVertex.push_back(rvecVertex[ExtremePointData.uIndex]);

  list<Surface3D> listAddition(0);
  HalfEdge3D* pHalfEdge3DAnchor = ExtremePointData.fDistance > 0.0 ? listHullSurface.begin()->pHalfEdge3DPointer : prev(listHullSurface.end())->pHalfEdge3DPointer;

  for (int i = 0; i < ADDITIONAL_SURFACE_COUNT; i++) {
    listAddition.push_back(Surface3D{pHalfEdge3DAnchor, Vertex3D{0.0F, 0.0F, 0.0F}, ExteriorSet{list<unsigned>(0), ExtremePoint{numeric_limits<unsigned>::max(), 0.0F}}});

    HalfEdge3D* pHalfEdge3DCurrent = pHalfEdge3DAnchor;

    pHalfEdge3DAnchor = pHalfEdge3DAnchor->pHalfEdge3DNext;
    pHalfEdge3DCurrent->pHalfEdge3DPrevious = new HalfEdge3D;
    pHalfEdge3DCurrent->pHalfEdge3DNext = new HalfEdge3D;
    pHalfEdge3DCurrent->pSurface3DPointer = &*prev(listAddition.end());
    pHalfEdge3DCurrent->pHalfEdge3DPrevious->pVertex3DPointer = &*prev(listHullVertex.end());
    pHalfEdge3DCurrent->pHalfEdge3DPrevious->pHalfEdge3DPrevious = pHalfEdge3DCurrent->pHalfEdge3DNext;
    pHalfEdge3DCurrent->pHalfEdge3DPrevious->pHalfEdge3DNext = pHalfEdge3DCurrent;
    pHalfEdge3DCurrent->pHalfEdge3DPrevious->pSurface3DPointer = &*prev(listAddition.end());
    pHalfEdge3DCurrent->pHalfEdge3DNext->pVertex3DPointer = pHalfEdge3DAnchor->pVertex3DPointer;
    pHalfEdge3DCurrent->pHalfEdge3DNext->pHalfEdge3DPrevious = pHalfEdge3DCurrent;
    pHalfEdge3DCurrent->pHalfEdge3DNext->pHalfEdge3DNext = pHalfEdge3DCurrent->pHalfEdge3DPrevious;
    pHalfEdge3DCurrent->pHalfEdge3DNext->pSurface3DPointer = &*prev(listAddition.end());

    prev(listAddition.end())->Vertex3DSurfaceNormal = SurfaceNormal(*prev(listAddition.end()), fEpsilon);
  }

  for (list<Surface3D>::iterator i = listAddition.begin(), j = prev(listAddition.end()); i != listAddition.end(); j = i++) {
    HalfEdge3D (*pHalfEdge3DFront) = i->pHalfEdge3DPointer, (*pHalfEdge3DBack) = j->pHalfEdge3DPointer;
    
    pHalfEdge3DFront->pHalfEdge3DPrevious->pHalfEdge3DPair = pHalfEdge3DBack->pHalfEdge3DNext;
    pHalfEdge3DBack->pHalfEdge3DNext->pHalfEdge3DPair = pHalfEdge3DFront->pHalfEdge3DPrevious;
  }

  listHullSurface.splice(next(listHullSurface.begin()), move(listAddition));

  if (ExtremePointData.fDistance > 0.0F)
    listHullSurface.pop_front();
  else
    listHullSurface.pop_back();

  for (int i = 0; i < rvecVertex.size(); i++) {
    list<Surface3D>::iterator iterPosition = listHullSurface.begin();
    float fDistance = 0.0F;

    do {
      pVertex3DAnchor = iterPosition->pHalfEdge3DPointer->pVertex3DPointer;
      fDistance = DotProduct(Vertex3D{rvecVertex[i].fX - pVertex3DAnchor->fX, rvecVertex[i].fY - pVertex3DAnchor->fY, rvecVertex[i].fZ - pVertex3DAnchor->fZ}, iterPosition->Vertex3DSurfaceNormal);

      if (fDistance > fEpsilon) {
        iterPosition->ExteriorSetData.listIndexSet.push_back(static_cast<unsigned>(i));

        if (iterPosition->ExteriorSetData.ExtremePointData.fDistance < fDistance) {
          iterPosition->ExteriorSetData.ExtremePointData.uIndex = static_cast<unsigned>(i);
          iterPosition->ExteriorSetData.ExtremePointData.fDistance = fDistance;
        }

        break;
      } else
        iterPosition++;
    } while (iterPosition != listHullSurface.end());
  }
}

void ConvexHull3D::FindExtremePoint(const vector<Vertex3D>& rvecVertex, unsigned* const& rpuIndexArray) const
{
  float Vertex3D::*pfVertex3DMember[ELEMENT_COUNT_3D] = {&Vertex3D::fX, &Vertex3D::fY, &Vertex3D::fZ};
  ExtremePoint ExtremePointMatrix[ELEMENT_COUNT_3D][EXTREME_ELEMENT_COUNT];

  for (int i = 0; i < ELEMENT_COUNT_3D; i++)
    for (int j = 0; j < EXTREME_ELEMENT_COUNT; j++) {
      ExtremePointMatrix[i][j].uIndex = 0U;
      ExtremePointMatrix[i][j].fDistance = *rvecVertex.begin().*pfVertex3DMember[i];
    }

  for (int i = 0; i < rvecVertex.size(); i++)
    for (int j = 0; j < ELEMENT_COUNT_3D; j++)
      if (rvecVertex[i].*pfVertex3DMember[j] < ExtremePointMatrix[j][INDICATOR_MINIMUM].fDistance) {
        ExtremePointMatrix[j][INDICATOR_MINIMUM].uIndex = static_cast<unsigned>(i);
        ExtremePointMatrix[j][INDICATOR_MINIMUM].fDistance = rvecVertex[i].*pfVertex3DMember[j];
      } else if (rvecVertex[i].*pfVertex3DMember[j] > ExtremePointMatrix[j][INDICATOR_MAXIMUM].fDistance) {
        ExtremePointMatrix[j][INDICATOR_MAXIMUM].uIndex = static_cast<unsigned>(i);
        ExtremePointMatrix[j][INDICATOR_MAXIMUM].fDistance = rvecVertex[i].*pfVertex3DMember[j];
      }

  if (ExtremePointMatrix[INDICATOR_X][INDICATOR_MAXIMUM].fDistance - ExtremePointMatrix[INDICATOR_X][INDICATOR_MINIMUM].fDistance >= ExtremePointMatrix[INDICATOR_Y][INDICATOR_MAXIMUM].fDistance - ExtremePointMatrix[INDICATOR_Y][INDICATOR_MINIMUM].fDistance)
    if (ExtremePointMatrix[INDICATOR_X][INDICATOR_MAXIMUM].fDistance - ExtremePointMatrix[INDICATOR_X][INDICATOR_MINIMUM].fDistance >= ExtremePointMatrix[INDICATOR_Z][INDICATOR_MAXIMUM].fDistance - ExtremePointMatrix[INDICATOR_Z][INDICATOR_MINIMUM].fDistance) {
      rpuIndexArray[INDICATOR_MINIMUM] = ExtremePointMatrix[INDICATOR_X][INDICATOR_MINIMUM].uIndex;
      rpuIndexArray[INDICATOR_MAXIMUM] = ExtremePointMatrix[INDICATOR_X][INDICATOR_MAXIMUM].uIndex;
    } else {
      rpuIndexArray[INDICATOR_MINIMUM] = ExtremePointMatrix[INDICATOR_Z][INDICATOR_MINIMUM].uIndex;
      rpuIndexArray[INDICATOR_MAXIMUM] = ExtremePointMatrix[INDICATOR_Z][INDICATOR_MAXIMUM].uIndex;
    }
  else if (ExtremePointMatrix[INDICATOR_Y][INDICATOR_MAXIMUM].fDistance - ExtremePointMatrix[INDICATOR_Y][INDICATOR_MINIMUM].fDistance >= ExtremePointMatrix[INDICATOR_Z][INDICATOR_MAXIMUM].fDistance - ExtremePointMatrix[INDICATOR_Z][INDICATOR_MINIMUM].fDistance) {
    rpuIndexArray[INDICATOR_MINIMUM] = ExtremePointMatrix[INDICATOR_Y][INDICATOR_MINIMUM].uIndex;
    rpuIndexArray[INDICATOR_MAXIMUM] = ExtremePointMatrix[INDICATOR_Y][INDICATOR_MAXIMUM].uIndex;
  } else {
    rpuIndexArray[INDICATOR_MINIMUM] = ExtremePointMatrix[INDICATOR_Z][INDICATOR_MINIMUM].uIndex;
    rpuIndexArray[INDICATOR_MAXIMUM] = ExtremePointMatrix[INDICATOR_Z][INDICATOR_MAXIMUM].uIndex;
  }
}

void ConvexHull3D::InitializeExteriorSet(const vector<Vertex3D>& rvecVertex, ExteriorSet& rExteriorSetDataA, ExteriorSet& rExteriorSetDataB, ExteriorSet& rExteriorSetDataC)
{
  vector<Vertex3D> vecEdgeNormal(TRIANGLE_ELEMENT_COUNT);
  int iCount = TRIANGLE_ELEMENT_COUNT - 1;

  for (list<Vertex3D>::iterator i = listHullVertex.begin(), j = prev(listHullVertex.end()); i != listHullVertex.end(); j = i++) {
    vecEdgeNormal[iCount] = CrossProduct(Vertex3D{i->fX - j->fX, i->fY - j->fY, i->fZ - j->fZ}, listHullSurface.begin()->Vertex3DSurfaceNormal);
    ++iCount %= TRIANGLE_ELEMENT_COUNT;
  }

  for (int i = 0; i < rvecVertex.size(); i++) {
    float fDistance = 0.0F;

    if (fDistance = DotProduct(Vertex3D{rvecVertex[i].fX - listHullVertex.begin()->fX, rvecVertex[i].fY - listHullVertex.begin()->fY, rvecVertex[i].fZ - listHullVertex.begin()->fZ}, *vecEdgeNormal.begin()), fDistance > fEpsilon) {
      rExteriorSetDataA.listIndexSet.push_back(static_cast<unsigned>(i));

      if (rExteriorSetDataA.ExtremePointData.fDistance < fDistance) {
        rExteriorSetDataA.ExtremePointData.uIndex = static_cast<unsigned>(i);
        rExteriorSetDataA.ExtremePointData.fDistance = fDistance;
      }
    } else if (fDistance = DotProduct(Vertex3D{rvecVertex[i].fX - next(listHullVertex.begin())->fX, rvecVertex[i].fY - next(listHullVertex.begin())->fY, rvecVertex[i].fZ - next(listHullVertex.begin())->fZ}, *next(vecEdgeNormal.begin())), fDistance > fEpsilon) {
      rExteriorSetDataB.listIndexSet.push_back(static_cast<unsigned>(i));

      if (rExteriorSetDataB.ExtremePointData.fDistance < fDistance) {
        rExteriorSetDataB.ExtremePointData.uIndex = static_cast<unsigned>(i);
        rExteriorSetDataB.ExtremePointData.fDistance = fDistance;
      }
    } else if (fDistance = DotProduct(Vertex3D{rvecVertex[i].fX - prev(listHullVertex.end())->fX, rvecVertex[i].fY - prev(listHullVertex.end())->fY, rvecVertex[i].fZ - prev(listHullVertex.end())->fZ}, *prev(vecEdgeNormal.end())), fDistance > fEpsilon) {
      rExteriorSetDataC.listIndexSet.push_back(static_cast<unsigned>(i));

      if (rExteriorSetDataC.ExtremePointData.fDistance < fDistance) {
        rExteriorSetDataC.ExtremePointData.uIndex = static_cast<unsigned>(i);
        rExteriorSetDataC.ExtremePointData.fDistance = fDistance;
      }
    } else
      ;
  }

  rExteriorSetDataA.ExtremePointData.fDistance = rExteriorSetDataB.ExtremePointData.fDistance = rExteriorSetDataC.ExtremePointData.fDistance = 0.0F;
}

void ConvexHull3D::ExpandHull(const vector<Vertex3D>& rvecVertex, list<Surface3D>::iterator& riterPosition)
{
  listHullVertex.push_back(rvecVertex[riterPosition->ExteriorSetData.ExtremePointData.uIndex]);

  list<Surface3D*> listProcessed(0);
  list<HalfEdge3D*> listHorizon(0);

  BuildHorizon(rvecVertex, listProcessed, listHorizon, riterPosition, &*prev(listHullVertex.end()));

  list<Surface3D> listAddition(0);
  list<unsigned> listComposedIndexSet(0);
  HalfEdge3D* pHalfEdge3DAnchor = nullptr;

  for (list<HalfEdge3D*>::iterator i = listHorizon.begin(); i != listHorizon.end(); i++) {
    listAddition.push_back(Surface3D{new HalfEdge3D{*(pHalfEdge3DAnchor = *i)}, Vertex3D{0.0F, 0.0F, 0.0F}, ExteriorSet{list<unsigned>(0), ExtremePoint{numeric_limits<unsigned>::max(), 0.0F}}});

    HalfEdge3D* pHalfEdge3DCurrent = prev(listAddition.end())->pHalfEdge3DPointer;

    pHalfEdge3DAnchor = pHalfEdge3DAnchor->pHalfEdge3DNext;
    pHalfEdge3DCurrent->pHalfEdge3DPrevious = new HalfEdge3D;
    pHalfEdge3DCurrent->pHalfEdge3DNext = new HalfEdge3D;
    pHalfEdge3DCurrent->pSurface3DPointer = &*prev(listAddition.end());
    pHalfEdge3DCurrent->pHalfEdge3DPrevious->pVertex3DPointer = &*prev(listHullVertex.end());
    pHalfEdge3DCurrent->pHalfEdge3DPrevious->pHalfEdge3DPrevious = pHalfEdge3DCurrent->pHalfEdge3DNext;
    pHalfEdge3DCurrent->pHalfEdge3DPrevious->pHalfEdge3DNext = pHalfEdge3DCurrent;
    pHalfEdge3DCurrent->pHalfEdge3DPrevious->pSurface3DPointer = &*prev(listAddition.end());
    pHalfEdge3DCurrent->pHalfEdge3DNext->pVertex3DPointer = pHalfEdge3DAnchor->pVertex3DPointer;
    pHalfEdge3DCurrent->pHalfEdge3DNext->pHalfEdge3DPrevious = pHalfEdge3DCurrent;
    pHalfEdge3DCurrent->pHalfEdge3DNext->pHalfEdge3DNext = pHalfEdge3DCurrent->pHalfEdge3DPrevious;
    pHalfEdge3DCurrent->pHalfEdge3DNext->pSurface3DPointer = &*prev(listAddition.end());
    pHalfEdge3DCurrent->pHalfEdge3DPair->pHalfEdge3DPair = pHalfEdge3DCurrent;

    prev(listAddition.end())->Vertex3DSurfaceNormal = SurfaceNormal(*prev(listAddition.end()), fEpsilon);
  }

  for (list<Surface3D>::iterator i = listAddition.begin(), j = prev(listAddition.end()); i != listAddition.end(); j = i++) {
    HalfEdge3D (*pHalfEdge3DFront) = i->pHalfEdge3DPointer, (*pHalfEdge3DBack) = j->pHalfEdge3DPointer;

    pHalfEdge3DFront->pHalfEdge3DPrevious->pHalfEdge3DPair = pHalfEdge3DBack->pHalfEdge3DNext;
    pHalfEdge3DBack->pHalfEdge3DNext->pHalfEdge3DPair = pHalfEdge3DFront->pHalfEdge3DPrevious;
  }

  for (list<Surface3D*>::iterator i = listProcessed.begin(); i != listProcessed.end(); i++) {
    list<Surface3D>::iterator iterMark = listHullSurface.begin();

    while (&*iterMark != *i)
      iterMark++;

    listComposedIndexSet.splice(listComposedIndexSet.end(), move(iterMark->ExteriorSetData.listIndexSet));

    pHalfEdge3DAnchor = iterMark->pHalfEdge3DPointer;

    do {
      HalfEdge3D* pHalfEdge3DCurrent = pHalfEdge3DAnchor;

      pHalfEdge3DAnchor = pHalfEdge3DAnchor->pHalfEdge3DNext;

      delete pHalfEdge3DCurrent;
    } while (pHalfEdge3DAnchor != iterMark->pHalfEdge3DPointer);

    if (riterPosition == iterMark)
      riterPosition = listHullSurface.erase(riterPosition);
    else
      listHullSurface.erase(iterMark);
  }

  for (list<unsigned>::iterator i = listComposedIndexSet.begin(); i != listComposedIndexSet.end(); i++)
    for (list<Surface3D>::iterator j = listAddition.begin(); j != listAddition.end(); j++) {
      Vertex3D* pVertex3DAnchor = &*prev(listHullVertex.end());
      float fDistance = DotProduct(Vertex3D{rvecVertex[*i].fX - pVertex3DAnchor->fX, rvecVertex[*i].fY - pVertex3DAnchor->fY, rvecVertex[*i].fZ - pVertex3DAnchor->fZ}, j->Vertex3DSurfaceNormal);

      if (fDistance > fEpsilon) {
        j->ExteriorSetData.listIndexSet.push_back(*i);

        if (j->ExteriorSetData.ExtremePointData.fDistance < fDistance) {
          j->ExteriorSetData.ExtremePointData.uIndex = *i;
          j->ExteriorSetData.ExtremePointData.fDistance = fDistance;
        }

        break;
      }
    }

  listHullSurface.splice(riterPosition--, move(listAddition));

  riterPosition++;
}

void ConvexHull3D::ExpandHull(const vector<Vertex3D>& rvecVertex, const list<Vertex3D>::iterator& riterStart, const list<Vertex3D>::iterator& riterEnd, ExteriorSet& rExteriorSetData)
{
  list<Vertex3D>::iterator iterIntermediate = listHullVertex.insert(next(riterStart), rvecVertex[rExteriorSetData.ExtremePointData.uIndex]);
  ExteriorSet ExteriorSetAddition{list<unsigned>(0), ExtremePoint{numeric_limits<unsigned>::max(), 0.0F}};
  Vertex3D Vertex3DEdgeNormalA, Vertex3DEdgeNormalB;

  Vertex3DEdgeNormalA = CrossProduct(Vertex3D{iterIntermediate->fX - riterStart->fX, iterIntermediate->fY - riterStart->fY, iterIntermediate->fZ - riterStart->fZ}, listHullSurface.begin()->Vertex3DSurfaceNormal);
  Vertex3DEdgeNormalB = CrossProduct(Vertex3D{riterEnd->fX - iterIntermediate->fX, riterEnd->fY - iterIntermediate->fY, riterEnd->fZ - iterIntermediate->fZ}, listHullSurface.begin()->Vertex3DSurfaceNormal);

  for (list<unsigned>::iterator i = rExteriorSetData.listIndexSet.begin(); i != rExteriorSetData.listIndexSet.end(); true) {
    float fDistance = 0.0F;

    if (fDistance = DotProduct(Vertex3D{rvecVertex[*i].fX - riterStart->fX, rvecVertex[*i].fY - riterStart->fY, rvecVertex[*i].fZ - riterStart->fZ}, Vertex3DEdgeNormalA), fDistance > fEpsilon) {
      if (rExteriorSetData.ExtremePointData.fDistance < fDistance) {
        rExteriorSetData.ExtremePointData.uIndex = *i;
        rExteriorSetData.ExtremePointData.fDistance = fDistance;
      }

      i++;
    } else if (fDistance = DotProduct(Vertex3D{rvecVertex[*i].fX - riterEnd->fX, rvecVertex[*i].fY - riterEnd->fY, rvecVertex[*i].fZ - riterEnd->fZ}, Vertex3DEdgeNormalB), fDistance > fEpsilon) {
      ExteriorSetAddition.listIndexSet.push_back(*i);

      if (ExteriorSetAddition.ExtremePointData.fDistance < fDistance) {
        ExteriorSetAddition.ExtremePointData.uIndex = *i;
        ExteriorSetAddition.ExtremePointData.fDistance = fDistance;
      }

      i = rExteriorSetData.listIndexSet.erase(i);
    } else
      i = rExteriorSetData.listIndexSet.erase(i);
  }

  rExteriorSetData.ExtremePointData.fDistance = ExteriorSetAddition.ExtremePointData.fDistance = 0.0F;

  if (rExteriorSetData.listIndexSet.size())
    ExpandHull(rvecVertex, riterStart, iterIntermediate, rExteriorSetData);

  if (ExteriorSetAddition.listIndexSet.size())
    ExpandHull(rvecVertex, iterIntermediate, riterEnd, ExteriorSetAddition);
}

void ConvexHull3D::BuildHorizon(const vector<Vertex3D>& rvecVertex, list<Surface3D*>& rlistProcessed, list<HalfEdge3D*>& rlistHorizon, const list<Surface3D>::iterator& riterPosition, const Vertex3D* const& rpVertex3DAnchor) const
{
  rlistProcessed.push_back(&*riterPosition);

  list<Surface3D*>::iterator iterMark = rlistProcessed.begin();
  HalfEdge3D* pHalfEdge3DCurrent = riterPosition->pHalfEdge3DPointer;

  do {
    const Vertex3D* const pVertex3DCurrent = pHalfEdge3DCurrent->pVertex3DPointer;

    if (find(rlistProcessed.begin(), rlistProcessed.end(), pHalfEdge3DCurrent->pHalfEdge3DPair->pSurface3DPointer) == rlistProcessed.end())
      if (DotProduct(Vertex3D{rpVertex3DAnchor->fX - pVertex3DCurrent->fX, rpVertex3DAnchor->fY - pVertex3DCurrent->fY, rpVertex3DAnchor->fZ - pVertex3DCurrent->fZ}, pHalfEdge3DCurrent->pHalfEdge3DPair->pSurface3DPointer->Vertex3DSurfaceNormal) > -fEpsilon) {
        iterMark = rlistProcessed.insert(next(iterMark), pHalfEdge3DCurrent->pHalfEdge3DPair->pSurface3DPointer);
        pHalfEdge3DCurrent = pHalfEdge3DCurrent->pHalfEdge3DPair->pHalfEdge3DNext;
      } else {
        rlistHorizon.push_back(pHalfEdge3DCurrent);

        pHalfEdge3DCurrent = pHalfEdge3DCurrent->pHalfEdge3DNext;
      }
    else if (*prev(iterMark) == pHalfEdge3DCurrent->pHalfEdge3DPair->pSurface3DPointer) {
      iterMark--;
      pHalfEdge3DCurrent = pHalfEdge3DCurrent->pHalfEdge3DPair->pHalfEdge3DNext;
    } else
      pHalfEdge3DCurrent = pHalfEdge3DCurrent->pHalfEdge3DNext;
  } while (rlistHorizon.empty() || rlistHorizon.front()->pVertex3DPointer != rlistHorizon.back()->pHalfEdge3DPair->pVertex3DPointer);
}

void ConvexHull3D::CopySurface(const ConvexHull3D& rConvexHull3DData)
{
  list<Surface3D>::const_iterator iterSourceSurfacePosition = rConvexHull3DData.listHullSurface.cbegin();
  list<Surface3D>::iterator iterDestinationSurfacePosition = listHullSurface.begin();

  do {
    const HalfEdge3D* pHalfEdge3DSourceAnchor = iterSourceSurfacePosition->pHalfEdge3DPointer;
    HalfEdge3D* pHalfEdge3DDestinationAnchor = iterDestinationSurfacePosition->pHalfEdge3DPointer = new HalfEdge3D{nullptr, nullptr, nullptr, nullptr, &*iterDestinationSurfacePosition};

    do {
      const Vertex3D* pVertex3DSourceAnchor = pHalfEdge3DSourceAnchor->pVertex3DPointer;
      list<Vertex3D>::const_iterator iterSourceVertexPosition = rConvexHull3DData.listHullVertex.cbegin();
      list<Vertex3D>::iterator iterDestinationVertexPosition = listHullVertex.begin();

      do {
        if (pVertex3DSourceAnchor == &*iterSourceVertexPosition)
          break;

        iterSourceVertexPosition++;
        iterDestinationVertexPosition++;
      } while (true);

      pHalfEdge3DDestinationAnchor->pVertex3DPointer = &*iterDestinationVertexPosition;

      if ((pHalfEdge3DSourceAnchor = pHalfEdge3DSourceAnchor->pHalfEdge3DNext) == iterSourceSurfacePosition->pHalfEdge3DPointer) {
        iterDestinationSurfacePosition->pHalfEdge3DPointer->pHalfEdge3DPrevious = pHalfEdge3DDestinationAnchor;
        pHalfEdge3DDestinationAnchor->pHalfEdge3DNext = iterDestinationSurfacePosition->pHalfEdge3DPointer;

        break;
      }

      pHalfEdge3DDestinationAnchor = pHalfEdge3DDestinationAnchor->pHalfEdge3DNext = new HalfEdge3D{nullptr, nullptr, pHalfEdge3DDestinationAnchor, nullptr, &*iterDestinationSurfacePosition};
    } while (true);

    iterSourceSurfacePosition++;
    iterDestinationSurfacePosition++;
  } while (iterSourceSurfacePosition != rConvexHull3DData.listHullSurface.cend());

  for (list<Surface3D>::iterator i = listHullSurface.begin(); i != listHullSurface.end(); i++) {
    HalfEdge3D* pHalfEdge3DAnchor = i->pHalfEdge3DPointer;

    do
      if (pHalfEdge3DAnchor->pHalfEdge3DPair == nullptr)
        for (list<Surface3D>::iterator j = listHullSurface.begin(); j != listHullSurface.end(); j++) {
          HalfEdge3D* pHalfEdge3DCurrent = j->pHalfEdge3DPointer;
          bool bResult = false;

          if (pHalfEdge3DAnchor == pHalfEdge3DCurrent)
            continue;

          do
            if (pHalfEdge3DAnchor->pVertex3DPointer == pHalfEdge3DCurrent->pHalfEdge3DNext->pVertex3DPointer && pHalfEdge3DAnchor->pHalfEdge3DNext->pVertex3DPointer == pHalfEdge3DCurrent->pVertex3DPointer) {
              pHalfEdge3DAnchor->pHalfEdge3DPair = pHalfEdge3DCurrent;
              pHalfEdge3DCurrent->pHalfEdge3DPair = pHalfEdge3DAnchor;
              bResult = true;

              break;
            }
          while ((pHalfEdge3DCurrent = pHalfEdge3DCurrent->pHalfEdge3DNext) != j->pHalfEdge3DPointer);

          if (bResult)
            break;
        }
    while ((pHalfEdge3DAnchor = pHalfEdge3DAnchor->pHalfEdge3DNext) != i->pHalfEdge3DPointer);
  }
}

void ConvexHull3D::MoveSurface(ConvexHull3D&& rConvexHull3DData)
{
  for (list<Surface3D>::iterator i = listHullSurface.begin(), j = rConvexHull3DData.listHullSurface.begin(); i != listHullSurface.end() && j != rConvexHull3DData.listHullSurface.end(); i++) {
    i->pHalfEdge3DPointer = move(j->pHalfEdge3DPointer);
    j->pHalfEdge3DPointer = nullptr;
  }
}

void ConvexHull3D::DeleteSurface(void)
{
  for (list<Surface3D>::iterator i = listHullSurface.begin(); i != listHullSurface.end(); i++) {
    HalfEdge3D* pHalfEdge3DAnchor = i->pHalfEdge3DPointer;

    if (pHalfEdge3DAnchor == nullptr)
      return ;

    do {
      HalfEdge3D* pHalfEdge3DCurrent = pHalfEdge3DAnchor;

      pHalfEdge3DAnchor = pHalfEdge3DAnchor->pHalfEdge3DNext;

      delete pHalfEdge3DCurrent;
    } while (pHalfEdge3DAnchor != i->pHalfEdge3DPointer);
  }
}

Mesh2D::Mesh2D(const vector<vector<Vertex2D>>& rvecVertex) : vecComposition()
{
  if (rvecVertex.empty())
    return ;

  for (int i = 0; i < rvecVertex.size() - 1; i++)
    vecComposition.push_back(ConvexHull2D(rvecVertex[i + 1]));
}

Mesh2D::Mesh2D(const Mesh2D& rMesh2DData) : vecComposition(rMesh2DData.vecComposition)
{
}

Mesh2D::Mesh2D(Mesh2D&& rMesh2DData) : vecComposition(move(rMesh2DData.vecComposition))
{
}

Mesh2D::~Mesh2D(void)
{
}

Mesh2D& Mesh2D::operator=(const Mesh2D& rMesh2DData)
{
  vecComposition = rMesh2DData.vecComposition;

  return *this;
}

Mesh2D& Mesh2D::operator=(Mesh2D&& rMesh2DData)
{
  vecComposition = move(rMesh2DData.vecComposition);

  return *this;
}

ConvexHull2D& Mesh2D::operator[](const int& riIndex)
{
  return vecComposition[riIndex];
}

const ConvexHull2D& Mesh2D::operator[](const int& riIndex) const
{
  return vecComposition[riIndex];
}

int Mesh2D::Size(void) const
{
  return vecComposition.size();
}

void Mesh2D::Push(const ConvexHull2D& rConvexHull2DData)
{
  vecComposition.push_back(rConvexHull2DData);
}

void Mesh2D::Pop(void)
{
  vecComposition.pop_back();
}

void Mesh2D::ApplyTransformation(const float (*const& rpfTransformationMatrix)[ELEMENT_COUNT_2D + 1])
{
  for (ConvexHull2D& rConvexHull2DData : vecComposition)
    rConvexHull2DData.ApplyTransformation(rpfTransformationMatrix);
}

Mesh3D::Mesh3D(const vector<vector<Vertex3D>>& rvecVertex) : vecComposition()
{
  if (rvecVertex.empty())
    return ;

  for (int i = 0; i < rvecVertex.size() - 1; i++)
    vecComposition.push_back(ConvexHull3D(rvecVertex[i + 1]));
}

Mesh3D::Mesh3D(const Mesh3D& rMesh3DData) : vecComposition(rMesh3DData.vecComposition)
{
}

Mesh3D::Mesh3D(Mesh3D&& rMesh3DData) : vecComposition(move(rMesh3DData.vecComposition))
{
}

Mesh3D::~Mesh3D(void)
{
}

Mesh3D& Mesh3D::operator=(const Mesh3D& rMesh3DData)
{
  vecComposition = rMesh3DData.vecComposition;

  return *this;
}

Mesh3D& Mesh3D::operator=(Mesh3D&& rMesh3DData)
{
  vecComposition = move(rMesh3DData.vecComposition);

  return *this;
}

ConvexHull3D& Mesh3D::operator[](const int& riIndex)
{
  return vecComposition[riIndex];
}

const ConvexHull3D& Mesh3D::operator[](const int& riIndex) const
{
  return vecComposition[riIndex];
}

int Mesh3D::Size(void) const
{
  return vecComposition.size();
}

void Mesh3D::Push(const ConvexHull3D& rConvexHull3DData)
{
  vecComposition.push_back(rConvexHull3DData);
}

void Mesh3D::Pop(void)
{
  vecComposition.pop_back();
}

void Mesh3D::ApplyTransformation(const float (*const& rpfTransformationMatrix)[ELEMENT_COUNT_3D + 1])
{
  for (ConvexHull3D& rConvexHull3DData : vecComposition)
    rConvexHull3DData.ApplyTransformation(rpfTransformationMatrix);
}