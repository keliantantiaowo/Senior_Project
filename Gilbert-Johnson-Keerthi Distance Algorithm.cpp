#include <iterator>
#include <limits>
#include <cstdlib>
#include <ctime>
#include "Geometric Operation.h"
#include "Gilbert-Johnson-Keerthi Distance Algorithm.h"

using namespace std;

bool GJKIntersection2D(const ConvexHull2D& rConvexHull2DDataA, const ConvexHull2D& rConvexHull2DDataB)
{
  vector<Vertex2D> vecSimplex(TRIANGLE_ELEMENT_COUNT, Vertex2D{0.0F, 0.0F});
  Vertex2D Vertex2DDirection{0.0F, 0.0F};

  switch (srand(static_cast<unsigned>(time(NULL))), rand() / DIRECTION_COUNT_2D) {
  case XAxisPositive :
    Vertex2DDirection.fX = 1.0F;

    break;
  case XAxisNegative :
    Vertex2DDirection.fX = -1.0F;

    break;
  case YAxisPositive :
    Vertex2DDirection.fY = 1.0F;

    break;
  case YAxisNegative :
    Vertex2DDirection.fY = -1.0F;

    break;
  default : break;
  }

  Vertex2D Vertex2DSupportPointA, Vertex2DSupportPointB;

  Vertex2DSupportPointA = SupportMapping2D(rConvexHull2DDataA, Vertex2DDirection);
  Vertex2DSupportPointB = SupportMapping2D(rConvexHull2DDataB, Vertex2D{-Vertex2DDirection.fX, -Vertex2DDirection.fY});
  *vecSimplex.begin() = Vertex2D{Vertex2DSupportPointA.fX - Vertex2DSupportPointB.fX, Vertex2DSupportPointA.fY - Vertex2DSupportPointB.fY};

  if (DotProduct(*vecSimplex.begin(), Vertex2DDirection) < 0.0F)
    return false;
  else if (vecSimplex.begin()->fX == 0.0F && vecSimplex.begin()->fY == 0.0F)
    return true;

  Vertex2DDirection = Vertex2D{-vecSimplex.begin()->fX, -vecSimplex.begin()->fY};
  Vertex2DSupportPointA = SupportMapping2D(rConvexHull2DDataA, Vertex2DDirection);
  Vertex2DSupportPointB = SupportMapping2D(rConvexHull2DDataB, Vertex2D{-Vertex2DDirection.fX, -Vertex2DDirection.fY});
  *next(vecSimplex.begin()) = Vertex2D{Vertex2DSupportPointA.fX - Vertex2DSupportPointB.fX, Vertex2DSupportPointA.fY - Vertex2DSupportPointB.fY};

  if (DotProduct(*next(vecSimplex.begin()), Vertex2DDirection) < 0.0F)
    return false;
  else if (vecSimplex.begin()->fX * next(vecSimplex.begin())->fY == vecSimplex.begin()->fY * next(vecSimplex.begin())->fX)
    return true;

  vector<Vertex2D>::iterator iterStart, iterEnd, iterMark;

  if (vecSimplex.begin()->fX * next(vecSimplex.begin())->fY - vecSimplex.begin()->fY * next(vecSimplex.begin())->fX < 0.0F) {
    iterStart = vecSimplex.begin();
    iterEnd = next(vecSimplex.begin());
    iterMark = prev(vecSimplex.end());
  } else {
    iterStart = next(vecSimplex.begin());
    iterEnd = vecSimplex.begin();
    iterMark = prev(vecSimplex.end());
  }

  while (true) {
    Vertex2DDirection = EdgeNormal(Vertex2D{iterEnd->fX - iterStart->fX, iterEnd->fY - iterStart->fY});
    Vertex2DSupportPointA = SupportMapping2D(rConvexHull2DDataA, Vertex2DDirection);
    Vertex2DSupportPointB = SupportMapping2D(rConvexHull2DDataB, Vertex2D{-Vertex2DDirection.fX, -Vertex2DDirection.fY});
    *iterMark = Vertex2D{Vertex2DSupportPointA.fX - Vertex2DSupportPointB.fX, Vertex2DSupportPointA.fY - Vertex2DSupportPointB.fY};

    if (DotProduct(*iterMark, Vertex2DDirection) < 0.0F)
      return false;
    else if (iterStart->fX * iterMark->fY - iterStart->fY * iterMark->fX < 0.0F) {
      vector<Vertex2D>::iterator iterTemp;

      iterTemp = iterEnd;
      iterEnd = iterMark;
      iterMark = iterTemp;
    } else if (iterMark->fX * iterEnd->fY - iterMark->fY * iterEnd->fX < 0.0F) {
      vector<Vertex2D>::iterator iterTemp;

      iterTemp = iterStart;
      iterStart = iterMark;
      iterMark = iterTemp;
    } else
      return true;
  }
}

bool GJKIntersection3D(const ConvexHull3D& rConvexHull3DDataA, const ConvexHull3D& rConvexHull3DDataB)
{
  vector<Vertex3D> vecSimplex(TETRAHEDRON_ELEMENT_COUNT, Vertex3D{0.0F, 0.0F, 0.0F});
  Vertex3D Vertex3DDirection{0.0F, 0.0F, 0.0F};

  switch (srand(static_cast<unsigned>(time(NULL))), rand() % DIRECTION_COUNT_3D) {
  case XAxisPositive :
    Vertex3DDirection.fX = 1.0F;

    break;
  case XAxisNegative :
    Vertex3DDirection.fX = -1.0F;

    break;
  case YAxisPositive :
    Vertex3DDirection.fY = 1.0F;

    break;
  case YAxisNegative :
    Vertex3DDirection.fY = -1.0F;

    break;
  case ZAxisPositive :
    Vertex3DDirection.fZ = 1.0F;

    break;
  case ZAxisNegative :
    Vertex3DDirection.fZ = -1.0F;

    break;
  default : break;
  }

  Vertex3D Vertex3DSupportPointA, Vertex3DSupportPointB;

  Vertex3DSupportPointA = SupportMapping3D(rConvexHull3DDataA, Vertex3DDirection);
  Vertex3DSupportPointB = SupportMapping3D(rConvexHull3DDataB, Vertex3D{-Vertex3DDirection.fX, -Vertex3DDirection.fY, -Vertex3DDirection.fZ});
  *vecSimplex.begin() = Vertex3D{Vertex3DSupportPointA.fX - Vertex3DSupportPointB.fX, Vertex3DSupportPointA.fY - Vertex3DSupportPointB.fY, Vertex3DSupportPointA.fZ - Vertex3DSupportPointB.fZ};

  if (DotProduct(*vecSimplex.begin(), Vertex3DDirection) < 0.0F)
    return false;
  else if (vecSimplex.begin()->fX == 0.0F && vecSimplex.begin()->fY == 0.0F)
    return true;

  Vertex3DDirection = Vertex3D{-vecSimplex.begin()->fX, -vecSimplex.begin()->fY, -vecSimplex.begin()->fZ};
  Vertex3DSupportPointA = SupportMapping3D(rConvexHull3DDataA, Vertex3DDirection);
  Vertex3DSupportPointB = SupportMapping3D(rConvexHull3DDataB, Vertex3D{-Vertex3DDirection.fX, -Vertex3DDirection.fY, -Vertex3DDirection.fZ});
  *next(vecSimplex.begin()) = Vertex3D{Vertex3DSupportPointA.fX - Vertex3DSupportPointB.fX, Vertex3DSupportPointA.fY - Vertex3DSupportPointB.fY, Vertex3DSupportPointA.fZ - Vertex3DSupportPointB.fZ};

  if (DotProduct(*next(vecSimplex.begin()), Vertex3DDirection) < 0.0F)
    return false;
  else if (vecSimplex.begin()->fX * next(vecSimplex.begin())->fY == vecSimplex.begin()->fY * next(vecSimplex.begin())->fX && vecSimplex.begin()->fX * next(vecSimplex.begin())->fZ == vecSimplex.begin()->fZ * next(vecSimplex.begin())->fX)
    return true;

  do {
    Vertex3D Vertex3DEdge{next(vecSimplex.begin())->fX - vecSimplex.begin()->fX, next(vecSimplex.begin())->fY - vecSimplex.begin()->fY, next(vecSimplex.begin())->fZ - vecSimplex.begin()->fZ};
    float fFactor = -DotProduct(*vecSimplex.begin(), Vertex3DEdge) / DotProduct(Vertex3DEdge, Vertex3DEdge);

    Vertex3DDirection = Vertex3D{-(vecSimplex.begin()->fX * (1.0F - fFactor) + next(vecSimplex.begin())->fX * fFactor), -(vecSimplex.begin()->fY * (1.0F - fFactor) + next(vecSimplex.begin())->fY * fFactor), -(vecSimplex.begin()->fZ * (1.0F - fFactor) + next(vecSimplex.begin())->fZ * fFactor)};
    Vertex3DSupportPointA = SupportMapping3D(rConvexHull3DDataA, Vertex3DDirection);
    Vertex3DSupportPointB = SupportMapping3D(rConvexHull3DDataB, Vertex3D{-Vertex3DDirection.fX, -Vertex3DDirection.fY, -Vertex3DDirection.fZ});
    *next(next(vecSimplex.begin())) = Vertex3D{Vertex3DSupportPointA.fX - Vertex3DSupportPointB.fX, Vertex3DSupportPointA.fY - Vertex3DSupportPointB.fY, Vertex3DSupportPointA.fZ - Vertex3DSupportPointB.fZ};

    if (DotProduct(*next(next(vecSimplex.begin())), Vertex3DDirection) < 0.0F)
      return false;
  } while (false);

  vector<Vertex3D>::iterator iterStart, iterIntermediate, iterEnd, iterMark;

  if (DotProduct(*vecSimplex.begin(), CrossProduct(*next(vecSimplex.begin()), *next(next(vecSimplex.begin())))) < 0.0F) {
    iterStart = vecSimplex.begin();
    iterIntermediate = next(vecSimplex.begin());
    iterEnd = next(next(vecSimplex.begin()));
    iterMark = prev(vecSimplex.end());
  } else {
    iterStart = vecSimplex.begin();
    iterIntermediate = next(next(vecSimplex.begin()));
    iterEnd = next(vecSimplex.begin());
    iterMark = prev(vecSimplex.end());
  }

  while (true) {
    Vertex3DDirection = CrossProduct(Vertex3D{iterEnd->fX - iterIntermediate->fX, iterEnd->fY - iterIntermediate->fY, iterEnd->fZ - iterIntermediate->fZ}, Vertex3D{iterStart->fX - iterIntermediate->fX, iterStart->fY - iterIntermediate->fY, iterStart->fZ - iterIntermediate->fZ});
    Vertex3DSupportPointA = SupportMapping3D(rConvexHull3DDataA, Vertex3DDirection);
    Vertex3DSupportPointB = SupportMapping3D(rConvexHull3DDataB, Vertex3D{-Vertex3DDirection.fX, -Vertex3DDirection.fY, -Vertex3DDirection.fZ});
    *iterMark = Vertex3D{Vertex3DSupportPointA.fX - Vertex3DSupportPointB.fX, Vertex3DSupportPointA.fY - Vertex3DSupportPointB.fY, Vertex3DSupportPointA.fZ - Vertex3DSupportPointB.fZ};

    if (DotProduct(*iterMark, Vertex3DDirection) < 0.0F)
      return false;
    else if (DotProduct(*iterStart, CrossProduct(*iterIntermediate, *iterMark)) < 0.0F) {
      vector<Vertex3D>::iterator iterTemp;

      iterTemp = iterEnd;
      iterEnd = iterMark;
      iterMark = iterTemp;
    } else if (DotProduct(*iterIntermediate, CrossProduct(*iterEnd, *iterMark)) < 0.0F) {
      vector<Vertex3D>::iterator iterTemp;

      iterTemp = iterStart;
      iterStart = iterMark;
      iterMark = iterTemp;
    } else if (DotProduct(*iterEnd, CrossProduct(*iterStart, *iterMark)) < 0.0F) {
      vector<Vertex3D>::iterator iterTemp;

      iterTemp = iterIntermediate;
      iterIntermediate = iterMark;
      iterMark = iterTemp;
    } else
      return true;
  }
}

Vertex2D SupportMapping2D(const ConvexHull2D& rConvexHull2DData, const Vertex2D& rVertex2DDirection)
{
  list<Vertex2D>::const_iterator iterPosition = rConvexHull2DData.End();
  float fExtremeValue = numeric_limits<float>::lowest();

  for (list<Vertex2D>::const_iterator i = rConvexHull2DData.Begin(); i != rConvexHull2DData.End(); i++) {
    float fDistance = DotProduct(rVertex2DDirection, *i);

    if (fExtremeValue < fDistance) {
      iterPosition = i;
      fExtremeValue = fDistance;
    }
  }

  return *iterPosition;
}

Vertex3D SupportMapping3D(const ConvexHull3D& rConvexHull3DData, const Vertex3D& rVertex3DDirection)
{
  list<Vertex3D>::const_iterator iterPosition = rConvexHull3DData.VertexEnd();
  float fExtremeValue = numeric_limits<float>::lowest();

  for (list<Vertex3D>::const_iterator i = rConvexHull3DData.VertexBegin(); i != rConvexHull3DData.VertexEnd(); i++) {
    float fDistance = DotProduct(rVertex3DDirection, *i);

    if (fExtremeValue < fDistance) {
      iterPosition = i;
      fExtremeValue = fDistance;
    }
  }

  return *iterPosition;
}