#ifndef GILBERT_JOHNSON_KEERTHI_DISTANCE_ALGORITHM
#define GILBERT_JOHNSON_KEERTHI_DISTANCE_ALGORITHM

#include <vector>
#include "Definition.h"
#include "Data Structure.h"
#include "Bounding Volume.h"

bool GJKIntersection2D(const ConvexHull2D& rConvexHull2DDataA, const ConvexHull2D& rConvexHull2DDataB);

bool GJKIntersection3D(const ConvexHull3D& rConvexHull3DDataA, const ConvexHull3D& rConvexHull3DDataB);

Vertex2D SupportMapping2D(const ConvexHull2D& rConvexHull2DData, const Vertex2D& rVertex2DDirection);

Vertex3D SupportMapping3D(const ConvexHull3D& rConvexHull3DData, const Vertex3D& rVertex3DDirection);

#endif