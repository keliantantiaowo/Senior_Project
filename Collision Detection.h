#ifndef COLLISION_DETECTION
#define COLLISION_DETECTION

#include "Definition.h"
#include "Geometric Operation.h"
#include "Piecewise Polynomial.h"
#include "Bounding Volume.h"

bool DetectCollision2D(const AABB2D& rAABB2DDataA, const AABB2D& rAABB2DDataB);

bool DetectCollision2D(const ConvexHull2D& rConvexHull2DDataA, const ConvexHull2D& rConvexHull2DDataB);

bool DetectCollision2D(const Mesh2D& rMesh2DDataA, const Mesh2D& rMesh2DDataB);

bool DetectCollision2D(const std::vector<Vertex2D>& rvecVertexA, const std::vector<Vertex2D>& rvecVertexB, const PiecewisePolynomial2D& rPiecewisePolynomial2DDataA, const PiecewisePolynomial2D& rPiecewisePolynomial2DDataB, const float (*const& rpfTransformationMatrixA)[ELEMENT_COUNT_2D + 1], const float (*const& rpfTransformationMatrixB)[ELEMENT_COUNT_2D + 1]);

bool DetectCollision3D(const AABB3D& rAABB3DDataA, const AABB3D& rAABB3DDataB);

bool DetectCollision3D(const ConvexHull3D& rConvexHull3DDataA, const ConvexHull3D& rConvexHull3DDataB);

bool DetectCollision3D(const Mesh3D& rMesh3DDataA, const Mesh3D& rMesh3DDataB);

bool DetectCollision3D(const std::vector<Vertex3D>& rvecVertexA, const std::vector<Vertex3D>& rvecVertexB, const PiecewisePolynomial3D& rPiecewisePolynomial3DDataA, const PiecewisePolynomial3D& rPiecewisePolynomial3DDataB, const float (*const& rpfTransformationMatrixA)[ELEMENT_COUNT_3D + 1], const float (*const& rpfTransformationMatrixB)[ELEMENT_COUNT_3D + 1]);

#endif