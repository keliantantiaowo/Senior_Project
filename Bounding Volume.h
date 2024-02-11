#ifndef BOUNDING_VOLUME
#define BOUNDING_VOLUME

#include <vector>
#include "Definition.h"
#include "Data Structure.h"

class AABB2D {
public:
  struct RangeInfo2D {
    float fWidthX;
    float fWidthY;
  };

private:
  Vertex2D Vertex2DCenter;
  RangeInfo2D RangeInfo2DData;

public:
  AABB2D(const std::vector<Vertex2D>& rvecVertex);
  AABB2D(const AABB2D& rAABB2DData);
  AABB2D(AABB2D&& rAABB2DData);
  ~AABB2D(void);

  AABB2D& operator=(const AABB2D& rAABB2DData);
  AABB2D& operator=(AABB2D&& rAABB2DData);

  Vertex2D Center(void) const;
  RangeInfo2D Size(void) const;
};

class AABB3D {
public:
  struct RangeInfo3D {
    float fWidthX;
    float fWidthY;
    float fWidthZ;
  };

private:
  Vertex3D Vertex3DCenter;
  RangeInfo3D RangeInfo3DData;

public:
  AABB3D(const std::vector<Vertex3D>& rvecVertex);
  AABB3D(const AABB3D& rAABB3DData);
  AABB3D(AABB3D&& rAABB3DData);
  ~AABB3D(void);

  AABB3D& operator=(const AABB3D& rAABB3DData);
  AABB3D& operator=(AABB3D&& rAABB3DData);

  Vertex3D Center(void) const;
  RangeInfo3D Size(void) const;
};

class ConvexHull2D {
  std::list<Vertex2D> listHullVertex;
  float fEpsilon;

public:
  ConvexHull2D(const std::vector<Vertex2D>& rvecVertex);
  ConvexHull2D(const ConvexHull2D& rConvexHull2DData);
  ConvexHull2D(ConvexHull2D&& rConvexHull2DData);
  ~ConvexHull2D(void);

  ConvexHull2D& operator=(const ConvexHull2D& rConvexHull2DData);
  ConvexHull2D& operator=(ConvexHull2D&& rConvexHull2DData);

  Vertex2D& operator[](const int& riIndex);
  const Vertex2D& operator[](const int& riIndex) const;

  std::list<Vertex2D>::iterator Begin(void);
  std::list<Vertex2D>::const_iterator Begin(void) const;
  std::list<Vertex2D>::iterator End(void);
  std::list<Vertex2D>::const_iterator End(void) const;

  int VertexCount(void) const;

  float Epsilon(void) const;

  void ApplyTransformation(const float (*const& rpfTransformationMatrix)[ELEMENT_COUNT_2D + 1]);

private:
  void GenerateInitialHull(const std::vector<Vertex2D>& rvecVertex, ExteriorSet& rExteriorSetDataA, ExteriorSet& rExteriorSetDataB);

  void FindExtremePoint(const std::vector<Vertex2D>& rvecVertex, unsigned* const& rpuIndexArray) const;

  void ExpandHull(const std::vector<Vertex2D>& rvecVertex, const std::list<Vertex2D>::iterator& riterStart, const std::list<Vertex2D>::iterator& riterEnd, ExteriorSet& rExteriorSetData);
};

class ConvexHull3D {
  std::list<Vertex3D> listHullVertex;
  std::list<Surface3D> listHullSurface;
  ShapeType ShapeTypeData;
  float fEpsilon;

public:
  ConvexHull3D(const std::vector<Vertex3D>& rvecVertex);
  ConvexHull3D(const ConvexHull3D& rConvexHull3DData);
  ConvexHull3D(ConvexHull3D&& rConvexHull3DData);
  ~ConvexHull3D(void);

  ConvexHull3D& operator=(const ConvexHull3D& rConvexHull3DData);
  ConvexHull3D& operator=(ConvexHull3D&& rConvexHull3DData);

  std::list<Vertex3D>::iterator VertexBegin(void);
  std::list<Vertex3D>::const_iterator VertexBegin(void) const;
  std::list<Vertex3D>::iterator VertexEnd(void);
  std::list<Vertex3D>::const_iterator VertexEnd(void) const;

  std::list<Surface3D>::iterator SurfaceBegin(void);
  std::list<Surface3D>::const_iterator SurfaceBegin(void) const;
  std::list<Surface3D>::iterator SurfaceEnd(void);
  std::list<Surface3D>::const_iterator SurfaceEnd(void) const;

  int VertexCount(void) const;
  int SurfaceCount(void) const;

  float Epsilon(void) const;

  const ShapeType ShapeState(void) const;

  void ApplyTransformation(const float (*const& rpfTransformationMatrix)[ELEMENT_COUNT_3D + 1]);

private:
  void GenerateInitialHull(const std::vector<Vertex3D>& rvecVertex);
  void GenerateInitialTriangle(const std::vector<Vertex3D>& rvecVertex);
  void GenerateInitialTetrahedron(const std::vector<Vertex3D>& rvecVertex);

  void FindExtremePoint(const std::vector<Vertex3D>& rvecVertex, unsigned* const& rpuIndexArray) const;

  void InitializeExteriorSet(const std::vector<Vertex3D>& rvecVertex, ExteriorSet& rExteriorSetDataA, ExteriorSet& rExteriorSetDataB, ExteriorSet& rExteriorSetDataC);

  void ExpandHull(const std::vector<Vertex3D>& rvecVertex, std::list<Surface3D>::iterator& riterPosition);
  void ExpandHull(const std::vector<Vertex3D>& rvecVertex, const std::list<Vertex3D>::iterator& riterStart, const std::list<Vertex3D>::iterator& riterEnd, ExteriorSet& rExteriorSetData);

  void BuildHorizon(const std::vector<Vertex3D>& rvecVertex, std::list<Surface3D*>& rlistProcessed, std::list<HalfEdge3D*>& rlistHorizon, const std::list<Surface3D>::iterator& riterPosition, const Vertex3D* const& rpVertex3DAnchor) const;

  void CopySurface(const ConvexHull3D& rConvexHull3DData);
  void MoveSurface(ConvexHull3D&& rConvexHull3DData);
  void DeleteSurface(void);
};

class Mesh2D {
private:
  std::vector<ConvexHull2D> vecComposition;

public:
  Mesh2D(const std::vector<std::vector<Vertex2D>>& rvecVertex);
  Mesh2D(const Mesh2D& rMesh2DData);
  Mesh2D(Mesh2D&& rMesh2DData);
  ~Mesh2D(void);

  Mesh2D& operator=(const Mesh2D& rMesh2DData);
  Mesh2D& operator=(Mesh2D&& rMesh2DData);

  ConvexHull2D& operator[](const int& riIndex);
  const ConvexHull2D& operator[](const int& riIndex) const;

  int Size(void) const;

  void Push(const ConvexHull2D& rConvexHull2DData);
  void Pop(void);

  void ApplyTransformation(const float (*const& rpfTransformationMatrix)[ELEMENT_COUNT_2D + 1]);
};

class Mesh3D {
private:
  std::vector<ConvexHull3D> vecComposition;

public:
  Mesh3D(const std::vector<std::vector<Vertex3D>>& rvecVertex);
  Mesh3D(const Mesh3D& rMesh3DData);
  Mesh3D(Mesh3D&& rMesh3DData);
  ~Mesh3D(void);

  Mesh3D& operator=(const Mesh3D& rMesh3DData);
  Mesh3D& operator=(Mesh3D&& rMesh3DData);

  ConvexHull3D& operator[](const int& riIndex);
  const ConvexHull3D& operator[](const int& riIndex) const;

  int Size(void) const;

  void Push(const ConvexHull3D& rConvexHull3DData);
  void Pop(void);

  void ApplyTransformation(const float (*const& rpfTransformationMatrix)[ELEMENT_COUNT_3D + 1]);
};

#endif