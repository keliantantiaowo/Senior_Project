#ifndef DATA_STRUCTURE
#define DATA_STRUCTURE

#include <list>

struct WeightedNormal;
struct ExtremePoint;
struct ExteriorSet;
struct Vertex2D;
struct HalfEdge2D;
struct Surface2D;
struct Vertex3D;
struct HalfEdge3D;
struct Surface3D;

enum DisplayMode {Newline, NoLineBreak};
enum NormalStyle {Normalized, Weighted};
enum Sign {Plus, Minus, None};
enum ShapeType {Plane, Solid, Undefined};
enum CoordinateAxisDirection {XAxisPositive, XAxisNegative, YAxisPositive, YAxisNegative, ZAxisPositive, ZAxisNegative};

struct WeightedNormal {
  unsigned uIndex;
  float fAngle;
};

struct ExtremePoint {
  unsigned uIndex;
  float fDistance;
};

struct ExteriorSet {
  std::list<unsigned> listIndexSet;
  ExtremePoint ExtremePointData;
};

struct Vertex2D {
  float fX;
  float fY;
  HalfEdge2D* pHalfEdge2DPointer;
};

struct HalfEdge2D {
  Vertex2D* pVertex2DPointer;
  HalfEdge2D* pHalfEdge2DPair;
  HalfEdge2D* pHalfEdge2DPrevious;
  HalfEdge2D* pHalfEdge2DNext;
  Surface2D* pSurface2DPointer;
};

struct Surface2D {
  HalfEdge2D* pHalfEdge2DPointer;
  ExteriorSet ExteriorSetData;
};

struct Vertex3D {
  float fX;
  float fY;
  float fZ;
  HalfEdge3D* pHalfEdge3DPointer;
};

struct HalfEdge3D {
  Vertex3D* pVertex3DPointer;
  HalfEdge3D* pHalfEdge3DPair;
  HalfEdge3D* pHalfEdge3DPrevious;
  HalfEdge3D* pHalfEdge3DNext;
  Surface3D* pSurface3DPointer;
};

struct Surface3D {
  HalfEdge3D* pHalfEdge3DPointer;
  Vertex3D Vertex3DSurfaceNormal;
  ExteriorSet ExteriorSetData;
};

#endif