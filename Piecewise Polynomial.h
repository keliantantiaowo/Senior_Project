#ifndef PIECEWISE_POLYNOMIAL
#define PIECEWISE_POLYNOMIAL

#include <vector>
#include "Data Structure.h"

class Polynomial2D {
  float** pfCoefficient;

public:
  Polynomial2D(const float& rfInitializer = 0.0F);
  Polynomial2D(const Polynomial2D& rPolynomial2DData);
  Polynomial2D(Polynomial2D&& rPolynomial2DData);
  ~Polynomial2D(void);

  Polynomial2D& operator=(const Polynomial2D& rPolynomial2DData);
  Polynomial2D& operator=(Polynomial2D&& rPolynomial2DData);

  float*& operator[](const int& riIndex);
  const float* const& operator[](const int& riIndex) const;
};

class Polynomial3D {
  float*** pfCoefficient;

public:
  Polynomial3D(const float& rfInitializer = 0.0F);
  Polynomial3D(const Polynomial3D& rPolynomial3DData);
  Polynomial3D(Polynomial3D&& rPolynomial3DData);
  ~Polynomial3D(void);

  Polynomial3D& operator=(const Polynomial3D& rPolynomial3DData);
  Polynomial3D& operator=(Polynomial3D&& rPolynomial3DData);

  float**& operator[](const int& riIndex);
  const float* const* const& operator[](const int& riIndex) const;
};

class PiecewisePolynomial2D {
public:
  struct DivisionInfo2D {
    int iCountX;
    int iCountY;
  };

  struct CoordinateInfo2D {
    float fMinimumX;
    float fMaximumX;
    float fMinimumY;
    float fMaximumY;
  };

  struct CellInfo2D {
    float fWidthX;
    float fWidthY;
  };

private:
  Polynomial2D** pPolynomial2DMatrix;
  float** pfCell;
  DivisionInfo2D DivisionInfo2DData;
  CoordinateInfo2D CoordinateInfo2DData;
  CellInfo2D CellInfo2DData;

public:
  PiecewisePolynomial2D(const std::vector<Vertex2D>& rvecVertex, const std::vector<Vertex2D>& rvecVertexNormal, const int& riInitializerX, const int& riInitializerY, const float& rfPadding = 0.1F);
  PiecewisePolynomial2D(const std::vector<Vertex2D>& rvecVertex, const std::vector<float>& rvecCoefficient, const int& riInitializerX, const int& riInitializerY, const float& rfPadding = 0.1F);
  PiecewisePolynomial2D(const PiecewisePolynomial2D& rPiecewisePolynomialData);
  PiecewisePolynomial2D(PiecewisePolynomial2D&& rPiecewisePolynomialData);
  ~PiecewisePolynomial2D(void);

  PiecewisePolynomial2D& operator=(const PiecewisePolynomial2D& rPiecewisePolynomialData);
  PiecewisePolynomial2D& operator=(PiecewisePolynomial2D&& rPiecewisePolynomialData);

  Vertex2D MinimumCorner(void) const;
  Vertex2D MaximumCorner(void) const;

  float At(const float* const& rpfCoordinate) const;
  float At(const std::vector<float>& rvecCoordinate) const;
  float At(const Vertex2D& rVertex2DData) const;

private:
  void InitializeCell(const std::vector<Vertex2D>& rvecVertex, const std::vector<Vertex2D>& rvecVertexNormal, const float& rfPadding);
  void InitializeCell(const std::vector<Vertex2D>& rvecVertex, const std::vector<float>& rvecCoefficient, const float& rfPadding);

  void GeneratePolynomial(void);

  void ObtainRange(const std::vector<Vertex2D>& rvecVertex, const float& rfPadding);

  void CalculateCoefficient(const float* const* const& rpfSolution);

  void Clear(void);
};

class PiecewisePolynomial3D {
public:
  struct DivisionInfo3D {
    int iCountX;
    int iCountY;
    int iCountZ;
  };

  struct CoordinateInfo3D {
    float fMinimumX;
    float fMaximumX;
    float fMinimumY;
    float fMaximumY;
    float fMinimumZ;
    float fMaximumZ;
  };

  struct CellInfo3D {
    float fWidthX;
    float fWidthY;
    float fWidthZ;
  };

private:
  Polynomial3D*** pPolynomial3DMatrix;
  float*** pfCell;
  DivisionInfo3D DivisionInfo3DData;
  CoordinateInfo3D CoordinateInfo3DData;
  CellInfo3D CellInfo3DData;

public:
  PiecewisePolynomial3D(const std::vector<Vertex3D>& rvecVertex, const std::vector<Vertex3D>& rvecVertexNormal, const std::vector<Vertex3D>& rvecSurfaceNormal, const std::vector<std::vector<unsigned>>& rvecTriangleVertexIndex, const std::vector<std::vector<unsigned>>& rvecTriangleNormalIndex, const int& riInitializerX, const int& riInitializerY, const int& riInitializerZ, const float& rfPadding = 0.1F);
  PiecewisePolynomial3D(const std::vector<Vertex3D>& rvecVertex, const std::vector<float>& rvecCoefficient, const int& riInitializerX, const int& riInitializerY, const int& riInitializerZ, const float& rfPadding = 0.1F);
  PiecewisePolynomial3D(const PiecewisePolynomial3D& rPiecewisePolynomialData);
  PiecewisePolynomial3D(PiecewisePolynomial3D&& rPiecewisePolynomialData);
  ~PiecewisePolynomial3D(void);

  PiecewisePolynomial3D& operator=(const PiecewisePolynomial3D& rPiecewisePolynomialData);
  PiecewisePolynomial3D& operator=(PiecewisePolynomial3D&& rPiecewisePolynomialData);

  Vertex3D MinimumCorner(void) const;
  Vertex3D MaximumCorner(void) const;

  float At(const float* const& rpfCoordinate) const;
  float At(const std::vector<float>& rvecCoordinate) const;
  float At(const Vertex3D& rVertex3DData) const;

private:
  void InitializeCell(const std::vector<Vertex3D>& rvecVertex, const std::vector<Vertex3D>& rvecVertexNormal, const std::vector<Vertex3D>& rvecSurfaceNormal, const std::vector<std::vector<unsigned>>& rvecTriangleVertexIndex, const std::vector<std::vector<unsigned>>& rvecTriangleNormalIndex, const float& rfPadding);
  void InitializeCell(const std::vector<Vertex3D>& rvecVertex, const std::vector<float>& rvecCoefficient, const float& rfPadding);

  void GeneratePolynomial(void);

  void ObtainRange(const std::vector<Vertex3D>& rvecVertex, const float& rfPadding);

  void CalculateCoefficient(const float* const* const* const& rpfSolution);

  void Clear(void);
};

#endif