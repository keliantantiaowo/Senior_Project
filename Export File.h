#ifndef EXPORT_FILE
#define EXPORT_FILE

#include <string>
#include <limits>
#include "Piecewise Polynomial.h"
#include "Bounding Volume.h"

void ExportObjectFile2D(std::string strFileName, const ConvexHull2D& rConvexHull2DData);

void ExportObjectFile2D(std::string strFileName, const Mesh2D& rMesh2DData);

void ExportObjectFile2D(std::string strFileName, const std::vector<Mesh2D>& rvecPolygonSet);

void ExportObjectFile2D(std::string strFileName, const PiecewisePolynomial2D& rPiecewisePolynomial2DData, const float (*const& rpfTransformationMatrix)[ELEMENT_COUNT_2D + 1], const int& riDivisionCountX, const int& riDivisionCountY);

void ExportObjectFile2D(std::string strFileName, const std::vector<PiecewisePolynomial2D>& rvecImplicitSurface, const float (*const& rpfTransformationMatrix)[ELEMENT_COUNT_2D + 1][ELEMENT_COUNT_2D + 1], const int& riDivisionCountX, const int& riDivisionCountY);

void ExportObjectFile3D(std::string strFileName, const ConvexHull3D& rConvexHull3DData);

void ExportObjectFile3D(std::string strFileName, const Mesh3D& rMesh3DData);

void ExportObjectFile3D(std::string strFileName, const std::vector<Mesh3D>& rvecPolygonSet);

void ExportObjectFile3D(std::string strFileName, const PiecewisePolynomial3D& rPiecewisePolynomial3DData, const float (*const& rpfTransformationMatrix)[ELEMENT_COUNT_3D + 1], const int& riDivisionCountX, const int& riDivisionCountY, const int& riDivisionCountZ);

void ExportObjectFile3D(std::string strFileName, const std::vector<PiecewisePolynomial3D>& rvecImplicitSurface, const float (*const& rpfTransformationMatrix)[ELEMENT_COUNT_3D + 1][ELEMENT_COUNT_3D + 1], const int& riDivisionCountX, const int& riDivisionCountY, const int& riDivisionCountZ);

void ExportDataFile2D(std::string strFileName, const std::vector<Vertex2D>& rvecVertex, const std::vector<float>& rvecCoefficient, const float (*const& rpfTransformationMatrix)[ELEMENT_COUNT_2D + 1], const int& riDivisionCountX, const int& riDivisionCountY, const Vertex2D& rVertex2DMinimumCorner = Vertex2D{std::numeric_limits<float>::max(), std::numeric_limits<float>::max()}, const Vertex2D& rVertex2DMaximumCorner = Vertex2D{std::numeric_limits<float>::lowest(), std::numeric_limits<float>::lowest()}, const float& rfPadding = 0.1F);

void ExportDataFile2D(std::string strFileName, const PiecewisePolynomial2D& rPiecewisePolynomial2DData, const float (*const& rpfTransformationMatrix)[ELEMENT_COUNT_2D + 1], const int& riDivisionCountX, const int& riDivisionCountY, const Vertex2D& rVertex2DMinimumCorner = Vertex2D{std::numeric_limits<float>::max(), std::numeric_limits<float>::max()}, const Vertex2D& rVertex2DMaximumCorner = Vertex2D{std::numeric_limits<float>::lowest(), std::numeric_limits<float>::lowest()});

void ExportDataFile3D(std::string strFileName, const std::vector<Vertex3D>& rvecVertex, const std::vector<float>& rvecCoefficient, const float (*const& rpfTransformationMatrix)[ELEMENT_COUNT_3D + 1], const int& riDivisionCountX, const int& riDivisionCountY, const int& riDivisionCountZ, const Vertex3D& rVertex3DMinimumCorner = Vertex3D{std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max()}, const Vertex3D& rVertex3DMaximumCorner = Vertex3D{std::numeric_limits<float>::lowest(), std::numeric_limits<float>::lowest(), std::numeric_limits<float>::lowest()}, const float& rfPadding = 0.1F);

void ExportDataFile3D(std::string strFileName, const PiecewisePolynomial3D& rPiecewisePolynomial3DData, const float (*const& rpfTransformationMatrix)[ELEMENT_COUNT_3D + 1], const int& riDivisionCountX, const int& riDivisionCountY, const int& riDivisionCountZ, const Vertex3D& rVertex3DMinimumCorner = Vertex3D{std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max()}, const Vertex3D& rVertex3DMaximumCorner = Vertex3D{std::numeric_limits<float>::lowest(), std::numeric_limits<float>::lowest(), std::numeric_limits<float>::lowest()});

void ExportPartialConvexHull2D(const std::string& rstrFileName, const ConvexHull2D& rConvexHull2DData);

void ExportPartialConvexHull3D(const std::string& rstrFileName, const ConvexHull3D& rConvexHull3DData);

#endif