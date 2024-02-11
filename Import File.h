#ifndef IMPORT_FILE
#define IMPORT_FILE

#include <sstream>
#include <string>
#include <vector>
#include "Definition.h"
#include "Data Structure.h"

void ImportObjectFile2D(std::string strFileName, std::vector<Vertex2D>& rvecVertex);

void ImportObjectFile2D(std::string strFileName, std::vector<std::vector<Vertex2D>>& rvecVertex);

void ImportObjectFile3D(std::string strFileName, std::vector<Vertex3D>& rvecVertex, std::vector<Vertex3D>& rvecSurfaceNormal, std::vector<std::vector<unsigned>>& rvecTriangleVertexIndex, std::vector<std::vector<unsigned>>& rvecTriangleNormalIndex);

void ImportObjectFile3D(std::string strFileName, std::vector<std::vector<Vertex3D>>& rvecVertex, std::vector<std::vector<Vertex3D>>& rvecSurfaceNormal, std::vector<std::vector<std::vector<unsigned>>>& rvecTriangleVertexIndex, std::vector<std::vector<std::vector<unsigned>>>& rvecTriangleNormalIndex);

void SetValue2D(std::stringstream& rssData, Vertex2D& rVertex2DData);

void SetValue2D(std::stringstream& rssData, std::vector<Vertex2D>& rvecData);

void SetValue2D(std::stringstream& rssData, std::vector<std::vector<unsigned>>& rvecDataA, std::vector<std::vector<unsigned>>& rvecDataB);

void SetValue3D(std::stringstream& rssData, Vertex3D& rVertex3DData);

void SetValue3D(std::stringstream& rssData, std::vector<Vertex3D>& rvecData);

void SetValue3D(std::stringstream& rssData, std::vector<std::vector<unsigned>>& rvecDataA, std::vector<std::vector<unsigned>>& rvecDataB);

#endif