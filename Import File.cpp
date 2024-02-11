#include <fstream>
#include <iomanip>
#include <cstdlib>
#include "Display Result.h"
#include "Import File.h"

using namespace std;

void ImportObjectFile2D(string strFileName, vector<Vertex2D>& rvecVertex)
{
  ifstream ifsFile;

  if (strFileName.empty()) {
    cout << "読み込むファイル名を入力してください：" << flush;

    getline(cin, strFileName);

    cout << endl;
  }

  ifsFile.open(strFileName + OBJECT_FILE_EXTENSION);

  if (!ifsFile.is_open()) {
    cerr << "ファイル「" << strFileName << OBJECT_FILE_EXTENSION << "」を開けませんでした。" << endl;

    exit(EXIT_FAILURE);
  }

  rvecVertex.clear();

  string strLine;

  do {
    getline(ifsFile, strLine);

    if (ifsFile.fail())
      if (ifsFile.eof())
        break;
      else {
        cerr << "ファイル「" << strFileName << OBJECT_FILE_EXTENSION << "」を読み込めませんでした。" << endl;

        exit(EXIT_FAILURE);
      }

    stringstream ssData(strLine);
    char chCharacter;

    ssData >> chCharacter;

    switch (chCharacter) {
    case VERTEX_CHARACTER :
      ssData >> chCharacter;

      switch (chCharacter) {
      case NORMAL_CHARACTER :
      case TEXTURE_CHARACTER :
      case PARAMETER_SPACE_CHARACTER : break;
      default :
        ssData.putback(chCharacter);

        SetValue2D(ssData, rvecVertex);

        break;
      }

      break;
    case FACE_ELEMENT_CHARACTER :
    default : break;
    }
  } while (true);

  ShowElement(rvecVertex, "rvecVertex", Newline);

  ifsFile.close();
}

void ImportObjectFile2D(string strFileName, vector<vector<Vertex2D>>& rvecVertex)
{
  ifstream ifsFile;

  if (strFileName.empty()) {
    cout << "読み込むファイル名を入力してください：" << flush;

    getline(cin, strFileName);

    cout << endl;
  }

  rvecVertex.assign({vector<Vertex2D>(0)});

  ImportObjectFile2D(strFileName, *rvecVertex.begin());

  int iCount = 0;

  while (true) {
    ostringstream ossData(strFileName);
    ios::fmtflags othFlag = cout.flags();
    char chCharacter = cout.fill();
    string::size_type othPosition = strFileName.rfind(DELIMITER_CHARACTER);

    ossData.seekp(0, ios::end);

    if (othPosition == string::npos)
      ossData << DELIMITER_CHARACTER << strFileName;
    else
      ossData << strFileName.substr(othPosition);

    ossData << ' ' << setw(DIGIT) << setfill(FILL_CHARACTER) << iCount + 1;
    cout << setiosflags(othFlag) << setfill(chCharacter);

    ifsFile.open(ossData.str() + OBJECT_FILE_EXTENSION);

    if (!ifsFile.is_open())
      if (iCount)
        return ;
      else {
        cerr << "ファイル「" << strFileName << OBJECT_FILE_EXTENSION << "」の凸分解データが見付かりませんでした。" << endl;

        exit(EXIT_FAILURE);
      }

    rvecVertex.push_back(vector<Vertex2D>(0));

    string strLine;

    do {
      getline(ifsFile, strLine);

      if (ifsFile.fail())
        if (ifsFile.eof())
          break;
        else {
          cerr << "ファイル「" << ossData.str() << OBJECT_FILE_EXTENSION << "を読み込めませんでした。" << endl;

          exit(EXIT_FAILURE);
        }

      stringstream ssData(strLine);

      ssData >> chCharacter;

      switch (chCharacter) {
      case VERTEX_CHARACTER :
        ssData >> chCharacter;

        switch (chCharacter) {
        case NORMAL_CHARACTER :
        case TEXTURE_CHARACTER :
        case PARAMETER_SPACE_CHARACTER : break;
        default :
          ssData.putback(chCharacter);

          SetValue2D(ssData, rvecVertex[iCount + 1]);

          break;
        }

        break;
      case FACE_ELEMENT_CHARACTER :
      default : break;
      }
    } while (true);

    ShowElement(rvecVertex[iCount + 1], "rvecVertex[" + to_string(iCount + 1) + ']', Newline);

    iCount++;

    ifsFile.close();
  }
}

void ImportObjectFile3D(string strFileName, vector<Vertex3D>& rvecVertex, vector<Vertex3D>& rvecSurfaceNormal, vector<vector<unsigned>>& rvecTriangleVertexIndex, vector<vector<unsigned>>& rvecTriangleNormalIndex)
{
  ifstream ifsFile;

  if (strFileName.empty()) {
    cout << "読み込むファイル名を入力してください：" << flush;

    getline(cin, strFileName);

    cout << endl;
  }

  ifsFile.open(strFileName + OBJECT_FILE_EXTENSION);

  if (!ifsFile.is_open()) {
    cerr << "ファイル「" << strFileName << OBJECT_FILE_EXTENSION << "」を開けませんでした。" << endl;

    exit(EXIT_FAILURE);
  }

  rvecVertex.clear();
  rvecSurfaceNormal.clear();
  rvecTriangleVertexIndex.clear();
  rvecTriangleNormalIndex.clear();

  string strLine;

  do {
    getline(ifsFile, strLine);

    if (ifsFile.fail())
      if (ifsFile.eof())
        break;
      else {
        cerr << "ファイル「" << strFileName << OBJECT_FILE_EXTENSION << "」を読み込めませんでした。" << endl;

        exit(EXIT_FAILURE);
      }

    stringstream ssData(strLine);
    char chCharacter;

    ssData >> chCharacter;

    switch (chCharacter) {
    case VERTEX_CHARACTER :
      ssData >> chCharacter;

      switch (chCharacter) {
      case NORMAL_CHARACTER :
        SetValue3D(ssData, rvecSurfaceNormal);

        break;
      case TEXTURE_CHARACTER :
      case PARAMETER_SPACE_CHARACTER : break;
      default :
        ssData.putback(chCharacter);

        SetValue3D(ssData, rvecVertex);

        break;
      }

      break;
    case FACE_ELEMENT_CHARACTER :
      SetValue3D(ssData, rvecTriangleVertexIndex, rvecTriangleNormalIndex);

      break;
    default : break;
    }
  } while (true);

  ShowElement(rvecVertex, "rvecVertex", Newline);
  ShowElement(rvecSurfaceNormal, "rvecSurfaceNormal", Newline);
  ShowElement(rvecTriangleVertexIndex, "rvecTriangleVertexIndex", Newline);
  ShowElement(rvecTriangleNormalIndex, "rvecTriangleNormalIndex", Newline);

  ifsFile.close();
}

void ImportObjectFile3D(string strFileName, vector<vector<Vertex3D>>& rvecVertex, vector<vector<Vertex3D>>& rvecSurfaceNormal, vector<vector<vector<unsigned>>>& rvecTriangleVertexIndex, vector<vector<vector<unsigned>>>& rvecTriangleNormalIndex)
{
  ifstream ifsFile;

  if (strFileName.empty()) {
    cout << "読み込むファイル名を入力してください：" << flush;

    getline(cin, strFileName);

    cout << endl;
  }

  rvecVertex.assign({vector<Vertex3D>(0)});
  rvecSurfaceNormal.assign({vector<Vertex3D>(0)});
  rvecTriangleVertexIndex.assign({vector<vector<unsigned>>(0)});
  rvecTriangleNormalIndex.assign({vector<vector<unsigned>>(0)});

  ImportObjectFile3D(strFileName, *rvecVertex.begin(), *rvecSurfaceNormal.begin(), *rvecTriangleVertexIndex.begin(), *rvecTriangleNormalIndex.begin());

  int iCount = 0;

  while (true) {
    ostringstream ossData(strFileName);
    ios::fmtflags othFlag = cout.flags();
    char chCharacter = cout.fill();
    string::size_type othPosition = strFileName.rfind(DELIMITER_CHARACTER);

    ossData.seekp(0, ios::end);

    if (othPosition == string::npos)
      ossData << DELIMITER_CHARACTER << strFileName;
    else
      ossData << strFileName.substr(othPosition);

    ossData << ' ' << setw(DIGIT) << setfill(FILL_CHARACTER) << iCount + 1;
    cout << setiosflags(othFlag) << setfill(chCharacter);

    ifsFile.open(ossData.str() + OBJECT_FILE_EXTENSION);

    if (!ifsFile.is_open())
      if (iCount)
        return ;
      else {
        cerr << "ファイル「" << strFileName << OBJECT_FILE_EXTENSION << "」の凸分解データが見付かりませんでした。" << endl;

        exit(EXIT_FAILURE);
      }

    rvecVertex.push_back(vector<Vertex3D>(0));
    rvecSurfaceNormal.push_back(vector<Vertex3D>(0));
    rvecTriangleVertexIndex.push_back(vector<vector<unsigned>>(0));
    rvecTriangleNormalIndex.push_back(vector<vector<unsigned>>(0));

    string strLine;

    do {
      getline(ifsFile, strLine);

      if (ifsFile.fail())
        if (ifsFile.eof())
          break;
        else {
          cerr << "ファイル「" << ossData.str() << OBJECT_FILE_EXTENSION << "」を読み込めませんでした。" << endl;

          exit(EXIT_FAILURE);
        }

      stringstream ssData(strLine);

      ssData >> chCharacter;

      switch (chCharacter) {
      case VERTEX_CHARACTER :
        ssData >> chCharacter;

        switch (chCharacter) {
        case NORMAL_CHARACTER :
          SetValue3D(ssData, rvecSurfaceNormal[iCount + 1]);

          break;
        case TEXTURE_CHARACTER :
        case PARAMETER_SPACE_CHARACTER : break;
        default :
          ssData.putback(chCharacter);

          SetValue3D(ssData, rvecVertex[iCount + 1]);

          break;
        }

        break;
      case FACE_ELEMENT_CHARACTER :
        SetValue3D(ssData, rvecTriangleVertexIndex[iCount + 1], rvecTriangleNormalIndex[iCount + 1]);

        break;
      default : break;
      }
    } while (true);

    ShowElement(rvecVertex[iCount + 1], "rvecVertex[" + to_string(iCount) + ']', Newline);
    ShowElement(rvecSurfaceNormal[iCount + 1], "rvecSurfaceNormal[" + to_string(iCount + 1) + ']', Newline);
    ShowElement(rvecTriangleVertexIndex[iCount + 1], "rvecTriangleVertexIndex[" + to_string(iCount + 1) + ']', Newline);
    ShowElement(rvecTriangleNormalIndex[iCount + 1], "rvecTriangleNormalIndex[" + to_string(iCount + 1) + ']', Newline);

    iCount++;

    ifsFile.close();
  }
}

void SetValue2D(stringstream& rssData, Vertex2D& rVertex2DData)
{
  rssData >> rVertex2DData.fX >> rVertex2DData.fY;

  rssData.clear();
  rssData.str("");
}

void SetValue2D(stringstream& rssData, vector<Vertex2D>& rvecData)
{
  rvecData.push_back(Vertex2D{0.0F, 0.0F});
  
  rssData >> rvecData[rvecData.size() - 1].fX >> rvecData[rvecData.size() - 1].fY;

  rssData.clear();
  rssData.str("");
}

void SetValue2D(stringstream& rssData, vector<vector<unsigned>>& rvecDataA, vector<vector<unsigned>>& rvecDataB)
{
  rvecDataA.push_back(vector<unsigned>(EDGE_ELEMENT_COUNT, 0U));
  rvecDataB.push_back(vector<unsigned>(EDGE_ELEMENT_COUNT, 0U));

  for (int i = 0; i < EDGE_ELEMENT_COUNT; i++) {
    rssData >> rvecDataA[rvecDataA.size() - 1][i];

    rvecDataA[rvecDataA.size() - 1][i]--;

    rssData.ignore(DELIMITER_COUNT);

    rssData >> rvecDataB[rvecDataB.size() - 1][i];

    rvecDataB[rvecDataB.size() - 1][i]--;
  }

  rssData.clear();
  rssData.str("");
}

void SetValue3D(stringstream& rssData, Vertex3D& rVertex3DData)
{
  rssData >> rVertex3DData.fX >> rVertex3DData.fY >> rVertex3DData.fZ;

  rssData.clear();
  rssData.str("");
}

void SetValue3D(stringstream& rssData, vector<Vertex3D>& rvecData)
{
  rvecData.push_back(Vertex3D{0.0F, 0.0F, 0.0F});
  
  rssData >> rvecData[rvecData.size() - 1].fX >> rvecData[rvecData.size() - 1].fY >> rvecData[rvecData.size() - 1].fZ;

  rssData.clear();
  rssData.str("");
}

void SetValue3D(stringstream& rssData, vector<vector<unsigned>>& rvecDataA, vector<vector<unsigned>>& rvecDataB)
{
  rvecDataA.push_back(vector<unsigned>(TRIANGLE_ELEMENT_COUNT, 0U));
  rvecDataB.push_back(vector<unsigned>(TRIANGLE_ELEMENT_COUNT, 0U));

  for (int i = 0; i < TRIANGLE_ELEMENT_COUNT; i++) {
    rssData >> rvecDataA[rvecDataA.size() - 1][i];

    rvecDataA[rvecDataA.size() - 1][i]--;

    rssData.ignore(DELIMITER_COUNT);

    rssData >> rvecDataB[rvecDataB.size() - 1][i];

    rvecDataB[rvecDataB.size() - 1][i]--;
  }

  rssData.clear();
  rssData.str("");
}