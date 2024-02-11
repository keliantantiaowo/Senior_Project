#include "Display Result.h"

using namespace std;

void ShowElement(const Vertex2D& rVertex2DData, const string& rstrLine, const DisplayMode& rDisplayModeSetting)
{
  cout << rstrLine << " = {" << rVertex2DData.fX << ", " << rVertex2DData.fY << '}' << endl;

  switch (rDisplayModeSetting) {
  case Newline :
    cout << endl;

    break;
  case NoLineBreak :
  default : break;
  }
}

void ShowElement(const Vertex3D& rVertex3DData, const string& rstrLine, const DisplayMode& rDisplayModeSetting)
{
  cout << rstrLine << " = {" << rVertex3DData.fX << ", " << rVertex3DData.fY << ", " << rVertex3DData.fZ << '}' << endl;

  switch (rDisplayModeSetting) {
  case Newline :
    cout << endl;

    break;
  case NoLineBreak :
  default : break;
  }
}