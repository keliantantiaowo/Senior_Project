#ifndef DISPLAY_RESULT
#define DISPLAY_RESULT

#include <iostream>
#include <string>
#include <vector>
#include "Definition.h"
#include "Data Structure.h"

template <typename Type> inline void ShowElement(const Type* const& rpvarArray, const std::size_t& rszSize, const std::string& rstrLine, const DisplayMode& rDisplayModeSetting)
{
  for (int i = 0; i < rszSize; i++)
    std::cout << rstrLine << '[' << i << "] = " << rpvarArray[i] << std::endl;

  switch (rDisplayModeSetting) {
  case Newline :
    std::cout << std::endl;

    break;
  case NoLineBreak :
  default : break;
  }
}

template <typename Type> inline void ShowElement(const std::vector<Type>& rvecData, const std::string& rstrLine, const DisplayMode& rDisplayModeSetting)
{
  for (int i = 0; i < rvecData.size(); i++)
    std::cout << rstrLine << '[' << i << "] = " << rvecData[i] << std::endl;

  switch (rDisplayModeSetting) {
  case Newline :
    std::cout << std::endl;

    break;
  case NoLineBreak :
  default : break;
  }
}

template <> inline void ShowElement<Vertex2D>(const std::vector<Vertex2D>& rvecData, const std::string& rstrLine, const DisplayMode& rDisplayModeSetting)
{
  for (int i = 0; i < rvecData.size(); i++)
    std::cout << rstrLine << '[' << i << "] = {" << rvecData[i].fX << ", " << rvecData[i].fY << '}' << std::endl;

  switch (rDisplayModeSetting) {
  case Newline :
    std::cout << std::endl;

    break;
  case NoLineBreak :
  default : break;
  }
}

template <> inline void ShowElement<Vertex3D>(const std::vector<Vertex3D>& rvecData, const std::string& rstrLine, const DisplayMode& rDisplayModeSetting)
{
  for (int i = 0; i < rvecData.size(); i++)
    std::cout << rstrLine << '[' << i << "] = {" << rvecData[i].fX << ", " << rvecData[i].fY << ", " << rvecData[i].fZ << '}' << std::endl;

  switch (rDisplayModeSetting) {
  case Newline :
    std::cout << std::endl;

    break;
  case NoLineBreak :
  default : break;
  }
}

template <> inline void ShowElement<WeightedNormal>(const std::vector<WeightedNormal>& rvecData, const std::string& rstrLine, const DisplayMode& rDisplayModeSetting)
{
  for (int i = 0; i < rvecData.size(); i++)
    std::cout << rstrLine << '[' << i << "] = {" << rvecData[i].uIndex << ", " << rvecData[i].fAngle << '}' << std::endl;

  switch (rDisplayModeSetting) {
  case Newline :
    std::cout << std::endl;

    break;
  case NoLineBreak :
  default : break;
  }
}

template <typename Type> inline void ShowElement(const Type* const* const& rpvarMatrix, const std::size_t& rszSize, const std::string& rstrLine, const DisplayMode& rDisplayModeSetting)
{
  for (int i = 0; i < rszSize; i++)
    for (int j = 0; j < rszSize; j++) {
      std::cout << rstrLine << '[' << i << "][" << j << "] = " << rpvarMatrix[i][j];

      if (rszSize - j - 1)
        std::cout << ", ";
      else
        std::cout << std::endl;
    }

  switch (rDisplayModeSetting) {
  case Newline :
    std::cout << std::endl;

    break;
  case NoLineBreak :
  default : break;
  }
}

template <typename Type, std::size_t szSize> inline void ShowElement(const Type (*const& rpvarMatrix)[szSize], const std::string& rstrLine, const DisplayMode& rDisplayModeSetting)
{
  for (int i = 0; i < szSize; i++)
    for (int j = 0; j < szSize; j++) {
      std::cout << rstrLine << '[' << i << "][" << j << "] = " << rpvarMatrix[i][j];

      if (szSize - j - 1)
        std::cout << ", ";
      else
        std::cout << std::endl;
    }

  switch (rDisplayModeSetting) {
  case Newline :
    std::cout << std::endl;

    break;
  case NoLineBreak :
  default : break;
  }
}

template <typename Type> inline void ShowElement(const std::vector<std::vector<Type>>& rvecData, const std::string& rstrLine, const DisplayMode& rDisplayModeSetting)
{
  for (int i = 0; i < rvecData.size(); i++)
    for (int j = 0; j < rvecData[i].size(); j++) {
      std::cout << rstrLine << '[' << i << "][" << j << "] = " << rvecData[i][j];

      if (rvecData[i].size() - j - 1)
        std::cout << ", ";
      else
        std::cout << std::endl;
    }

  switch (rDisplayModeSetting) {
  case Newline :
    std::cout << std::endl;

    break;
  case NoLineBreak :
  default : break;
  }
}

template <> inline void ShowElement<Vertex2D>(const std::vector<std::vector<Vertex2D>>& rvecData, const std::string& rstrLine, const DisplayMode& rDisplayModeSetting)
{
  for (int i = 0; i < rvecData.size(); i++)
    for (int j = 0; j < rvecData[i].size(); j++) {
      std::cout << rstrLine << '[' << i << "][" << j << "] = {" << rvecData[i][j].fX << ", " << rvecData[i][j].fY << '}';

      if (rvecData.size() - i - 1)
        std::cout << ", ";
      else
        std::cout << std::endl;
    }

  switch (rDisplayModeSetting) {
  case Newline :
    std::cout << std::endl;

    break;
  case NoLineBreak :
  default : break;
  }
}

template <> inline void ShowElement<Vertex3D>(const std::vector<std::vector<Vertex3D>>& rvecData, const std::string& rstrLine, const DisplayMode& rDisplayModeSetting)
{
  for (int i = 0; i < rvecData.size(); i++)
    for (int j = 0; j < rvecData[i].size(); j++) {
      std::cout << rstrLine << '[' << i << "][" << j << "] = {" << rvecData[i][j].fX << ", " << rvecData[i][j].fY << ", " << rvecData[i][j].fZ << '}';

      if (rvecData.size() - i - 1)
        std::cout << ", ";
      else
        std::cout << std::endl;
    }

  switch (rDisplayModeSetting) {
  case Newline :
    std::cout << std::endl;

    break;
  case NoLineBreak :
  default : break;
  }
}

template <> inline void ShowElement<WeightedNormal>(const std::vector<std::vector<WeightedNormal>>& rvecData, const std::string& rstrLine, const DisplayMode& rDisplayModeSetting)
{
  for (int i = 0; i < rvecData.size(); i++)
    for (int j = 0; j < rvecData[i].size(); j++) {
      std::cout << rstrLine << '[' << i << "][" << j << "] = {" << rvecData[i][j].uIndex << ", " << rvecData[i][j].fAngle << '}';

      if (rvecData[i].size() - j - 1)
        std::cout << ", ";
      else
        std::cout << std::endl;
    }

  switch (rDisplayModeSetting) {
  case Newline :
    std::cout << std::endl;

    break;
  case NoLineBreak :
  default : break;
  }
}

void ShowElement(const Vertex2D& rVertex2DData, const std::string& rstrLine, const DisplayMode& rDisplayModeSetting);

void ShowElement(const Vertex3D& rVertex3DData, const std::string& rstrLine, const DisplayMode& rDisplayModeSetting);

#endif