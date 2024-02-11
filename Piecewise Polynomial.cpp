#include <limits>
#include "Geometric Operation.h"
#include "Equation Solver.h"
#include "Generate Implicit Surface.h"
#include "Piecewise Polynomial.h"

using namespace std;

Polynomial2D::Polynomial2D(const float& rfInitializer) : pfCoefficient(new float*[DEGREE + 1])
{
  for (int i = 0; i <= DEGREE; i++) {
    pfCoefficient[i] = new float[DEGREE + 1];

    for (int j = 0; j <= DEGREE; j++)
      pfCoefficient[i][j] = rfInitializer;
  }
}

Polynomial2D::Polynomial2D(const Polynomial2D& rPolynomial2DData) : pfCoefficient(new float*[DEGREE + 1])
{
  for (int i = 0; i <= DEGREE; i++) {
    pfCoefficient[i] = new float[DEGREE + 1];

    for (int j = 0; j <= DEGREE; j++)
      pfCoefficient[i][j] = rPolynomial2DData.pfCoefficient[i][j];
  }
}

Polynomial2D::Polynomial2D(Polynomial2D&& rPolynomial2DData) : pfCoefficient(move(rPolynomial2DData.pfCoefficient))
{
  rPolynomial2DData.pfCoefficient = nullptr;
}

Polynomial2D::~Polynomial2D(void)
{
  if (pfCoefficient == nullptr)
    return ;

  for (int i = 0; i <= DEGREE; i++)
    delete[] pfCoefficient[i];

  delete[] pfCoefficient;
}

Polynomial2D& Polynomial2D::operator=(const Polynomial2D& rPolynomial2DData)
{
  for (int i = 0; i <= DEGREE; i++)
    for (int j = 0; j <= DEGREE; j++)
      pfCoefficient[i][j] = rPolynomial2DData.pfCoefficient[i][j];

  return *this;
}

Polynomial2D& Polynomial2D::operator=(Polynomial2D&& rPolynomial2DData)
{
  for (int i = 0; i <= DEGREE; i++)
    delete[] pfCoefficient[i];

  delete[] pfCoefficient;

  pfCoefficient = move(rPolynomial2DData.pfCoefficient);
  rPolynomial2DData.pfCoefficient = nullptr;

  return *this;
}

float*& Polynomial2D::operator[](const int& riIndex)
{
  return pfCoefficient[riIndex];
}

const float* const& Polynomial2D::operator[](const int& riIndex) const
{
  return const_cast<const float* const&>(pfCoefficient[riIndex]);
}

Polynomial3D::Polynomial3D(const float& rfInitializer) : pfCoefficient(new float**[DEGREE + 1])
{
  for (int i = 0; i <= DEGREE; i++) {
    pfCoefficient[i] = new float*[DEGREE + 1];

    for (int j = 0; j <= DEGREE; j++) {
      pfCoefficient[i][j] = new float[DEGREE + 1];

      for (int k = 0; k <= DEGREE; k++)
        pfCoefficient[i][j][k] = rfInitializer;
    }
  }
}

Polynomial3D::Polynomial3D(const Polynomial3D& rPolynomial3DData) : pfCoefficient(new float**[DEGREE + 1])
{
  for (int i = 0; i <= DEGREE; i++) {
    pfCoefficient[i] = new float*[DEGREE + 1];

    for (int j = 0; j <= DEGREE; j++) {
      pfCoefficient[i][j] = new float[DEGREE + 1];

      for (int k = 0; k <= DEGREE; k++)
        pfCoefficient[i][j][k] = rPolynomial3DData.pfCoefficient[i][j][k];
    }
  }
}

Polynomial3D::Polynomial3D(Polynomial3D&& rPolynomial3DData) : pfCoefficient(move(rPolynomial3DData.pfCoefficient))
{
  rPolynomial3DData.pfCoefficient = nullptr;
}

Polynomial3D::~Polynomial3D(void)
{
  if (pfCoefficient == nullptr)
    return ;

  for (int i = 0; i <= DEGREE; i++) {
    for (int j = 0; j <= DEGREE; j++)
      delete[] pfCoefficient[i][j];

    delete[] pfCoefficient[i];
  }

  delete[] pfCoefficient;
}

Polynomial3D& Polynomial3D::operator=(const Polynomial3D& rPolynomial3DData)
{
  for (int i = 0; i <= DEGREE; i++)
    for (int j = 0; j <= DEGREE; j++)
      for (int k = 0; k <= DEGREE; k++)
        pfCoefficient[i][j][k] = rPolynomial3DData.pfCoefficient[i][j][k];

  return *this;
}

Polynomial3D& Polynomial3D::operator=(Polynomial3D&& rPolynomial3DData)
{
  for (int i = 0; i <= DEGREE; i++) {
    for (int j = 0; j <= DEGREE; j++)
      delete[] pfCoefficient[i][j];

    delete[] pfCoefficient[i];
  }

  delete[] pfCoefficient;

  pfCoefficient = move(rPolynomial3DData.pfCoefficient);
  rPolynomial3DData.pfCoefficient = nullptr;

  return *this;
}

float**& Polynomial3D::operator[](const int& riIndex)
{
  return pfCoefficient[riIndex];
}

const float* const* const& Polynomial3D::operator[](const int& riIndex) const
{
  return const_cast<const float* const* const&>(pfCoefficient[riIndex]);
}

PiecewisePolynomial2D::PiecewisePolynomial2D(const vector<Vertex2D>& rvecVertex, const vector<Vertex2D>& rvecVertexNormal, const int& riInitializerX, const int& riInitializerY, const float& rfPadding) : pPolynomial2DMatrix(new Polynomial2D*[riInitializerX]), pfCell(new float*[riInitializerX]), DivisionInfo2DData{riInitializerX, riInitializerY}, CoordinateInfo2DData{numeric_limits<float>::max(), numeric_limits<float>::lowest(), numeric_limits<float>::max(), numeric_limits<float>::lowest()}, CellInfo2DData{0.0F, 0.0F}
{
  for (int i = 0; i < DivisionInfo2DData.iCountX; i++) {
    pPolynomial2DMatrix[i] = new Polynomial2D[DivisionInfo2DData.iCountY];
    pfCell[i] = new float[DivisionInfo2DData.iCountY];
  }

  InitializeCell(rvecVertex, rvecVertexNormal, rfPadding);

  GeneratePolynomial();
}

PiecewisePolynomial2D::PiecewisePolynomial2D(const vector<Vertex2D>& rvecVertex, const vector<float>& rvecCoefficient, const int& riInitializerX, const int& riInitializerY, const float& rfPadding) : pPolynomial2DMatrix(new Polynomial2D*[riInitializerX]), pfCell(new float*[riInitializerX]), DivisionInfo2DData{riInitializerX, riInitializerY}, CoordinateInfo2DData{numeric_limits<float>::max(), numeric_limits<float>::lowest(), numeric_limits<float>::max(), numeric_limits<float>::lowest()}, CellInfo2DData{0.0F, 0.0F}
{
  for (int i = 0; i < DivisionInfo2DData.iCountX; i++) {
    pPolynomial2DMatrix[i] = new Polynomial2D[DivisionInfo2DData.iCountY];
    pfCell[i] = new float[riInitializerY];
  }

  InitializeCell(rvecVertex, rvecCoefficient, rfPadding);

  GeneratePolynomial();
}

PiecewisePolynomial2D::PiecewisePolynomial2D(const PiecewisePolynomial2D& rPiecewisePolynomial2DData) : pPolynomial2DMatrix(new Polynomial2D*[rPiecewisePolynomial2DData.DivisionInfo2DData.iCountX]), pfCell(new float*[rPiecewisePolynomial2DData.DivisionInfo2DData.iCountX]), DivisionInfo2DData{rPiecewisePolynomial2DData.DivisionInfo2DData}, CoordinateInfo2DData{rPiecewisePolynomial2DData.CoordinateInfo2DData}, CellInfo2DData{rPiecewisePolynomial2DData.CellInfo2DData}
{
  for (int i = 0; i < DivisionInfo2DData.iCountX; i++) {
    pPolynomial2DMatrix[i] = new Polynomial2D[DivisionInfo2DData.iCountY];
    pfCell[i] = new float[DivisionInfo2DData.iCountY];

    for (int j = 0; j < DivisionInfo2DData.iCountY; j++) {
      pPolynomial2DMatrix[i][j] = rPiecewisePolynomial2DData.pPolynomial2DMatrix[i][j];
      pfCell[i][j] = rPiecewisePolynomial2DData.pfCell[i][j];
    }
  }
}

PiecewisePolynomial2D::PiecewisePolynomial2D(PiecewisePolynomial2D&& rPiecewisePolynomial2DData) : pPolynomial2DMatrix(move(rPiecewisePolynomial2DData.pPolynomial2DMatrix)), pfCell(move(rPiecewisePolynomial2DData.pfCell)), DivisionInfo2DData{move(rPiecewisePolynomial2DData.DivisionInfo2DData)}, CoordinateInfo2DData{move(rPiecewisePolynomial2DData.CoordinateInfo2DData)}, CellInfo2DData{move(rPiecewisePolynomial2DData.CellInfo2DData)}
{
  rPiecewisePolynomial2DData.pPolynomial2DMatrix = nullptr;
  rPiecewisePolynomial2DData.pfCell = nullptr;
}

PiecewisePolynomial2D::~PiecewisePolynomial2D(void)
{
  Clear();
}

PiecewisePolynomial2D& PiecewisePolynomial2D::operator=(const PiecewisePolynomial2D& rPiecewisePolynomial2DData)
{
  if (DivisionInfo2DData.iCountX == rPiecewisePolynomial2DData.DivisionInfo2DData.iCountX && DivisionInfo2DData.iCountY == rPiecewisePolynomial2DData.DivisionInfo2DData.iCountY) {
    for (int i = 0; i < DivisionInfo2DData.iCountX; i++)
      for (int j = 0; j < DivisionInfo2DData.iCountY; j++) {
        pPolynomial2DMatrix[i][j] = rPiecewisePolynomial2DData.pPolynomial2DMatrix[i][j];
        pfCell[i][j] = rPiecewisePolynomial2DData.pfCell[i][j];
      }

    CoordinateInfo2DData = rPiecewisePolynomial2DData.CoordinateInfo2DData;
    CellInfo2DData = rPiecewisePolynomial2DData.CellInfo2DData;
  } else
    throw length_error("代入式両辺の区分的多項式のサイズが異なります。");

  return *this;
}

PiecewisePolynomial2D& PiecewisePolynomial2D::operator=(PiecewisePolynomial2D&& rPiecewisePolynomial2DData)
{
  if (DivisionInfo2DData.iCountX == rPiecewisePolynomial2DData.DivisionInfo2DData.iCountX && DivisionInfo2DData.iCountY == rPiecewisePolynomial2DData.DivisionInfo2DData.iCountY) {
    Clear();

    pPolynomial2DMatrix = move(rPiecewisePolynomial2DData.pPolynomial2DMatrix);
    pfCell = move(rPiecewisePolynomial2DData.pfCell);
    CoordinateInfo2DData = move(rPiecewisePolynomial2DData.CoordinateInfo2DData);
    CellInfo2DData = move(rPiecewisePolynomial2DData.CellInfo2DData);
    rPiecewisePolynomial2DData.pPolynomial2DMatrix = nullptr;
    rPiecewisePolynomial2DData.pfCell = nullptr;
  } else
    throw length_error("代入式両辺の区分的多項式のサイズが異なります。");

  return *this;
}

Vertex2D PiecewisePolynomial2D::MinimumCorner(void) const
{
  return Vertex2D{CoordinateInfo2DData.fMinimumX, CoordinateInfo2DData.fMinimumY};
}

Vertex2D PiecewisePolynomial2D::MaximumCorner(void) const
{
  return Vertex2D{CoordinateInfo2DData.fMaximumX, CoordinateInfo2DData.fMaximumY};
}

float PiecewisePolynomial2D::At(const float* const& rpfCoordinate) const
{
  float fTransformedCoordinate[ELEMENT_COUNT_2D] = {(rpfCoordinate[INDICATOR_X] - CoordinateInfo2DData.fMinimumX) / CellInfo2DData.fWidthX, (rpfCoordinate[INDICATOR_Y] - CoordinateInfo2DData.fMinimumY) / CellInfo2DData.fWidthY};
  int iIndex[ELEMENT_COUNT_2D] = {static_cast<int>(floorf(fTransformedCoordinate[INDICATOR_X])), static_cast<int>(floorf(fTransformedCoordinate[INDICATOR_Y]))};

  if (signbit(iIndex[INDICATOR_X])) {
    fTransformedCoordinate[INDICATOR_X] = 0.0F;
    iIndex[INDICATOR_X] = 0;
  } else if (iIndex[INDICATOR_X] >= DivisionInfo2DData.iCountX) {
    fTransformedCoordinate[INDICATOR_X] = 1.0F;
    iIndex[INDICATOR_X] = DivisionInfo2DData.iCountX - 1;
  } else
    fTransformedCoordinate[INDICATOR_X] = modff(fTransformedCoordinate[INDICATOR_X], &fTransformedCoordinate[INDICATOR_X]);

  if (signbit(iIndex[INDICATOR_Y])) {
    fTransformedCoordinate[INDICATOR_Y] = 0.0F;
    iIndex[INDICATOR_Y] = 0;
  } else if (iIndex[INDICATOR_Y] >= DivisionInfo2DData.iCountY) {
    fTransformedCoordinate[INDICATOR_Y] = 1.0F;
    iIndex[INDICATOR_Y] = DivisionInfo2DData.iCountY - 1;
  } else
    fTransformedCoordinate[INDICATOR_Y] = modff(fTransformedCoordinate[INDICATOR_Y], &fTransformedCoordinate[INDICATOR_Y]);

  float fResult = 0.0F;

  for (int i = 0; i <= DEGREE; i++)
    for (int j = 0; j <= DEGREE; j++)
      fResult += pow(fTransformedCoordinate[INDICATOR_X], static_cast<float>(i)) * pow(fTransformedCoordinate[INDICATOR_Y], static_cast<float>(j)) * pPolynomial2DMatrix[iIndex[INDICATOR_X]][iIndex[INDICATOR_Y]][i][j];

  return fResult;
}

float PiecewisePolynomial2D::At(const vector<float>& rvecCoordinate) const
{
  vector<float> vecTransformedCoordinate{(rvecCoordinate[INDICATOR_X] - CoordinateInfo2DData.fMinimumX) / CellInfo2DData.fWidthX, (rvecCoordinate[INDICATOR_Y] - CoordinateInfo2DData.fMinimumY) / CellInfo2DData.fWidthY};
  int iIndex[ELEMENT_COUNT_2D] = {static_cast<int>(floorf(vecTransformedCoordinate[INDICATOR_X])), static_cast<int>(floorf(vecTransformedCoordinate[INDICATOR_Y]))};

  if (signbit(iIndex[INDICATOR_X])) {
    vecTransformedCoordinate[INDICATOR_X] = 0.0F;
    iIndex[INDICATOR_X] = 0;
  } else if (iIndex[INDICATOR_X] >= DivisionInfo2DData.iCountX) {
    vecTransformedCoordinate[INDICATOR_X] = 1.0F;
    iIndex[INDICATOR_X] = DivisionInfo2DData.iCountX - 1;
  } else
    vecTransformedCoordinate[INDICATOR_X] = modff(vecTransformedCoordinate[INDICATOR_X], &vecTransformedCoordinate[INDICATOR_X]);

  if (signbit(iIndex[INDICATOR_Y])) {
    vecTransformedCoordinate[INDICATOR_Y] = 0.0F;
    iIndex[INDICATOR_Y] = 0;
  } else if (iIndex[INDICATOR_Y] >= DivisionInfo2DData.iCountY) {
    vecTransformedCoordinate[INDICATOR_Y] = 1.0F;
    iIndex[INDICATOR_Y] = DivisionInfo2DData.iCountY - 1;
  } else
    vecTransformedCoordinate[INDICATOR_Y] = modff(vecTransformedCoordinate[INDICATOR_Y], &vecTransformedCoordinate[INDICATOR_Y]);

  float fResult = 0.0F;

  for (int i = 0; i <= DEGREE; i++)
    for (int j = 0; j <= DEGREE; j++)
      fResult += pow(vecTransformedCoordinate[INDICATOR_X], static_cast<float>(i)) * pow(vecTransformedCoordinate[INDICATOR_Y], static_cast<float>(j)) * pPolynomial2DMatrix[iIndex[INDICATOR_X]][iIndex[INDICATOR_Y]][i][j];

  return fResult;
}

float PiecewisePolynomial2D::At(const Vertex2D& rVertex2DData) const
{
  Vertex2D Vertex2DTransformedCoordinate{(rVertex2DData.fX - CoordinateInfo2DData.fMinimumX) / CellInfo2DData.fWidthX, (rVertex2DData.fY - CoordinateInfo2DData.fMinimumY) / CellInfo2DData.fWidthY};
  int iIndex[ELEMENT_COUNT_2D] = {static_cast<int>(floorf(Vertex2DTransformedCoordinate.fX)), static_cast<int>(floorf(Vertex2DTransformedCoordinate.fY))};

  if (signbit(iIndex[INDICATOR_X])) {
    Vertex2DTransformedCoordinate.fX = 0.0F;
    iIndex[INDICATOR_X] = 0;
  } else if (iIndex[INDICATOR_X] >= DivisionInfo2DData.iCountX) {
    Vertex2DTransformedCoordinate.fX = 1.0F;
    iIndex[INDICATOR_X] = DivisionInfo2DData.iCountX - 1;
  } else
    Vertex2DTransformedCoordinate.fX = modff(Vertex2DTransformedCoordinate.fX, &Vertex2DTransformedCoordinate.fX);

  if (signbit(iIndex[INDICATOR_Y])) {
    Vertex2DTransformedCoordinate.fY = 0.0F;
    iIndex[INDICATOR_Y] = 0;
  } else if (iIndex[INDICATOR_Y] >= DivisionInfo2DData.iCountY) {
    Vertex2DTransformedCoordinate.fY = 1.0F;
    iIndex[INDICATOR_Y] = DivisionInfo2DData.iCountY - 1;
  } else
    Vertex2DTransformedCoordinate.fY = modff(Vertex2DTransformedCoordinate.fY, &Vertex2DTransformedCoordinate.fY);

  float fResult = 0.0F;

  for (int i = 0; i <= DEGREE; i++)
    for (int j = 0; j <= DEGREE; j++)
      fResult += pow(Vertex2DTransformedCoordinate.fX, static_cast<float>(i)) * pow(Vertex2DTransformedCoordinate.fY, static_cast<float>(j)) * pPolynomial2DMatrix[iIndex[INDICATOR_X]][iIndex[INDICATOR_Y]][i][j];

  return fResult;
}

void PiecewisePolynomial2D::InitializeCell(const vector<Vertex2D>& rvecVertex, const vector<Vertex2D>& rvecVertexNormal, const float& rfPadding)
{
  ObtainRange(rvecVertex, rfPadding);

  for (int i = 0; i < DivisionInfo2DData.iCountX; i++)
    for (int j = 0; j < DivisionInfo2DData.iCountY; j++) {
      Vertex2D Vertex2DCenter{CoordinateInfo2DData.fMinimumX + CellInfo2DData.fWidthX * (static_cast<float>(i) + 0.5F), CoordinateInfo2DData.fMinimumY + CellInfo2DData.fWidthY * (static_cast<float>(j) + 0.5F)};
      float fExtremeValue = numeric_limits<float>::max();
      int iIndex = 0;

      do {
        float fDistance = Distance(rvecVertex, rvecVertexNormal, Vertex2DCenter, iIndex++);

        if (fabs(fExtremeValue) > fabs(fDistance))
          fExtremeValue = -fDistance;
      } while (rvecVertex.size() > iIndex);

      pfCell[i][j] = fExtremeValue;
    }
}

void PiecewisePolynomial2D::InitializeCell(const vector<Vertex2D>& rvecVertex, const vector<float>& rvecCoefficient, const float& rfPadding)
{
  ObtainRange(rvecVertex, rfPadding);

  for (int i = 0; i < DivisionInfo2DData.iCountX; i++)
    for (int j = 0; j < DivisionInfo2DData.iCountY; j++) {
      float fCellCenter[ELEMENT_COUNT_2D] = {CoordinateInfo2DData.fMinimumX + CellInfo2DData.fWidthX * (static_cast<float>(i) + 0.5F), CoordinateInfo2DData.fMinimumY + CellInfo2DData.fWidthY * (static_cast<float>(j) + 0.5F)};

      pfCell[i][j] = FunctionValue2D(fCellCenter, rvecVertex, rvecCoefficient);
    }
}

void PiecewisePolynomial2D::ObtainRange(const vector<Vertex2D>& rvecVertex, const float& rfPadding)
{
  for (int i = 0; i < rvecVertex.size(); i++) {
    if (rvecVertex[i].fX < CoordinateInfo2DData.fMinimumX)
      CoordinateInfo2DData.fMinimumX = rvecVertex[i].fX;

    if (rvecVertex[i].fX > CoordinateInfo2DData.fMaximumX)
      CoordinateInfo2DData.fMaximumX = rvecVertex[i].fX;

    if (rvecVertex[i].fY < CoordinateInfo2DData.fMinimumY)
      CoordinateInfo2DData.fMinimumY = rvecVertex[i].fY;

    if (rvecVertex[i].fY > CoordinateInfo2DData.fMaximumY)
      CoordinateInfo2DData.fMaximumY = rvecVertex[i].fY;
  }

  float fDifference;
  
  fDifference = CoordinateInfo2DData.fMaximumX - CoordinateInfo2DData.fMinimumX;
  CoordinateInfo2DData.fMinimumX -= rfPadding * fDifference;
  CoordinateInfo2DData.fMaximumX += rfPadding * fDifference;
  fDifference = CoordinateInfo2DData.fMaximumY - CoordinateInfo2DData.fMinimumY;
  CoordinateInfo2DData.fMinimumY -= rfPadding * fDifference;
  CoordinateInfo2DData.fMaximumY += rfPadding * fDifference;
  CellInfo2DData.fWidthX = (CoordinateInfo2DData.fMaximumX - CoordinateInfo2DData.fMinimumX) / static_cast<float>(DivisionInfo2DData.iCountX);
  CellInfo2DData.fWidthY = (CoordinateInfo2DData.fMaximumY - CoordinateInfo2DData.fMinimumY) / static_cast<float>(DivisionInfo2DData.iCountY);
}

void PiecewisePolynomial2D::GeneratePolynomial(void)
{
  float** pfSolution = new float*[DivisionInfo2DData.iCountX + BORDER_SIZE];
  float** pfMatrixX = new float*[DivisionInfo2DData.iCountX];
  float** pfMatrixY = new float*[DivisionInfo2DData.iCountY];
  vector<float> vecArrayX(DivisionInfo2DData.iCountX, 0.0F), vecArrayY(DivisionInfo2DData.iCountY, 0.0F);

  for (int i = 0; i < DivisionInfo2DData.iCountX + BORDER_SIZE; i++)
    pfSolution[i] = new float[DivisionInfo2DData.iCountY + BORDER_SIZE];

  for (int i = 0; i < DivisionInfo2DData.iCountX; i++) {
    pfMatrixX[i] = new float[DivisionInfo2DData.iCountX]();

    if (!i) {
      pfMatrixX[i][i] = 0.875F;
      pfMatrixX[i][i + 1] = 0.125F;
    } else if (!(DivisionInfo2DData.iCountX - i - 1)) {
      pfMatrixX[i][i] = 0.875F;
      pfMatrixX[i][i - 1] = 0.125F;
    } else {
      pfMatrixX[i][i] = 0.75F;
      pfMatrixX[i][i + 1] = pfMatrixX[i][i - 1] = 0.125F;
    }
  }

  for (int i = 0; i < DivisionInfo2DData.iCountY; i++) {
    pfMatrixY[i] = new float[DivisionInfo2DData.iCountY]();

    if (!i) {
      pfMatrixY[i][i] = 0.875F;
      pfMatrixY[i][i + 1] = 0.125F;
    } else if (!(DivisionInfo2DData.iCountY - i - 1)) {
      pfMatrixY[i][i] = 0.875F;
      pfMatrixY[i][i - 1] = 0.125F;
    } else {
      pfMatrixY[i][i] = 0.75F;
      pfMatrixY[i][i + 1] = pfMatrixY[i][i - 1] = 0.125F;
    }
  }

  for (int i = 0; i < DivisionInfo2DData.iCountX; i++) {
    for (int j = 0; j < DivisionInfo2DData.iCountY; j++)
      vecArrayY[j] = pfCell[i][j];

    SolveEquation(pfMatrixY, vecArrayY, DivisionInfo2DData.iCountY);

    for (int j = 0; j < DivisionInfo2DData.iCountY; j++)
      pfSolution[i + 1][j + 1] = vecArrayY[j];
  }

  for (int i = 0; i < DivisionInfo2DData.iCountY; i++) {
    for (int j = 0; j < DivisionInfo2DData.iCountX; j++)
      vecArrayX[j] = pfSolution[j + 1][i + 1];

    SolveEquation(pfMatrixX, vecArrayX, DivisionInfo2DData.iCountX);

    for (int j = 0; j < DivisionInfo2DData.iCountX; j++)
      pfSolution[j + 1][i + 1] = vecArrayX[j];
  }

  pfSolution[0][0] = pfSolution[1][1];
  pfSolution[DivisionInfo2DData.iCountX + 1][0] = pfSolution[DivisionInfo2DData.iCountX][1];
  pfSolution[0][DivisionInfo2DData.iCountY + 1] = pfSolution[1][DivisionInfo2DData.iCountY];
  pfSolution[DivisionInfo2DData.iCountX + 1][DivisionInfo2DData.iCountY + 1] = pfSolution[DivisionInfo2DData.iCountX][DivisionInfo2DData.iCountY];

  for (int i = 0; i < DivisionInfo2DData.iCountX; i++) {
    pfSolution[i + 1][0] = pfSolution[i + 1][1];
    pfSolution[i + 1][DivisionInfo2DData.iCountY + 1] = pfSolution[i + 1][DivisionInfo2DData.iCountY];
  }

  for (int i = 0; i < DivisionInfo2DData.iCountY; i++) {
    pfSolution[0][i + 1] = pfSolution[1][i + 1];
    pfSolution[DivisionInfo2DData.iCountX + 1][i + 1] = pfSolution[DivisionInfo2DData.iCountX][i + 1];
  }

  CalculateCoefficient(pfSolution);

  for (int i = 0; i < DivisionInfo2DData.iCountX + BORDER_SIZE; i++)
    delete[] pfSolution[i];

  for (int i = 0; i < DivisionInfo2DData.iCountX; i++)
    delete[] pfMatrixX[i];

  for (int i = 0; i < DivisionInfo2DData.iCountY; i++)
    delete[] pfMatrixY[i];

  delete[] pfSolution;
  delete[] pfMatrixX;
  delete[] pfMatrixY;
}

void PiecewisePolynomial2D::CalculateCoefficient(const float* const* const& rpfSolution)
{
  float fSubCoefficient[DEGREE + 1][DEGREE + 1] = {
    {0.5F, 0.5F, 0.0F},
    {-1.0F, 1.0F, 0.0F},
    {0.5F, -1.0F, 0.5F}
  };

  for (int i = 0; i < DivisionInfo2DData.iCountX; i++)
    for (int j = 0; j < DivisionInfo2DData.iCountY; j++)
      for (int k = 0; k <= DEGREE; k++)
        for (int l = 0; l <= DEGREE; l++)
          for (int m = 0; m <= DEGREE; m++)
            for (int n = 0; n <= DEGREE; n++)
              pPolynomial2DMatrix[i][j][k][l] += rpfSolution[i + m][j + n] * fSubCoefficient[k][m] * fSubCoefficient[l][n];
}

void PiecewisePolynomial2D::Clear(void)
{
  if (pPolynomial2DMatrix == nullptr && pfCell == nullptr)
    return ;

  for (int i = 0; i < DivisionInfo2DData.iCountX; i++) {
    delete[] pPolynomial2DMatrix[i];
    delete[] pfCell[i];
  }

  delete[] pPolynomial2DMatrix;
  delete[] pfCell;
}

PiecewisePolynomial3D::PiecewisePolynomial3D(const vector<Vertex3D>& rvecVertex, const vector<Vertex3D>& rvecVertexNormal, const vector<Vertex3D>& rvecSurfaceNormal, const vector<vector<unsigned>>& rvecTriangleVertexIndex, const vector<vector<unsigned>>& rvecTriangleNormalIndex, const int& riInitializerX, const int& riInitializerY, const int& riInitializerZ, const float& rfPadding) : pPolynomial3DMatrix(new Polynomial3D**[riInitializerX]), pfCell(new float**[riInitializerX]), DivisionInfo3DData{riInitializerX, riInitializerY, riInitializerZ}, CoordinateInfo3DData{numeric_limits<float>::max(), numeric_limits<float>::lowest(), numeric_limits<float>::max(), numeric_limits<float>::lowest(), numeric_limits<float>::max(), numeric_limits<float>::lowest()}, CellInfo3DData{0.0F, 0.0F, 0.0F}
{
  for (int i = 0; i < DivisionInfo3DData.iCountX; i++) {
    pPolynomial3DMatrix[i] = new Polynomial3D*[DivisionInfo3DData.iCountY];
    pfCell[i] = new float*[DivisionInfo3DData.iCountY];

    for (int j = 0; j < DivisionInfo3DData.iCountY; j++) {
      pPolynomial3DMatrix[i][j] = new Polynomial3D[DivisionInfo3DData.iCountZ];
      pfCell[i][j] = new float[DivisionInfo3DData.iCountZ];
    }
  }

  InitializeCell(rvecVertex, rvecVertexNormal, rvecSurfaceNormal, rvecTriangleVertexIndex, rvecTriangleNormalIndex, rfPadding);

  GeneratePolynomial();
}

PiecewisePolynomial3D::PiecewisePolynomial3D(const vector<Vertex3D>& rvecVertex, const vector<float>& rvecCoefficient, const int& riInitializerX, const int& riInitializerY, const int& riInitializerZ, const float& rfPadding) : pPolynomial3DMatrix(new Polynomial3D**[riInitializerX]), pfCell(new float**[riInitializerX]), DivisionInfo3DData{riInitializerX, riInitializerY, riInitializerZ}, CoordinateInfo3DData{numeric_limits<float>::max(), numeric_limits<float>::lowest(), numeric_limits<float>::max(), numeric_limits<float>::lowest(), numeric_limits<float>::max(), numeric_limits<float>::lowest()}, CellInfo3DData{0.0F, 0.0F, 0.0F}
{
  for (int i = 0; i < DivisionInfo3DData.iCountX; i++) {
    pPolynomial3DMatrix[i] = new Polynomial3D*[DivisionInfo3DData.iCountY];
    pfCell[i] = new float*[DivisionInfo3DData.iCountY];

    for (int j = 0; j < DivisionInfo3DData.iCountY; j++) {
      pPolynomial3DMatrix[i][j] = new Polynomial3D[DivisionInfo3DData.iCountZ];
      pfCell[i][j] = new float[DivisionInfo3DData.iCountZ];
    }
  }

  InitializeCell(rvecVertex, rvecCoefficient, rfPadding);

  GeneratePolynomial();
}

PiecewisePolynomial3D::PiecewisePolynomial3D(const PiecewisePolynomial3D& rPiecewisePolynomial3DData) : pPolynomial3DMatrix(new Polynomial3D**[rPiecewisePolynomial3DData.DivisionInfo3DData.iCountX]), pfCell(new float**[rPiecewisePolynomial3DData.DivisionInfo3DData.iCountX]), DivisionInfo3DData{rPiecewisePolynomial3DData.DivisionInfo3DData}, CoordinateInfo3DData{rPiecewisePolynomial3DData.CoordinateInfo3DData}, CellInfo3DData{rPiecewisePolynomial3DData.CellInfo3DData}
{
  for (int i = 0; i < DivisionInfo3DData.iCountX; i++) {
    pPolynomial3DMatrix[i] = new Polynomial3D*[DivisionInfo3DData.iCountY];
    pfCell[i] = new float*[DivisionInfo3DData.iCountY];

    for (int j = 0; j < DivisionInfo3DData.iCountY; j++) {
      pPolynomial3DMatrix[i][j] = new Polynomial3D[DivisionInfo3DData.iCountZ];
      pfCell[i][j] = new float[DivisionInfo3DData.iCountZ];

      for (int k = 0; k < DivisionInfo3DData.iCountZ; k++) {
        pPolynomial3DMatrix[i][j][k] = rPiecewisePolynomial3DData.pPolynomial3DMatrix[i][j][k];
        pfCell[i][j][k] = rPiecewisePolynomial3DData.pfCell[i][j][k];
      }
    }
  }
}

PiecewisePolynomial3D::PiecewisePolynomial3D(PiecewisePolynomial3D&& rPiecewisePolynomial3DData) : pPolynomial3DMatrix(move(rPiecewisePolynomial3DData.pPolynomial3DMatrix)), pfCell(move(rPiecewisePolynomial3DData.pfCell)), DivisionInfo3DData{move(rPiecewisePolynomial3DData.DivisionInfo3DData)}, CoordinateInfo3DData{move(rPiecewisePolynomial3DData.CoordinateInfo3DData)}, CellInfo3DData{move(rPiecewisePolynomial3DData.CellInfo3DData)}
{
  rPiecewisePolynomial3DData.pPolynomial3DMatrix = nullptr;
  rPiecewisePolynomial3DData.pfCell = nullptr;
}

PiecewisePolynomial3D::~PiecewisePolynomial3D(void)
{
  Clear();
}

PiecewisePolynomial3D& PiecewisePolynomial3D::operator=(const PiecewisePolynomial3D& rPiecewisePolynomial3DData)
{
  if (DivisionInfo3DData.iCountX == rPiecewisePolynomial3DData.DivisionInfo3DData.iCountX && DivisionInfo3DData.iCountY == rPiecewisePolynomial3DData.DivisionInfo3DData.iCountY && DivisionInfo3DData.iCountZ == rPiecewisePolynomial3DData.DivisionInfo3DData.iCountZ) {
    for (int i = 0; i < DivisionInfo3DData.iCountX; i++)
      for (int j = 0; j < DivisionInfo3DData.iCountY; j++)
        for (int k = 0; k < DivisionInfo3DData.iCountZ; k++) {
          pPolynomial3DMatrix[i][j][k] = rPiecewisePolynomial3DData.pPolynomial3DMatrix[i][j][k];
          pfCell[i][j][k] = rPiecewisePolynomial3DData.pfCell[i][j][k];
        }

    CoordinateInfo3DData = rPiecewisePolynomial3DData.CoordinateInfo3DData;
    CellInfo3DData = rPiecewisePolynomial3DData.CellInfo3DData;
  } else
    throw length_error("代入式両辺の区分的多項式のサイズが異なります。");

  return *this;
}

PiecewisePolynomial3D& PiecewisePolynomial3D::operator=(PiecewisePolynomial3D&& rPiecewisePolynomial3DData)
{
  if (DivisionInfo3DData.iCountX == rPiecewisePolynomial3DData.DivisionInfo3DData.iCountX && DivisionInfo3DData.iCountY == rPiecewisePolynomial3DData.DivisionInfo3DData.iCountY && DivisionInfo3DData.iCountZ == rPiecewisePolynomial3DData.DivisionInfo3DData.iCountZ) {
    Clear();

    pPolynomial3DMatrix = move(rPiecewisePolynomial3DData.pPolynomial3DMatrix);
    pfCell = move(rPiecewisePolynomial3DData.pfCell);
    CoordinateInfo3DData = move(rPiecewisePolynomial3DData.CoordinateInfo3DData);
    CellInfo3DData = move(rPiecewisePolynomial3DData.CellInfo3DData);
    rPiecewisePolynomial3DData.pPolynomial3DMatrix = nullptr;
    rPiecewisePolynomial3DData.pfCell = nullptr;
  } else
    throw length_error("代入式両辺の区分的多項式のサイズが異なります。");

  return *this;
}

Vertex3D PiecewisePolynomial3D::MinimumCorner(void) const
{
  return Vertex3D{CoordinateInfo3DData.fMinimumX, CoordinateInfo3DData.fMinimumY, CoordinateInfo3DData.fMinimumZ};
}

Vertex3D PiecewisePolynomial3D::MaximumCorner(void) const
{
  return Vertex3D{CoordinateInfo3DData.fMaximumX, CoordinateInfo3DData.fMaximumY, CoordinateInfo3DData.fMaximumZ};
}

float PiecewisePolynomial3D::At(const float* const& rpfCoordinate) const
{
  float fTransformedCoordinate[ELEMENT_COUNT_3D] = {(rpfCoordinate[INDICATOR_X] - CoordinateInfo3DData.fMinimumX) / CellInfo3DData.fWidthX, (rpfCoordinate[INDICATOR_Y] - CoordinateInfo3DData.fMinimumY) / CellInfo3DData.fWidthY, (rpfCoordinate[INDICATOR_Z] - CoordinateInfo3DData.fMinimumZ) / CellInfo3DData.fWidthZ};
  int iIndex[ELEMENT_COUNT_3D] = {static_cast<int>(floorf(fTransformedCoordinate[INDICATOR_X])), static_cast<int>(floorf(fTransformedCoordinate[INDICATOR_Y])), static_cast<int>(floorf(fTransformedCoordinate[INDICATOR_Z]))};

  if (signbit(iIndex[INDICATOR_X])) {
    fTransformedCoordinate[INDICATOR_X] = 0.0F;
    iIndex[INDICATOR_X] = 0;
  } else if (iIndex[INDICATOR_X] >= DivisionInfo3DData.iCountX) {
    fTransformedCoordinate[INDICATOR_X] = 1.0F;
    iIndex[INDICATOR_X] = DivisionInfo3DData.iCountX - 1;
  } else
    fTransformedCoordinate[INDICATOR_X] = modff(fTransformedCoordinate[INDICATOR_X], &fTransformedCoordinate[INDICATOR_X]);

  if (signbit(iIndex[INDICATOR_Y])) {
    fTransformedCoordinate[INDICATOR_Y] = 0.0F;
    iIndex[INDICATOR_Y] = 0;
  } else if (iIndex[INDICATOR_Y] >= DivisionInfo3DData.iCountY) {
    fTransformedCoordinate[INDICATOR_Y] = 1.0F;
    iIndex[INDICATOR_Y] = DivisionInfo3DData.iCountY - 1;
  } else
    fTransformedCoordinate[INDICATOR_Y] = modff(fTransformedCoordinate[INDICATOR_Y], &fTransformedCoordinate[INDICATOR_Y]);

  if (signbit(iIndex[INDICATOR_Z])) {
    fTransformedCoordinate[INDICATOR_Z] = 0.0F;
    iIndex[INDICATOR_Z] = 0;
  } else if (iIndex[INDICATOR_Z] >= DivisionInfo3DData.iCountZ) {
    fTransformedCoordinate[INDICATOR_Z] = 1.0F;
    iIndex[INDICATOR_Z] = DivisionInfo3DData.iCountZ - 1;
  } else
    fTransformedCoordinate[INDICATOR_Z] = modff(fTransformedCoordinate[INDICATOR_Z], &fTransformedCoordinate[INDICATOR_Z]);

  float fResult = 0.0F;

  for (int i = 0; i <= DEGREE; i++)
    for (int j = 0; j <= DEGREE; j++)
      for (int k = 0; k <= DEGREE; k++)
        fResult += pow(fTransformedCoordinate[INDICATOR_X], static_cast<float>(i)) * pow(fTransformedCoordinate[INDICATOR_Y], static_cast<float>(j)) * pow(fTransformedCoordinate[INDICATOR_Z], static_cast<float>(k)) * pPolynomial3DMatrix[iIndex[INDICATOR_X]][iIndex[INDICATOR_Y]][iIndex[INDICATOR_Z]][i][j][k];

  return fResult;
}

float PiecewisePolynomial3D::At(const vector<float>& rvecCoordinate) const
{
  vector<float> vecTransformedCoordinate{(rvecCoordinate[INDICATOR_X] - CoordinateInfo3DData.fMinimumX) / CellInfo3DData.fWidthX, (rvecCoordinate[INDICATOR_Y] - CoordinateInfo3DData.fMinimumY) / CellInfo3DData.fWidthY, (rvecCoordinate[INDICATOR_Z] - CoordinateInfo3DData.fMinimumZ) / CellInfo3DData.fWidthZ};
  int iIndex[ELEMENT_COUNT_3D] = {static_cast<int>(floorf(vecTransformedCoordinate[INDICATOR_X])), static_cast<int>(floorf(vecTransformedCoordinate[INDICATOR_Y])), static_cast<int>(floorf(vecTransformedCoordinate[INDICATOR_Z]))};

  if (signbit(iIndex[INDICATOR_X])) {
    vecTransformedCoordinate[INDICATOR_X] = 0.0F;
    iIndex[INDICATOR_X] = 0;
  } else if (iIndex[INDICATOR_X] >= DivisionInfo3DData.iCountX) {
    vecTransformedCoordinate[INDICATOR_X] = 1.0F;
    iIndex[INDICATOR_X] = DivisionInfo3DData.iCountX - 1;
  } else
    vecTransformedCoordinate[INDICATOR_X] = modff(vecTransformedCoordinate[INDICATOR_X], &vecTransformedCoordinate[INDICATOR_X]);

  if (signbit(iIndex[INDICATOR_Y])) {
    vecTransformedCoordinate[INDICATOR_Y] = 0.0F;
    iIndex[INDICATOR_Y] = 0;
  } else if (iIndex[INDICATOR_Y] >= DivisionInfo3DData.iCountY) {
    vecTransformedCoordinate[INDICATOR_Y] = 1.0F;
    iIndex[INDICATOR_Y] = DivisionInfo3DData.iCountY - 1;
  } else
    vecTransformedCoordinate[INDICATOR_Y] = modff(vecTransformedCoordinate[INDICATOR_Y], &vecTransformedCoordinate[INDICATOR_Y]);

  if (signbit(iIndex[INDICATOR_Z])) {
    vecTransformedCoordinate[INDICATOR_Z] = 0.0F;
    iIndex[INDICATOR_Z] = 0;
  } else if (iIndex[INDICATOR_Z] >= DivisionInfo3DData.iCountZ) {
    vecTransformedCoordinate[INDICATOR_Z] = 1.0F;
    iIndex[INDICATOR_Z] = DivisionInfo3DData.iCountZ - 1;
  } else
    vecTransformedCoordinate[INDICATOR_Z] = modff(vecTransformedCoordinate[INDICATOR_Z], &vecTransformedCoordinate[INDICATOR_Z]);

  float fResult = 0.0F;

  for (int i = 0; i <= DEGREE; i++)
    for (int j = 0; j <= DEGREE; j++)
      for (int k = 0; k <= DEGREE; k++)
        fResult += pow(vecTransformedCoordinate[INDICATOR_X], static_cast<float>(i)) * pow(vecTransformedCoordinate[INDICATOR_Y], static_cast<float>(j)) * pow(vecTransformedCoordinate[INDICATOR_Z], static_cast<float>(k)) * pPolynomial3DMatrix[iIndex[INDICATOR_X]][iIndex[INDICATOR_Y]][iIndex[INDICATOR_Z]][i][j][k];

  return fResult;
}

float PiecewisePolynomial3D::At(const Vertex3D& rVertex3DData) const
{
  Vertex3D Vertex3DTransformedCoordinate{(rVertex3DData.fX - CoordinateInfo3DData.fMinimumX) / CellInfo3DData.fWidthX, (rVertex3DData.fY - CoordinateInfo3DData.fMinimumY) / CellInfo3DData.fWidthY, (rVertex3DData.fZ - CoordinateInfo3DData.fMinimumZ) / CellInfo3DData.fWidthZ};
  int iIndex[ELEMENT_COUNT_3D] = {static_cast<int>(floorf(Vertex3DTransformedCoordinate.fX)), static_cast<int>(floorf(Vertex3DTransformedCoordinate.fY)), static_cast<int>(floorf(Vertex3DTransformedCoordinate.fZ))};

  if (signbit(iIndex[INDICATOR_X])) {
    Vertex3DTransformedCoordinate.fX = 0.0F;
    iIndex[INDICATOR_X] = 0;
  } else if (iIndex[INDICATOR_X] >= DivisionInfo3DData.iCountX) {
    Vertex3DTransformedCoordinate.fX = 1.0F;
    iIndex[INDICATOR_X] = DivisionInfo3DData.iCountX - 1;
  } else
    Vertex3DTransformedCoordinate.fX = modff(Vertex3DTransformedCoordinate.fX, &Vertex3DTransformedCoordinate.fX);

  if (signbit(iIndex[INDICATOR_Y])) {
    Vertex3DTransformedCoordinate.fY = 0.0F;
    iIndex[INDICATOR_Y] = 0;
  } else if (iIndex[INDICATOR_Y] >= DivisionInfo3DData.iCountY) {
    Vertex3DTransformedCoordinate.fY = 1.0F;
    iIndex[INDICATOR_Y] = DivisionInfo3DData.iCountY - 1;
  } else
    Vertex3DTransformedCoordinate.fY = modff(Vertex3DTransformedCoordinate.fY, &Vertex3DTransformedCoordinate.fY);

  if (signbit(iIndex[INDICATOR_Z])) {
    Vertex3DTransformedCoordinate.fZ = 0.0F;
    iIndex[INDICATOR_Z] = 0;
  } else if (iIndex[INDICATOR_Z] >= DivisionInfo3DData.iCountZ) {
    Vertex3DTransformedCoordinate.fZ = 1.0F;
    iIndex[INDICATOR_Z] = DivisionInfo3DData.iCountZ - 1;
  } else
    Vertex3DTransformedCoordinate.fZ = modff(Vertex3DTransformedCoordinate.fZ, &Vertex3DTransformedCoordinate.fZ);

  float fResult = 0.0F;

  for (int i = 0; i <= DEGREE; i++)
    for (int j = 0; j <= DEGREE; j++)
      for (int k = 0; k <= DEGREE; k++)
        fResult += pow(Vertex3DTransformedCoordinate.fX, static_cast<float>(i)) * pow(Vertex3DTransformedCoordinate.fY, static_cast<float>(j)) * pow(Vertex3DTransformedCoordinate.fZ, static_cast<float>(k)) * pPolynomial3DMatrix[iIndex[INDICATOR_X]][iIndex[INDICATOR_Y]][iIndex[INDICATOR_Z]][i][j][k];

  return fResult;
}

void PiecewisePolynomial3D::InitializeCell(const vector<Vertex3D>& rvecVertex, const vector<Vertex3D>& rvecVertexNormal, const vector<Vertex3D>& rvecSurfaceNormal, const vector<vector<unsigned>>& rvecTriangleVertexIndex, const vector<vector<unsigned>>& rvecTriangleNormalIndex, const float& rfPadding)
{
  ObtainRange(rvecVertex, rfPadding);

  for (int i = 0; i < DivisionInfo3DData.iCountX; i++)
    for (int j = 0; j < DivisionInfo3DData.iCountY; j++)
      for (int k = 0; k < DivisionInfo3DData.iCountZ; k++) {
        Vertex3D Vertex3DCenter{CoordinateInfo3DData.fMinimumX + CellInfo3DData.fWidthX * (static_cast<float>(i) + 0.5F), CoordinateInfo3DData.fMinimumY + CellInfo3DData.fWidthY * (static_cast<float>(j) + 0.5F), CoordinateInfo3DData.fMinimumZ + CellInfo3DData.fWidthZ * (static_cast<float>(k) + 0.5F)};
        float fExtremeValue = numeric_limits<float>::max();
        int iIndex = 0;

        do {
          float fDistance = Distance(rvecVertex, rvecVertexNormal, rvecSurfaceNormal, rvecTriangleVertexIndex, rvecTriangleNormalIndex, Vertex3DCenter, iIndex++);

          if (fabs(fExtremeValue) > fabs(fDistance))
            fExtremeValue = -fDistance;
        } while (rvecTriangleVertexIndex.size() > iIndex);

        pfCell[i][j][k] = fExtremeValue;
      }
}

void PiecewisePolynomial3D::InitializeCell(const vector<Vertex3D>& rvecVertex, const vector<float>& rvecCoefficient, const float& rfPadding)
{
  ObtainRange(rvecVertex, rfPadding);

  for (int i = 0; i < DivisionInfo3DData.iCountX; i++)
    for (int j = 0; j < DivisionInfo3DData.iCountY; j++)
      for (int k = 0; k < DivisionInfo3DData.iCountZ; k++) {
        float fCellCenter[ELEMENT_COUNT_3D] = {CoordinateInfo3DData.fMinimumX + CellInfo3DData.fWidthX * (static_cast<float>(i) + 0.5F), CoordinateInfo3DData.fMinimumY + CellInfo3DData.fWidthY * (static_cast<float>(j) + 0.5F), CoordinateInfo3DData.fMinimumZ + CellInfo3DData.fWidthZ * (static_cast<float>(k) + 0.5F)};

        pfCell[i][j][k] = FunctionValue3D(fCellCenter, rvecVertex, rvecCoefficient);
      }
}

void PiecewisePolynomial3D::GeneratePolynomial(void)
{
  float*** pfSolution = new float**[DivisionInfo3DData.iCountX + BORDER_SIZE];
  float** pfMatrixX = new float*[DivisionInfo3DData.iCountX];
  float** pfMatrixY = new float*[DivisionInfo3DData.iCountY];
  float** pfMatrixZ = new float*[DivisionInfo3DData.iCountZ];
  vector<float> vecArrayX(DivisionInfo3DData.iCountX, 0.0F), vecArrayY(DivisionInfo3DData.iCountY, 0.0F), vecArrayZ(DivisionInfo3DData.iCountZ, 0.0F);

  for (int i = 0; i < DivisionInfo3DData.iCountX + BORDER_SIZE; i++) {
    pfSolution[i] = new float*[DivisionInfo3DData.iCountY + BORDER_SIZE];

    for (int j = 0; j < DivisionInfo3DData.iCountY + BORDER_SIZE; j++)
      pfSolution[i][j] = new float[DivisionInfo3DData.iCountZ + BORDER_SIZE];
  }

  for (int i = 0; i < DivisionInfo3DData.iCountX; i++) {
    pfMatrixX[i] = new float[DivisionInfo3DData.iCountX]();

    if (!i) {
      pfMatrixX[i][i] = 0.875F;
      pfMatrixX[i][i + 1] = 0.125F;
    } else if (!(DivisionInfo3DData.iCountX - i - 1)) {
      pfMatrixX[i][i] = 0.875F;
      pfMatrixX[i][i - 1] = 0.125F;
    } else {
      pfMatrixX[i][i] = 0.75F;
      pfMatrixX[i][i + 1] = pfMatrixX[i][i - 1] = 0.125F;
    }
  }

  for (int i = 0; i < DivisionInfo3DData.iCountY; i++) {
    pfMatrixY[i] = new float[DivisionInfo3DData.iCountY]();

    if (!i) {
      pfMatrixY[i][i] = 0.875F;
      pfMatrixY[i][i + 1] = 0.125F;
    } else if (!(DivisionInfo3DData.iCountY - i - 1)) {
      pfMatrixY[i][i] = 0.875F;
      pfMatrixY[i][i - 1] = 0.125F;
    } else {
      pfMatrixY[i][i] = 0.75F;
      pfMatrixY[i][i + 1] = pfMatrixY[i][i - 1] = 0.125F;
    }
  }

  for (int i = 0; i < DivisionInfo3DData.iCountZ; i++) {
    pfMatrixZ[i] = new float[DivisionInfo3DData.iCountZ]();

    if (!i) {
      pfMatrixZ[i][i] = 0.875F;
      pfMatrixZ[i][i + 1] = 0.125F;
    } else if (!(DivisionInfo3DData.iCountZ - i - 1)) {
      pfMatrixZ[i][i] = 0.875F;
      pfMatrixZ[i][i - 1] = 0.125F;
    } else {
      pfMatrixZ[i][i] = 0.75F;
      pfMatrixZ[i][i + 1] = pfMatrixZ[i][i - 1] = 0.125F;
    }
  }

  for (int i = 0; i < DivisionInfo3DData.iCountX; i++)
    for (int j = 0; j < DivisionInfo3DData.iCountY; j++) {
      for (int k = 0; k < DivisionInfo3DData.iCountZ; k++)
        vecArrayZ[k] = pfCell[i][j][k];

      SolveEquation(pfMatrixZ, vecArrayZ, DivisionInfo3DData.iCountZ);

      for (int k = 0; k < DivisionInfo3DData.iCountZ; k++)
        pfSolution[i + 1][j + 1][k + 1] = vecArrayZ[k];
    }

  for (int i = 0; i < DivisionInfo3DData.iCountY; i++)
    for (int j = 0; j < DivisionInfo3DData.iCountZ; j++) {
      for (int k = 0; k < DivisionInfo3DData.iCountX; k++)
        vecArrayX[k] = pfSolution[k + 1][i + 1][j + 1];

      SolveEquation(pfMatrixX, vecArrayX, DivisionInfo3DData.iCountX);

      for (int k = 0; k < DivisionInfo3DData.iCountX; k++)
        pfSolution[k + 1][i + 1][j + 1] = vecArrayX[k];
    }

  for (int i = 0; i < DivisionInfo3DData.iCountZ; i++)
    for (int j = 0; j < DivisionInfo3DData.iCountX; j++) {
      for (int k = 0; k < DivisionInfo3DData.iCountY; k++)
        vecArrayY[k] = pfSolution[j + 1][k + 1][i + 1];

      SolveEquation(pfMatrixY, vecArrayY, DivisionInfo3DData.iCountY);

      for (int k = 0; k < DivisionInfo3DData.iCountY; k++)
        pfSolution[j + 1][k + 1][i + 1] = vecArrayY[k];
    }

  pfSolution[0][0][0] = pfSolution[1][1][1];
  pfSolution[DivisionInfo3DData.iCountX + 1][0][0] = pfSolution[DivisionInfo3DData.iCountX][1][1];
  pfSolution[0][DivisionInfo3DData.iCountY + 1][0] = pfSolution[1][DivisionInfo3DData.iCountY][1];
  pfSolution[0][0][DivisionInfo3DData.iCountZ + 1] = pfSolution[1][1][DivisionInfo3DData.iCountZ];
  pfSolution[DivisionInfo3DData.iCountX + 1][DivisionInfo3DData.iCountY + 1][0] = pfSolution[DivisionInfo3DData.iCountX][DivisionInfo3DData.iCountY][1];
  pfSolution[DivisionInfo3DData.iCountX + 1][0][DivisionInfo3DData.iCountZ + 1] = pfSolution[DivisionInfo3DData.iCountX][1][DivisionInfo3DData.iCountZ];
  pfSolution[0][DivisionInfo3DData.iCountY + 1][DivisionInfo3DData.iCountZ + 1] = pfSolution[1][DivisionInfo3DData.iCountY][DivisionInfo3DData.iCountZ];
  pfSolution[DivisionInfo3DData.iCountX + 1][DivisionInfo3DData.iCountY + 1][DivisionInfo3DData.iCountZ + 1] = pfSolution[DivisionInfo3DData.iCountX][DivisionInfo3DData.iCountY][DivisionInfo3DData.iCountZ];

  for (int i = 0; i < DivisionInfo3DData.iCountX; i++) {
    pfSolution[i + 1][0][0] = pfSolution[i + 1][1][1];
    pfSolution[i + 1][DivisionInfo3DData.iCountY + 1][0] = pfSolution[i + 1][DivisionInfo3DData.iCountY][1];
    pfSolution[i + 1][0][DivisionInfo3DData.iCountZ + 1] = pfSolution[i + 1][1][DivisionInfo3DData.iCountZ];
    pfSolution[i + 1][DivisionInfo3DData.iCountY + 1][DivisionInfo3DData.iCountZ + 1] = pfSolution[i + 1][DivisionInfo3DData.iCountY][DivisionInfo3DData.iCountZ];

    for (int j = 0; j < DivisionInfo3DData.iCountY; j++) {
      pfSolution[i + 1][j + 1][0] = pfSolution[i + 1][j + 1][1];
      pfSolution[i + 1][j + 1][DivisionInfo3DData.iCountZ + 1] = pfSolution[i + 1][j + 1][DivisionInfo3DData.iCountZ];
    }
  }

  for (int i = 0; i < DivisionInfo3DData.iCountY; i++) {
    pfSolution[0][i + 1][0] = pfSolution[1][i + 1][1];
    pfSolution[DivisionInfo3DData.iCountX + 1][i + 1][0] = pfSolution[DivisionInfo3DData.iCountX][i + 1][1];
    pfSolution[0][i + 1][DivisionInfo3DData.iCountZ + 1] = pfSolution[1][i + 1][DivisionInfo3DData.iCountZ];
    pfSolution[DivisionInfo3DData.iCountX + 1][i + 1][DivisionInfo3DData.iCountZ + 1] = pfSolution[DivisionInfo3DData.iCountX][i + 1][DivisionInfo3DData.iCountZ];

    for (int j = 0; j < DivisionInfo3DData.iCountZ; j++) {
      pfSolution[0][i + 1][j + 1] = pfSolution[1][i + 1][j + 1];
      pfSolution[DivisionInfo3DData.iCountX + 1][i + 1][j + 1] = pfSolution[DivisionInfo3DData.iCountX][i + 1][j + 1];
    }
  }

  for (int i = 0; i < DivisionInfo3DData.iCountZ; i++) {
    pfSolution[0][0][i + 1] = pfSolution[1][1][i + 1];
    pfSolution[DivisionInfo3DData.iCountX + 1][0][i + 1] = pfSolution[DivisionInfo3DData.iCountX][1][i + 1];
    pfSolution[0][DivisionInfo3DData.iCountY + 1][i + 1] = pfSolution[1][DivisionInfo3DData.iCountY][i + 1];
    pfSolution[DivisionInfo3DData.iCountX + 1][DivisionInfo3DData.iCountY + 1][i + 1] = pfSolution[DivisionInfo3DData.iCountX][DivisionInfo3DData.iCountY][i + 1];

    for (int j = 0; j < DivisionInfo3DData.iCountX; j++) {
      pfSolution[j + 1][0][i + 1] = pfSolution[j + 1][1][i + 1];
      pfSolution[j + 1][DivisionInfo3DData.iCountY + 1][i + 1] = pfSolution[j + 1][DivisionInfo3DData.iCountY][i + 1];
    }
  }

  CalculateCoefficient(pfSolution);

  for (int i = 0; i < DivisionInfo3DData.iCountX + BORDER_SIZE; i++) {
    for (int j = 0; j < DivisionInfo3DData.iCountY + BORDER_SIZE; j++)
      delete[] pfSolution[i][j];

    delete[] pfSolution[i];
  }

  for (int i = 0; i < DivisionInfo3DData.iCountX; i++)
    delete[] pfMatrixX[i];

  for (int i = 0; i < DivisionInfo3DData.iCountY; i++)
    delete[] pfMatrixY[i];

  for (int i = 0; i < DivisionInfo3DData.iCountZ; i++)
    delete[] pfMatrixZ[i];

  delete[] pfSolution;
  delete[] pfMatrixX;
  delete[] pfMatrixY;
  delete[] pfMatrixZ;
}

void PiecewisePolynomial3D::ObtainRange(const vector<Vertex3D>& rvecVertex, const float& rfPadding)
{
  for (int i = 0; i < rvecVertex.size(); i++) {
    if (rvecVertex[i].fX < CoordinateInfo3DData.fMinimumX)
      CoordinateInfo3DData.fMinimumX = rvecVertex[i].fX;

    if (rvecVertex[i].fX > CoordinateInfo3DData.fMaximumX)
      CoordinateInfo3DData.fMaximumX = rvecVertex[i].fX;

    if (rvecVertex[i].fY < CoordinateInfo3DData.fMinimumY)
      CoordinateInfo3DData.fMinimumY = rvecVertex[i].fY;

    if (rvecVertex[i].fY > CoordinateInfo3DData.fMaximumY)
      CoordinateInfo3DData.fMaximumY = rvecVertex[i].fY;

    if (rvecVertex[i].fZ < CoordinateInfo3DData.fMinimumZ)
      CoordinateInfo3DData.fMinimumZ = rvecVertex[i].fZ;

    if (rvecVertex[i].fZ > CoordinateInfo3DData.fMaximumZ)
      CoordinateInfo3DData.fMaximumZ = rvecVertex[i].fZ;
  }

  float fDifference;

  fDifference = CoordinateInfo3DData.fMaximumX - CoordinateInfo3DData.fMinimumX;
  CoordinateInfo3DData.fMinimumX -= rfPadding * fDifference;
  CoordinateInfo3DData.fMaximumX += rfPadding * fDifference;
  fDifference = CoordinateInfo3DData.fMaximumY - CoordinateInfo3DData.fMinimumY;
  CoordinateInfo3DData.fMinimumY -= rfPadding * fDifference;
  CoordinateInfo3DData.fMaximumY += rfPadding * fDifference;
  fDifference = CoordinateInfo3DData.fMaximumZ - CoordinateInfo3DData.fMinimumZ;
  CoordinateInfo3DData.fMinimumZ -= rfPadding * fDifference;
  CoordinateInfo3DData.fMaximumZ += rfPadding * fDifference;
  CellInfo3DData.fWidthX = (CoordinateInfo3DData.fMaximumX - CoordinateInfo3DData.fMinimumX) / static_cast<float>(DivisionInfo3DData.iCountX);
  CellInfo3DData.fWidthY = (CoordinateInfo3DData.fMaximumY - CoordinateInfo3DData.fMinimumY) / static_cast<float>(DivisionInfo3DData.iCountY);
  CellInfo3DData.fWidthZ = (CoordinateInfo3DData.fMaximumZ - CoordinateInfo3DData.fMinimumZ) / static_cast<float>(DivisionInfo3DData.iCountZ);
}

void PiecewisePolynomial3D::CalculateCoefficient(const float* const* const* const& rpfSolution)
{
  float fSubCoefficient[DEGREE + 1][DEGREE + 1] = {
    {0.5F, 0.5F, 0.0F},
    {-1.0F, 1.0F, 0.0F},
    {0.5F, -1.0F, 0.5F}
  };

  for (int i = 0; i < DivisionInfo3DData.iCountX; i++)
    for (int j = 0; j < DivisionInfo3DData.iCountY; j++)
      for (int k = 0; k < DivisionInfo3DData.iCountZ; k++)
        for (int l = 0; l <= DEGREE; l++)
          for (int m = 0; m <= DEGREE; m++)
            for (int n = 0; n <= DEGREE; n++)
              for (int o = 0; o <= DEGREE; o++)
                for (int p = 0; p <= DEGREE; p++)
                  for (int q = 0; q <= DEGREE; q++)
                    pPolynomial3DMatrix[i][j][k][l][m][n] += rpfSolution[i + o][j + p][k + q] * fSubCoefficient[l][o] * fSubCoefficient[m][p] * fSubCoefficient[n][q];
}

void PiecewisePolynomial3D::Clear(void)
{
  if (pPolynomial3DMatrix == nullptr && pfCell == nullptr)
    return ;

  for (int i = 0; i < DivisionInfo3DData.iCountX; i++) {
    for (int j = 0; j < DivisionInfo3DData.iCountY; j++) {
      delete[] pPolynomial3DMatrix[i][j];
      delete[] pfCell[i][j];
    }

    delete[] pPolynomial3DMatrix[i];
    delete[] pfCell[i];
  }

  delete[] pPolynomial3DMatrix;
  delete[] pfCell;
}