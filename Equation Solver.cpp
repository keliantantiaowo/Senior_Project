#include <cmath>
#include "Display Result.h"
#include "Equation Solver.h"

using namespace std;

void LUDecomposition(const float* const* const& rpfMatrix, float** const& rpfLowerTriangularMatrix, float** const& rpfUpperTriangularMatrix, const size_t& rszSize)
{
  for (int i = 0; i < rszSize; i++)
    for (int j = 0; j < rszSize; j++)
      rpfUpperTriangularMatrix[i][j] = rpfMatrix[i][j];

  for (int i = 0; i < rszSize; i++) {
    rpfLowerTriangularMatrix[i][i] = 1.0F;

    for (int j = i + 1; j < rszSize; j++) {
      rpfLowerTriangularMatrix[i][j] = 0.0F;
      rpfLowerTriangularMatrix[j][i] = rpfUpperTriangularMatrix[j][i] / rpfUpperTriangularMatrix[i][i];
      rpfUpperTriangularMatrix[j][i] = 0.0F;

      for (int k = i + 1; k < rszSize; k++)
        rpfUpperTriangularMatrix[j][k] -= rpfLowerTriangularMatrix[j][i] * rpfUpperTriangularMatrix[i][k];
    }
  }
}

void LUPDecomposition(const float* const* const& rpfMatrix, float** const& rpfLowerTriangularMatrix, float** const& rpfUpperTriangularMatrix, unsigned* const& rpuPivot, const size_t& rszSize)
{
  for (int i = 0; i < rszSize; i++) {
    for (int j = 0; j < rszSize; j++)
      rpfUpperTriangularMatrix[i][j] = rpfMatrix[i][j];

    rpuPivot[i] = static_cast<unsigned>(i);
  }

  for (int i = 0; i < rszSize; i++) {
    float fMaximum = 0.0F;
    int iIndex = i;

    for (int j = i; j < rszSize; j++)
      if (fabs(rpfUpperTriangularMatrix[j][i]) > fMaximum) {
        fMaximum = fabs(rpfUpperTriangularMatrix[j][i]);
        iIndex = j;
      }

    float (*pfTempA), (*pfTempB);

    pfTempA = rpfLowerTriangularMatrix[i];
    pfTempB = rpfUpperTriangularMatrix[i];
    rpfLowerTriangularMatrix[i] = rpfLowerTriangularMatrix[iIndex];
    rpfLowerTriangularMatrix[iIndex] = pfTempA;
    rpfUpperTriangularMatrix[i] = rpfUpperTriangularMatrix[iIndex];
    rpfUpperTriangularMatrix[iIndex] = pfTempB;
    rpfLowerTriangularMatrix[i][i] = 1.0F;

    unsigned uTemp;
    
    uTemp = rpuPivot[i];
    rpuPivot[i] = rpuPivot[iIndex];
    rpuPivot[iIndex] = uTemp;

    for (int j = i + 1; j < rszSize; j++) {
      rpfLowerTriangularMatrix[i][j] = 0.0F;
      rpfLowerTriangularMatrix[j][i] = rpfUpperTriangularMatrix[j][i] / rpfUpperTriangularMatrix[i][i];
      rpfUpperTriangularMatrix[j][i] = 0.0F;

      for (int k = i + 1; k < rszSize; k++)
        rpfUpperTriangularMatrix[j][k] -= rpfLowerTriangularMatrix[j][i] * rpfUpperTriangularMatrix[i][k];
    }
  }
}

void ForwardSubstitution(const float* const* const& rpfMatrix, vector<float>& rvecSolution, const size_t& rszSize)
{
  rvecSolution[0] /= rpfMatrix[0][0];

  for (int i = 1; i < rszSize; i++) {
    for (int j = 0; j < i; j++)
      rvecSolution[i] -= rpfMatrix[i][j] * rvecSolution[j];

    rvecSolution[i] /= rpfMatrix[i][i];
  }
}

void BackSubstitution(const float* const* const& rpfMatrix, vector<float>& rvecSolution, const size_t& rszSize)
{
  rvecSolution[rszSize - 1] /= rpfMatrix[rszSize - 1][rszSize - 1];

  for (int i = 1; i < rszSize; i++) {
    for (int j = 0; j < i; j++)
      rvecSolution[rszSize - i - 1] -= rpfMatrix[rszSize - i - 1][rszSize - j - 1] * rvecSolution[rszSize - j - 1];

    rvecSolution[rszSize - i - 1] /= rpfMatrix[rszSize - i - 1][rszSize - i - 1];
  }
}

void SolveEquation(const float* const* const& rpfMatrix, vector<float>& rvecSolution, const size_t& rszSize)
{
  float** pfLowerTriangularMatrix = new float*[rszSize];
  float** pfUpperTriangularMatrix = new float*[rszSize];

  for (int i = 0; i < rszSize; i++) {
    pfLowerTriangularMatrix[i] = new float[rszSize];
    pfUpperTriangularMatrix[i] = new float[rszSize];
  }

  LUDecomposition(rpfMatrix, pfLowerTriangularMatrix, pfUpperTriangularMatrix, rszSize);

  ForwardSubstitution(pfLowerTriangularMatrix, rvecSolution, rszSize);
  BackSubstitution(pfUpperTriangularMatrix, rvecSolution, rszSize);

  for (int i = 0; i < rszSize; i++) {
    delete[] pfLowerTriangularMatrix[i];
    delete[] pfUpperTriangularMatrix[i];
  }

  delete[] pfLowerTriangularMatrix;
  delete[] pfUpperTriangularMatrix;
}

void SolveEquation(const float* const* const& rpfMatrix, unsigned* const& rpuPivot, vector<float>& rvecSolution, const size_t& rszSize)
{
  float** pfLowerTriangularMatrix = new float*[rszSize];
  float** pfUpperTriangularMatrix = new float*[rszSize];

  for (int i = 0; i < rszSize; i++) {
    pfLowerTriangularMatrix[i] = new float[rszSize];
    pfUpperTriangularMatrix[i] = new float[rszSize];
  }

  LUPDecomposition(rpfMatrix, pfLowerTriangularMatrix, pfUpperTriangularMatrix, rpuPivot, rszSize);

  float* pfOrder = new float[rszSize];

  for (int i = 0; i < rszSize; i++)
    pfOrder[i] = rvecSolution[rpuPivot[i]];

  for (int i = 0; i < rszSize; i++)
    rvecSolution[i] = pfOrder[i];

  delete[] pfOrder;

  ForwardSubstitution(pfLowerTriangularMatrix, rvecSolution, rszSize);
  BackSubstitution(pfUpperTriangularMatrix, rvecSolution, rszSize);

  for (int i = 0; i < rszSize; i++) {
    delete[] pfLowerTriangularMatrix[i];
    delete[] pfUpperTriangularMatrix[i];
  }

  delete[] pfLowerTriangularMatrix;
  delete[] pfUpperTriangularMatrix;
}