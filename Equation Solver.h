#ifndef EQUATION_SOLVER
#define EQUATION_SOLVER

#include <vector>
#include <cstddef>

void LUDecomposition(const float* const* const& rpfMatrix, float** const& rpfLowerTriangularMatrix, float** const& rpfUpperTriangularMatrix, const std::size_t& rszSize);

void LUPDecomposition(const float* const* const& rpfMatrix, float** const& rpfLowerTriangularMatrix, float** const& rpfUpperTriangularMatrix, unsigned* const& rpuPivot, const std::size_t& rszSize);

void ForwardSubstitution(const float* const* const& rpfMatrix, std::vector<float>& rvecSolution, const std::size_t& rszSize);

void BackSubstitution(const float* const* const& rpfMatrix, std::vector<float>& rvecSolution, const std::size_t& rszSize);

void SolveEquation(const float* const* const& rpfMatrix, std::vector<float>& rvecSolution, const std::size_t& rszSize);

void SolveEquation(const float* const* const& rpfMatrix, unsigned* const& rpuPivot, std::vector<float>& rvecSolution, const std::size_t& rszSize);

#endif