#pragma once
#include "BigMatrix.hpp"

template<size_t BLOCK_SIZE = 512>
bool equalsByEps(const struct SimpleBlockMatrix<BLOCK_SIZE>& lhs, const struct SimpleBlockMatrix<BLOCK_SIZE>& rhs)
{
    for (size_t ind = 0; ind < BLOCK_SIZE * BLOCK_SIZE; ++ind)
    {
        double differenceAbs = std::abs(
            lhs.elements[ind] - rhs.elements[ind]);

        double minAbs = std::min(std::abs(lhs.elements[ind], std::abs(lhs.elements[ind])));

        if (differenceAbs > std::numeric_limits<double>::epsilon * minAbs)
            return false;
    }
    return true;
}

template<size_t BLOCK_SIZE = 512>
struct BigSquareMatrix<BLOCK_SIZE> getToeplizMatrix(size_t n, size_t m)
{
    struct BigSquareMatrix<BLOCK_SIZE> mm(n, m);

    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < m; ++j)
            mm.at(i, j) = i + j;

    return mm;
}

template<size_t BLOCK_SIZE = 512>
struct BigSquareMatrix<BLOCK_SIZE> getIdentityMatrix(size_t n, size_t m)
{
    struct BigSquareMatrix<BLOCK_SIZE> mm(n, m);

    for (size_t i = 0; i < std::min(n, m); ++i)
        mm.at(i, i) = 1;

    return mm;
}