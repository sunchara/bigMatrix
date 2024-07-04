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