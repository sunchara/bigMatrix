#include "BigMatrix.hpp"

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


int main()
{
    try
    {
        std::cout << getToeplizMatrix<512>(10, 9) * getToeplizMatrix<512>(9, 10);
    }
    catch (...)
    {
        return -1;
    }
    return 0;
}