#include "BigMatrix.hpp"
#include "BigMatrixUtils.hpp"




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