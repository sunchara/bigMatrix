#pragma once
#include <iostream>
#include <cstring>
#include <sstream>
#include <fstream>
#include <future>
#include <vector>
#include <optional>
#include <iomanip>
#include <algorithm>

template<size_t BLOCK_SIZE = 512>
struct SimpleBlockMatrix
{
    std::vector<double> elements = std::vector<double>(BLOCK_SIZE * BLOCK_SIZE, 0.);

    double& at(size_t x, size_t y)
    {
        return elements[x * BLOCK_SIZE + y];
    }
    double at(size_t x, size_t y) const
    {
        return elements[x * BLOCK_SIZE + y];
    }
    void multiply(double alpha)
    {
        for (double& x : elements)
            x *= alpha;
    }
    struct SimpleBlockMatrix<BLOCK_SIZE> operator+(const struct SimpleBlockMatrix<BLOCK_SIZE>& rhs) const
    {
        struct SimpleBlockMatrix<BLOCK_SIZE> m;

        for (size_t ind = 0; ind < BLOCK_SIZE * BLOCK_SIZE; ++ind)
            m.elements[ind] = elements[ind] + rhs.elements[ind];

        return m;
    }

    struct SimpleBlockMatrix<BLOCK_SIZE> operator-(const struct SimpleBlockMatrix<BLOCK_SIZE>& rhs) const
    {

        struct SimpleBlockMatrix<BLOCK_SIZE> m;

        for (size_t ind = 0; ind < BLOCK_SIZE * BLOCK_SIZE; ++ind)
            m.elements[ind] = elements[ind] - rhs.elements[ind];

        return m;
    }

    struct SimpleBlockMatrix<BLOCK_SIZE> operator*(const struct SimpleBlockMatrix<BLOCK_SIZE>& rhs)const
    {
        struct SimpleBlockMatrix<BLOCK_SIZE> dest;
        for (size_t i = 0; i < BLOCK_SIZE; ++i)
            for (size_t j = 0; j < BLOCK_SIZE; ++j)
            {
                double sum = 0, c = 0, y, t;
                for (size_t k = 0; k < BLOCK_SIZE; ++k)
                {
                    y = at(i, k) * rhs.at(k, j) - c;
                    t = sum + y;
                    sum = t;
                }
                dest.at(i, j) = sum;
            }

        return dest;
    }

    bool operator==(const struct SimpleBlockMatrix<BLOCK_SIZE>& rhs) const
    {
        return elements == rhs.elements;
    }

    double infNorm() const
    {
        return std::abs (*std::max_element(
            elements.begin(), elements.end(), [](double rhs, double lhs) { return std::abs(rhs) < std::abs(lhs); }));
    }

    bool isZero()
    {
        return std::all_of(elements.begin(), elements.end(),
            [](double x) { return std::abs(x) < std::numeric_limits<double>::epsilon; });
    }
};

template<size_t BLOCK_SIZE = 512>
struct BigSquareMatrix
{
    std::vector <std::optional<struct SimpleBlockMatrix<BLOCK_SIZE>>> elements;

    size_t dim = 0;
    size_t max_x = 0;
    size_t max_y = 0;
    BigSquareMatrix() = default;
    BigSquareMatrix(size_t n, size_t m, std::optional<size_t> dimOpt = std::nullopt) : max_x(n), max_y(m)
    {
        if (std::min(n, m) == 0) return;

        dim = 1;
        while (std::max(max_x, max_y) > dim * BLOCK_SIZE) dim *= 2;

        dim = std::max(dim, dimOpt.value_or(0));
        elements = std::vector <std::optional<struct SimpleBlockMatrix<BLOCK_SIZE>>>(dim * dim, std::nullopt);
    }

    double& at(size_t x, size_t y)
    {
        if (x >= max_x || y >= max_y)
            throw std::runtime_error("Out of range");

        size_t block_x = x / BLOCK_SIZE;
        size_t block_y = y / BLOCK_SIZE;

        if (!blockAt(block_x, block_y).has_value())
        {
            blockAt(block_x, block_y) = SimpleBlockMatrix<BLOCK_SIZE>();
        }

        return blockAt(block_x, block_y)->at(x % BLOCK_SIZE, y % BLOCK_SIZE);
    }

    double at(size_t x, size_t y) const
    {
        if (x >= max_x || y >= max_y)
            throw std::runtime_error("Out of range");

        size_t block_x = x / BLOCK_SIZE;
        size_t block_y = y / BLOCK_SIZE;

        if (!elements[block_x * dim + block_y].has_value())
        {
            return 0;
        }

        return elements[block_x * dim + block_y]->at(x % BLOCK_SIZE, y % BLOCK_SIZE);
    }
    const std::optional<struct SimpleBlockMatrix<BLOCK_SIZE>>& blockAt(size_t x, size_t y) const
    {
        return elements[x * dim + y];
    }
    std::optional<struct SimpleBlockMatrix<BLOCK_SIZE>>& blockAt(size_t x, size_t y)
    {
        return elements[x * dim + y];
    }

    struct BigSquareMatrix<BLOCK_SIZE> operator+(struct BigSquareMatrix<BLOCK_SIZE>& rhs) const
    {
        if (max_x != rhs.max_x || max_y != rhs.max_y || dim != rhs.dim)
            throw std::runtime_error("Bad size");

        struct BigSquareMatrix<BLOCK_SIZE> m(max_x, max_y);

        for (size_t ind = 0; ind < dim * dim; ++ind)
        {
            if (!elements[ind].has_value())
            {
                m.elements[ind] = rhs.elements[ind];
                continue;
            }
            if (!rhs.elements[ind].has_value())
            {
                m.elements[ind] = elements[ind];
                continue;
            }
            m.elements[ind] = elements[ind].value() + rhs.elements[ind].value();
        }

        return m;
    }


    void resize(size_t newDim)
    {
        if (dim >= newDim) return;

        auto newElements = std::vector <std::optional<struct SimpleBlockMatrix<BLOCK_SIZE>>>(newDim * newDim, std::nullopt);

        for (size_t x = 0; x < dim; ++x)
            for (size_t y = 0; y < dim; ++y)
                newElements[x * newDim + y] = std::move(elements[x * dim + y]);

        elements = newElements;

        dim = newDim;
    }

    struct BigSquareMatrix<BLOCK_SIZE> operator*(struct BigSquareMatrix<BLOCK_SIZE> rhs) const
    {
        if (max_y != rhs.max_x)
            throw std::runtime_error("Bad size");

        struct BigSquareMatrix<BLOCK_SIZE> m(max_x, rhs.max_y, std::max(dim, rhs.dim));
        if (dim == 1)
        {
            if (!blockAt(0, 0).has_value() || !rhs.blockAt(0, 0).has_value())
                return m;

            m.elements[0] = blockAt(0, 0).value() * rhs.blockAt(0, 0).value();
            return m;
        }

        size_t MAX_X = max_x, MAX_Y = rhs.max_y;

        struct BigSquareMatrix<BLOCK_SIZE> lb[2][2], rb[2][2], lhs = *this;
        lhs.resize(std::max(lhs.dim, rhs.dim));
        rhs.resize(std::max(lhs.dim, rhs.dim));

        lhs.toMinors(lb);
        rhs.toMinors(rb);
        size_t minorDim;


        static std::atomic_int countAsyncLaunches = 0;

        auto getPolicy = [&](bool& f)
            {
                if (countAsyncLaunches < 64)
                {
                    f = true;
                    ++countAsyncLaunches;
                    return std::launch::async;
                }
                f = false;
                return std::launch::deferred;
            };

        struct BigSquareMatrix<BLOCK_SIZE> d, d_1, d_2, v1, v2, h1, h2;
        {
            bool isLaunchedAsync[7] = { false,false,false,false,false,false,false };
            std::future<void> futures[7] = {
                std::async(getPolicy(isLaunchedAsync[0]), [&, A = (lb[0][0] + lb[1][1]), B = (rb[0][0] + rb[1][1])] {
                    d = A * B; if (isLaunchedAsync[0]) --countAsyncLaunches;
                    }),
                std::async(getPolicy(isLaunchedAsync[1]), [&, A = (lb[0][1] - lb[1][1]), B = (rb[1][0] + rb[1][1])] {d_1 = A * B; if (isLaunchedAsync[1]) --countAsyncLaunches; }),
                std::async(getPolicy(isLaunchedAsync[2]), [&, A = (lb[1][0] - lb[0][0]), B = (rb[0][0] + rb[0][1])] {d_2 = A * B; if (isLaunchedAsync[2]) --countAsyncLaunches; }),

                std::async(getPolicy(isLaunchedAsync[3]), [&, B = rb[1][0] - rb[0][0]] {v1 = lb[1][1] * B; if (isLaunchedAsync[3]) --countAsyncLaunches; }),
                std::async(getPolicy(isLaunchedAsync[4]), [&, B = rb[0][1] - rb[1][1]] {v2 = lb[0][0] * B; if (isLaunchedAsync[4]) --countAsyncLaunches; }),

                std::async(getPolicy(isLaunchedAsync[5]), [&, A = lb[0][0] + lb[0][1]] {h1 = A * rb[1][1]; if (isLaunchedAsync[5]) --countAsyncLaunches; }),
                std::async(getPolicy(isLaunchedAsync[6]), [&, A = lb[1][0] + lb[1][1]] {h2 = A * rb[0][0]; if (isLaunchedAsync[6]) --countAsyncLaunches; })
            };
            try
            {
                for (size_t i = 0; i < 7; ++i)
                    if (!isLaunchedAsync[i])
                        futures[i].get();

                for (size_t i = 0; i < 7; ++i)
                    if (isLaunchedAsync[i])
                        futures[i].get();

            }
            catch (std::runtime_error& what)
            {
                std::cerr << what.what() << std::endl;
                throw what;
            }
        }

        struct BigSquareMatrix<BLOCK_SIZE> mb[2][2] = {
            {d + d_1 + v1 - h1, v2 + h1},
            {v1 + h2, d + d_2 + v2 - h2}
        };

        m.fromMinors(std::move(mb));
        m.max_x = MAX_X;
        m.max_y = MAX_Y;
        return m;
    }
    struct BigSquareMatrix<BLOCK_SIZE> operator-(const struct BigSquareMatrix<BLOCK_SIZE>& rhs) const
    {
        if (max_x != rhs.max_x || max_y != rhs.max_y || dim != rhs.dim)
            throw std::runtime_error("Bad size");

        struct BigSquareMatrix<BLOCK_SIZE> m(max_x, max_y, dim);

        for (size_t ind = 0; ind < dim * dim; ++ind)
        {
            if (!rhs.elements[ind].has_value())
            {
                m.elements[ind] = elements[ind];
                continue;
            }

            if (!elements[ind].has_value())
            {
                m.elements[ind] = rhs.elements[ind];
                if (m.elements[ind].has_value())
                    m.elements[ind]->multiply(-1);

                continue;
            }
            m.elements[ind] = elements[ind].value() - rhs.elements[ind].value();
        }

        return m;
    }

    void toMinors(struct BigSquareMatrix<BLOCK_SIZE> minors[2][2])
    {
        minors[0][0] = minors[0][1] = minors[1][0] = minors[1][1] = BigSquareMatrix<BLOCK_SIZE>(dim * BLOCK_SIZE / 2, dim * BLOCK_SIZE / 2, dim / 2);

        size_t minorDim = dim / 2;

        for (size_t i = 0; i < minorDim * 2; ++i)
            for (size_t j = 0; j < minorDim * 2; ++j)
            {
                minors[i / minorDim][j / minorDim].blockAt(i % minorDim, j % minorDim) = std::move(blockAt(i, j));
                blockAt(i, j) = std::nullopt;
            }
    }

    void fromMinors(struct BigSquareMatrix<BLOCK_SIZE> minors[2][2])
    {
        size_t half_dim = minors[0][0].dim;
        dim = half_dim * 2;
        max_x = minors[0][0].max_x * 2;
        max_y = minors[0][0].max_y * 2;

        elements.resize(dim * dim);

        for (size_t i = 0; i < dim; ++i)
            for (size_t j = 0; j < dim; ++j)
            {
                blockAt(i, j) = std::move(minors[i / half_dim][j / half_dim].blockAt(i % half_dim, j % half_dim));
            }

    }
};

template<size_t BLOCK_SIZE = 512>
std::ostream& operator<<(std::ostream& out, const struct BigSquareMatrix<BLOCK_SIZE>& m)
{
    static std::mutex mut;
    std::unique_lock lk(mut);
    for (size_t x = 0; x < m.max_x; ++x)
    {
        for (size_t y = 0; y < m.max_y; ++y)
            out << m.at(x, y) << ' ';
        out << std::endl;
    }
    return out;
}

template<size_t BLOCK_SIZE = 512>
std::istream& operator>>(std::istream& in, struct BigSquareMatrix<BLOCK_SIZE>& m)
{
    for (size_t x = 0; x < m.max_x; ++x)
        for (size_t y = 0; y < m.max_y; ++y)
            in >> m.at(x, y);

    return in;
}