#include <cassert>
#include <iostream>

#include "Sobol.hpp"

namespace
{
    inline unsigned int most_significant_bit(std::uint64_t i)
    {
       assert(i > 0);

       unsigned int index = 0;
       while (i)
       {
          ++index;
          i >>= 1;
       }
       return --index;
    }

    inline unsigned int rightmost_zero_bit(std::uint64_t i)
    {
        unsigned int index = 0;

        if (i > 0)
        {
            unsigned value = i-1;
            while (value & 1)
            {
                value >>= 1;
                ++index;
            }
        }
        return index;
    }
}

constexpr Sobol::table_t Sobol::P_TABLE[];
constexpr Sobol::table_t Sobol::MI_TABLE[];

Sobol::Sobol(unsigned dimension)
    : dimension_{dimension}
    , currentDimension_{dimension}
{
    initMiTable();
}

double Sobol::Get()
{
    if (currentDimension_ == dimension_)
    {
        nextPoint();
        currentDimension_ = 0;
    }

    return static_cast<double>(x_[currentDimension_++]) / POW2;
}

void Sobol::initMiTable()
{
    m_.resize(BIT_COUNT * dimension_);

    for (unsigned int i = 0; i != BIT_COUNT; ++i)
    {
        m(i, 0) = static_cast<point_t>(1);
    }

    for (unsigned int dim = 1; dim < dimension_; ++dim)
    {
        unsigned int const currentPoly = P_TABLE[dim-1];
        unsigned int const degree = most_significant_bit(currentPoly);

        // initial m_i from table
        for (unsigned int i = 0; i != degree; ++i)
        {
            m(i, dim) = m_initial(dim-1, i);
        }

        // other m_i. Bratley and Fox, TOMS 14, 88 (1988)
        for (unsigned int i = degree; i < BIT_COUNT; ++i)
        {
            unsigned int polyCopy = currentPoly;
            // m_i = 2*a_i*m_{i-1} ^ 2^2*a_2*m_{i-2} ^ ... ^ 2^degree*m_{i-degree} ^ m_{i-degree}
            m(i, dim) = m(i - degree, dim);
            for (unsigned int shift = degree; shift > 0; --shift)
            {
                unsigned int const a_i = polyCopy & 1;
                m(i, dim) ^= (a_i * m(i - shift, dim)) << shift;
                polyCopy >>= 1;
            }
        }
    }

    unsigned int shift = BIT_COUNT - 1;
    for (unsigned int i = 0; i != BIT_COUNT - 1; ++i)
    {
        for (unsigned int dim = 0; dim != dimension_; ++dim)
        {
            m(i, dim) <<= shift;
        }
        --shift;
    }
}

void Sobol::nextPoint()
{
    ++pointId_;
    unsigned int const c = rightmost_zero_bit(pointId_);

    for (unsigned int dim = 0; dim != dimension_; ++dim)
    {
        x_[dim] ^= m(c, dim);
    }
}
