#include <cassert>
#include <cmath>
#include <fstream>

#include "DebugUtils.hpp"
#include "Niederreiter.hpp"

namespace
{
    inline std::uint32_t most_significant_bit(std::uint64_t i)
    {
       assert(i > 0);

       std::uint32_t index = 0;
       while (i != 0U)
       {
          ++index;
          i >>= 1U;
       }
       return --index;
    }

    inline std::uint32_t rightmost_zero_bit(std::uint64_t i)
    {
        std::uint32_t index = 0;

        if (i > 0)
        {
            std::uint64_t value = i-1;
            while ((value & 1U) != 0U)
            {
                value >>= 1U;
                ++index;
            }
        }
        return index;
    }

    inline std::bitset<128> multiply(std::uint32_t p, std::bitset<128> bOld)
    {
        std::bitset<128> result;

        for (; p != 0U; p >>= 1U)
        {
            if ((p & 1U) != 0U)
            {
                result ^= bOld;
            }
            bOld <<= 1U;
        }

        return result;
    }
} // namespace

constexpr std::uint8_t Niederreiter::P_TABLE[];
double const Niederreiter::POW2 = std::pow(2.0, BIT_COUNT);

Niederreiter::Niederreiter(unsigned dimension)
    : dimension_{dimension > 0 ? dimension : 1}
    , currentDimension_{dimension_}
{
    DATA_ASSERT(dimension_ <= MAX_DIMENSION, "Niederreiter Generator dimension is too high.");
    initCTable();
}

double Niederreiter::Get()
{
    if (currentDimension_ == dimension_)
    {
        nextPoint();
        currentDimension_ = 0;
    }

    return static_cast<double>(x_[currentDimension_++]) / POW2;
}

void Niederreiter::Skip()
{
    if (currentDimension_ != dimension_)
    {
        nextPoint();
        currentDimension_ = 0;
    }
}

void Niederreiter::save() const
{
    if (!outputFile_.empty())
    {
        std::ofstream stream(outputFile_);
        save(stream);
    }
}

void Niederreiter::load(std::string const& filename)
{
    std::ifstream stream(filename);
    DATA_ASSERT(stream.is_open(), filename + " should exist.");
    load(stream);
}

void Niederreiter::setOutputFile(std::string const& filename)
{
    outputFile_ = filename;
}

std::string Niederreiter::GetConfiguration() const
{
    return std::string("Niederreiter. dimension = ") + std::to_string(dimension_);
}

void Niederreiter::save(std::ostream& stream) const
{
    stream << pointId_ << "\n";
    for (unsigned int dim = 0; dim != dimension_; ++dim)
    {
        stream << x_[dim] << "\t";
    }
}

void Niederreiter::load(std::istream& stream)
{
    stream >> pointId_;
    for (unsigned int dim = 0; dim != dimension_; ++dim)
    {
        stream >> x_[dim];
    }
}


// Compute V(q,r) values based on BFN1988 sections 2.3 and 3.3.
inline std::bitset<128> Niederreiter::compute_v(std::bitset<128> const& b, std::uint32_t e, std::uint32_t m, std::uint32_t mOld)
{
    std::bitset<128> v;

    // BFN1988 Section 2.3
    // (3) ... using v_i = 0 for 0 <= i <= m - 2
    // v_{m-1} = 1,
    // BFN1988 Section 3.3
    // Our program currently sets each K_q equal to eq.
    // This has the effect of setting all unrestricted values of v to 1."
    for (std::size_t i = mOld; i < m; ++i)
    {
        v.set(i);
    }

    // BFN1988 Section 2.3
    // (3) ... and v_i = ADD_{k=1}^m MUL(b_{m-k}, v_{i-k}) for m <= i <= R+e-2
    for (std::size_t i = m; i < e + BIT_COUNT - 1; ++i)
    {
        bool v_i = false;
        for (std::uint32_t k = 0; k < m; ++k)
        {
            v_i = (v_i != (b.test(k) && v[i + k - m]));
        }
        v[i] = v_i;
    }

    return v;
}

// We use algorithm described in Bratley, Fox, Niederreiter, ACM Trans. Model. Comp. Sim. 2, 195 (1992) - BFN1988
// We use Niederreiter base 2 sequence, hence we need only 1 bit to store polynomial coefficients
// Moreover, ADD and MUL operations for this coefficient are bitwise XOR and AND respectively
void Niederreiter::initCTable()
{
    c_.resize(BIT_COUNT * dimension_);

    for (std::uint32_t dim = 0; dim < dimension_; ++dim)
    {
        // BFN1988 Section 2.3
        // (1) Choose a suitable monic polynomial p(x) with coefficients in F_b
        // Let the degree of p(x) be e >= 1
        std::uint32_t const p = P_TABLE[dim];
        std::uint32_t const e = most_significant_bit(p);

        // bitset to hold b(x) = p(x)^{q+1}
        std::bitset<128> b;
        std::uint32_t bDegree = 0;
        b.set(bDegree);

        for (std::uint32_t j=1; j!=BIT_COUNT;)
        {
            // BFN1988 Section 2.3
            // (3) ... Calculate first
            // b(x) = p(x)^{q+1} = x^m - b_{m-1}x^{m-1} -  ... - b_0
            // NB: here b_i are negative but in pb they are positive.
            std::uint32_t bDegreeOld = bDegree;
            b = multiply(p, b);
            bDegree += e;
            DATA_ASSERT(bDegree < e + BIT_COUNT - 1, "Niederreiter::initCTable b(x)=p(x)^{q+1} degree is to high.");

            // DFN1988 Section 2.3
            // (3) ... and then the elements v_i
            // We need BIT_COUNT + maximum degree of poly from P_TABLE bits. 128 bits are enough
            std::bitset<128> v = compute_v(b, e, bDegree, bDegreeOld);

            //  DFN1988 Section 2.3
            // (4) ... Increment u. If j < R = BIT_COUNT, go to step (2)
            // (2) Increment j. If u = e do step (3); otherwise, go to step (4)
            for (std::uint32_t u = 0; j!=BIT_COUNT && u != e; ++u, ++j)
            {
                // mask to change j-th bit of c(dim, r) which holds c(dim, j, r) values
                // In DFN1988 zero bit of the number is the "leftmost" one
                std::uint64_t j_bit = static_cast<std::uint64_t>(1U) << (BIT_COUNT - j);

                // DFN1988 Section 2.3
                // (4) For 0 <= r <= R-1 set c_{jr} = v_{r+u}
                for (std::uint32_t r = 0; r != BIT_COUNT; ++r)
                {
                    std::uint64_t const value = c(r, dim);
                    c(r, dim) = (value & ~j_bit) | (static_cast<std::uint64_t>(v[r+u]) * j_bit);
                }
            }
        }
    }
}


void Niederreiter::nextPoint()
{
    ++pointId_;
    std::uint32_t const rzb = rightmost_zero_bit(pointId_);

    for (std::uint32_t dim = 0; dim != dimension_; ++dim)
    {
        x_[dim] ^= c(rzb, dim);
    }
}
