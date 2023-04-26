#include <cassert>
#include <fstream>

#include "DebugUtils.hpp"
#include "Faure.hpp"

namespace
{
    inline std::uint64_t log(std::uint64_t const base, std::uint64_t const x)
    {
        std::uint64_t res{ 0 };
        std::uint64_t power{ base };
        while (power <= x)
        {
            power *= base;
            ++res;
        }
        return res;
    }
} // namespace

constexpr std::uint64_t Faure::PRIME_TABLE[];

Faure::Faure(unsigned dimension)
    : dimension_{dimension > 0 ? dimension : 1}
    , currentDimension_{dimension_}
{
    DATA_ASSERT(dimension_ <= MAX_DIMENSION, "Faure Generator dimension is too high.");

    q_ = selectPrimeNumber();
    qs_ = q_;
    initCTable();
}

double Faure::Get()
{
    if (currentDimension_ == dimension_)
    {
        nextPoint();
        currentDimension_ = 0;
    }

    return x_[currentDimension_++];
}

void Faure::Skip()
{
    if (currentDimension_ != dimension_)
    {
        nextPoint();
        currentDimension_ = 0;
    }
}

void Faure::save() const
{
    if (!outputFile_.empty())
    {
        std::ofstream stream(outputFile_);
        save(stream);
    }
}

void Faure::load(std::string const& filename)
{
    std::ifstream stream(filename);
    DATA_ASSERT(stream.is_open(), filename + " should exist.");
    load(stream);
}

void Faure::setOutputFile(std::string const& filename)
{
    outputFile_ = filename;
}

std::string Faure::GetConfiguration() const
{
    return std::string("Faure. dimension = ") + std::to_string(dimension_);
}

void Faure::save(std::ostream& stream) const
{
    stream << pointId_ << "\t" << qs_/q_ << "\n";
    for (std::uint32_t dim = 0; dim != dimension_; ++dim)
    {
        stream << x_[dim] << "\t";
    }
}

void Faure::load(std::istream& stream)
{
    stream >> pointId_ >> qs_;
    for (std::uint32_t dim = 0; dim != dimension_; ++dim)
    {
        stream >> x_[dim];
    }

    logqPoint_ = log(q_, pointId_);
    initCTable();
}

std::uint64_t Faure::selectPrimeNumber() const
{
    for (std::size_t i=0; i!=PRIMES_NUMBER; ++i)
    {
        if (PRIME_TABLE[i] >= dimension_)
        {
            return PRIME_TABLE[i];
        }
    }
    DATA_ASSERT(false, std::to_string(PRIMES_NUMBER) + " prime numbers are not enough for dimension " + std::to_string(dimension_));
    return 0;
}

void Faure::initCTable()
{
    c_.resize((logqPoint_ + 1) * (logqPoint_ + 1));

    for (std::uint32_t n=0; n<=logqPoint_; ++n)
    {
        c(n, 0) = 1;
        c(n, n) = 1;
        for (std::uint32_t k=1; k<n; ++k)
        {
            c(n, k) = (c(n-1, k-1) + c(n-1, k)) % qs_;
        }
	}
    qs_ *= q_;
}

void Faure::nextPoint()
{
    std::uint64_t newLogqPoint = log(q_, pointId_);
    tmp_.resize(newLogqPoint + 1);

    if (newLogqPoint != logqPoint_)
    {
        logqPoint_ = newLogqPoint;
        initCTable();
    }

    std::uint64_t qi = qs_/q_;
    std::uint64_t idRem = pointId_;

    for (std::uint64_t i=0; i != newLogqPoint + 1; ++i)
    {
        qi /= q_;
        tmp_[newLogqPoint - i] = idRem / qi;
        idRem %= qi;
    }

    auto r = static_cast<double>(tmp_[newLogqPoint]);

    for (std::uint64_t i=0; i != newLogqPoint; ++i)
    {
        r = static_cast<double>(tmp_[newLogqPoint - i - 1]) + r / static_cast<double>(q_);
    }

    x_[0] = r / static_cast<double>(q_);

    for (std::uint32_t dim = 1; dim < dimension_; ++dim)
    {
        x_[dim] = 0.0;
        r = 1.0 / static_cast<double>(q_);

        for (std::uint64_t j = 0; j <= newLogqPoint; ++j)
        {
            // tmp[j] = sum_{i = j}^newLogqPoint (tmp(i) * c(i,j) ) mod q_.
            std::uint64_t sum = 0;
            for (std::uint64_t i = j; i <= newLogqPoint; ++i)
            {
                sum += tmp_[i] * c(i,j);
            }
            tmp_[j] = sum % q_;
            x_[dim] += static_cast<double>(tmp_[j]) * r;
            r = r / static_cast<double>(q_);
        }
    }
    ++pointId_;
}
