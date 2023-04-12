#include <cassert>
#include <fstream>

#include "DebugUtils.hpp"
#include "Halton.hpp"

constexpr uint64_t Halton::PRIME_NUMBERS[];

Halton::Halton(unsigned dimension)
    : dimension_{ dimension }
    , currentDimension_{ dimension }
{
    DATA_ASSERT(dimension_ <= MAX_DIMENSION, "Halton Generator dimension is too high.");

    for (unsigned int i=0; i!=dimension_; ++i)
    {
        n_[i] = 0;
        d_[i] = 1;
    }
}

double Halton::Get()
{
    if (currentDimension_ == dimension_)
    {
        nextPoint();
        currentDimension_ = 0;
    }
    double const result = static_cast<double>(n_[currentDimension_]) / static_cast<double>(d_[currentDimension_]);
    ++currentDimension_;
    return result;
}

void Halton::Skip()
{
    if (currentDimension_ != dimension_)
    {
        nextPoint();
        currentDimension_ = 0;
    }
}

void Halton::save() const
{
    if (!outputFile_.empty())
    {
        std::ofstream stream(outputFile_);
        save(stream);
    }
}

void Halton::load(std::string const& filename)
{
    std::ifstream stream(filename);
    DATA_ASSERT(stream.is_open(), filename + " should exist.");
    load(stream);
}

void Halton::setOutputFile(std::string const& filename)
{
    outputFile_ = filename;
}

std::string Halton::GetConfiguration() const
{
    return std::string("Halton. dimension = ") + std::to_string(dimension_);
}

void Halton::save(std::ostream& stream) const
{
    for (unsigned int dim = 0; dim != dimension_; ++dim)
    {
        stream << n_[dim] << "\t";
    }
    stream << std::endl;
    for (unsigned int dim = 0; dim != dimension_; ++dim)
    {
        stream << d_[dim] << "\t";
    }
}

void Halton::load(std::istream& stream)
{
    for (unsigned int dim = 0; dim != dimension_; ++dim)
    {
        stream >> n_[dim];
    }

    for (unsigned int dim = 0; dim != dimension_; ++dim)
    {
        stream >> d_[dim];
    }
}

void Halton::nextPoint()
{
    for (unsigned int i=0; i!=dimension_; ++i)
    {
        uint64_t current_base = PRIME_NUMBERS[i];
        uint64_t x = d_[i] - n_[i];
        if (x == 1)
        {
            n_[i] = 1;
            d_[i] *= current_base;
        }
        else
        {
            uint64_t y = d_[i] / current_base;
            while (x <= y)
            {
                 y /= current_base;
            }
            n_[i] = (current_base + 1) * y - x;
        }
    }
}
