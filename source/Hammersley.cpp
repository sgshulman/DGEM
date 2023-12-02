#include <cassert>
#include <fstream>
#include <iostream>

#include "DebugUtils.hpp"
#include "Hammersley.hpp"

constexpr uint64_t Hammersley::PRIME_NUMBERS[];

Hammersley::Hammersley(unsigned dimension, std::uint64_t numPoints)
    : curPoint_{ 0 }
    , numPoints_{ numPoints }
    , dimension_{ dimension }
    , currentDimension_{ dimension }
{
    DATA_ASSERT(dimension_ <= MAX_DIMENSION, "Hammersley Generator dimension is too high.");

    for (unsigned int i=0; i!=dimension_-1; ++i)
    {
        n_[i] = 0;
        d_[i] = 1;
    }
}

double Hammersley::Get()
{
    if (currentDimension_ == dimension_)
    {
        nextPoint();
        currentDimension_ = 0;
    }
    double const result = (currentDimension_ == 0)
        ? (static_cast<double>(curPoint_) - 0.5) / static_cast<double>(numPoints_)
        : static_cast<double>(n_[currentDimension_-1]) / static_cast<double>(d_[currentDimension_-1]);

    ++currentDimension_;
    return result;
}

void Hammersley::Skip()
{
    if (currentDimension_ != dimension_)
    {
        nextPoint();
        currentDimension_ = 0;
    }
}

void Hammersley::save() const
{
    if (!outputFile_.empty())
    {
        std::ofstream stream(outputFile_);
        save(stream);
    }
}

void Hammersley::load(std::string const& filename)
{
    std::ifstream stream(filename);
    DATA_ASSERT(stream.is_open(), filename + " should exist.");
    load(stream);
}

void Hammersley::setOutputFile(std::string const& filename)
{
    outputFile_ = filename;
}

std::string Hammersley::GetConfiguration() const
{
    return std::string("Hammersley. dimension = ") + std::to_string(dimension_) + ". Number of points = " + std::to_string(numPoints_);
}

void Hammersley::save(std::ostream& stream) const
{
    stream << curPoint_ << std::endl;
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

void Hammersley::load(std::istream& stream)
{
    stream >> curPoint_;
    for (unsigned int dim = 0; dim != dimension_; ++dim)
    {
        stream >> n_[dim];
    }

    for (unsigned int dim = 0; dim != dimension_; ++dim)
    {
        stream >> d_[dim];
    }
}

void Hammersley::nextPoint()
{
    ++curPoint_;

    for (unsigned int i=0; i != dimension_-1; ++i)
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
