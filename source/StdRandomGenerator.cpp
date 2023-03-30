#include <fstream>

#include "DebugUtils.hpp"
#include "StdRandomGenerator.hpp"

IRandomGenerator* CreateStdRandomGenerator(RandomGeneratorType type, std::int32_t seed)
{
    switch (type)
    {
        case MINIMUM_STANDARD:
            return new StdRandomGenerator<std::minstd_rand>(seed);
        case MERSENNE_TWISTER:
            return new StdRandomGenerator<std::mt19937_64>(seed);
        case RANLUX:
            return new StdRandomGenerator<std::ranlux48>(seed);
        default:
            return nullptr;
    }
    return nullptr;
}


template<typename StdGenerator>
StdRandomGenerator<StdGenerator>::StdRandomGenerator(std::int32_t const seed)
    : generator_(static_cast<uint32_t>(seed))
    , distribution_(0.0, 1.0)
{}

template<typename StdGenerator>
void StdRandomGenerator<StdGenerator>::save() const
{
    if (!outputFile_.empty())
    {
        std::ofstream stream(outputFile_);
        save(stream);
    }
}

template<typename StdGenerator>
void StdRandomGenerator<StdGenerator>::save(std::ostream& stream) const
{
    stream << generator_;
}

template<typename StdGenerator>
void StdRandomGenerator<StdGenerator>::load(std::string const& filename)
{
    std::ifstream stream(filename);
    DATA_ASSERT(stream.is_open(), filename + " should exist.");
    load(stream);
}

template<typename StdGenerator>
void StdRandomGenerator<StdGenerator>::load(std::istream& stream)
{
    stream >> generator_;
}

template<typename StdGenerator>
void StdRandomGenerator<StdGenerator>::setOutputFile(std::string const& filename)
{
    outputFile_ = filename;
}

template<typename StdGenerator>
double StdRandomGenerator<StdGenerator>::Get()
{
    return distribution_(generator_);
}

template class StdRandomGenerator<std::minstd_rand>;
template class StdRandomGenerator<std::mt19937_64>;
template class StdRandomGenerator<std::ranlux48>;
