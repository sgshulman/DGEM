#include "IRandomGenerator.hpp"
#include "Faure.hpp"
#include "Halton.hpp"
#include "LEcuyer.hpp"
#include "Niederreiter.hpp"
#include "Sobol.hpp"
#include "StdRandomGenerator.hpp"

IRandomGenerator* IRandomGenerator::create(RandomGeneratorType type, std::int32_t seed, std::uint32_t dimension)
{
    IRandomGenerator* rand{ nullptr };

    if (RandomGeneratorType::LECUYER == type)
    {
        rand = new LEcuyer(seed);
    }
    else if (RandomGeneratorType::HALTON == type)
    {
        rand = new Halton(dimension);
    }
    else if (RandomGeneratorType::FAURE == type)
    {
        rand = new Faure(dimension);
    }
    else if (RandomGeneratorType::SOBOL == type)
    {
        rand = new Sobol(dimension);
    }
    else if (RandomGeneratorType::NIEDERREITER == type)
    {
        rand = new Niederreiter(dimension);
    }
    else
    {
        rand = CreateStdRandomGenerator(type, seed);
    }

    return rand;
}

