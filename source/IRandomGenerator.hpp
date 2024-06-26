#ifndef I_RANDOM_GENERATOR_HPP_
#define I_RANDOM_GENERATOR_HPP_

#include <string>
#include <cstdint>

enum class RandomGeneratorType : std::uint8_t
{
    NONE,
    MINIMUM_STANDARD,
    MERSENNE_TWISTER,
    RANLUX,
    LECUYER,
    HALTON,
    FAURE,
    SOBOL,
    NIEDERREITER,
    HAMMERSLEY
};

class IRandomGenerator
{
public:
    static IRandomGenerator* create(RandomGeneratorType type, std::int32_t seed, std::uint32_t dimension, std::uint64_t numberOfPoints);

    IRandomGenerator() = default;
    virtual ~IRandomGenerator() =default;

    IRandomGenerator(IRandomGenerator const&) = delete;
    IRandomGenerator& operator=(IRandomGenerator const&) = delete;

    IRandomGenerator(IRandomGenerator&&) noexcept = default;
    IRandomGenerator& operator=(IRandomGenerator &&) noexcept = default;

    virtual double Get() = 0;
    // Skip random number if we must get them by sequences of a certain length
    virtual void Skip() {}

    virtual void save() const = 0;
    virtual void load(std::string const& filename) = 0;
    virtual void setOutputFile(std::string const& filename) = 0;
    virtual std::string GetConfiguration() const = 0;
};

#endif