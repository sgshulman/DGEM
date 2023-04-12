#ifndef STD_RANDOM_GENERATOR_HPP_
#define STD_RANDOM_GENERATOR_HPP_

#include "IRandomGenerator.hpp"
#include <cstdint>
#include <random>

IRandomGenerator* CreateStdRandomGenerator(RandomGeneratorType type, std::int32_t seed);

template<typename StdGenerator>
class StdRandomGenerator final: public IRandomGenerator
{
public:
    explicit StdRandomGenerator(std::int32_t seed);

    void save() const override;
    void load(std::string const& filename) override;
    void setOutputFile(std::string const& filename) override;
    std::string GetConfiguration() const override;

    void save(std::ostream& stream) const;
    void load(std::istream& stream);

    double Get() override;

private:
    StdGenerator generator_;
    std::uniform_real_distribution<double> distribution_;
    std::string outputFile_;
    std::int32_t seed_;
};

#endif // STD_RANDOM_GENERATOR_HPP_
