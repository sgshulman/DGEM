#ifndef HALTON_HPP_
#define HALTON_HPP_

#include <string>

#include "IRandomGenerator.hpp"

class Halton final: public IRandomGenerator
{
public:
    Halton(unsigned dimension);
    ~Halton() = default;

    void save() const override;
    void load(std::string const& filename) override;
    void setOutputFile(std::string const& filename) override;
    std::string GetConfiguration() const override;

    void save(std::ostream& stream) const;
    void load(std::istream& stream);

    double Get() override;
    void Skip() override;

private:
    constexpr static std::uint32_t MAX_DIMENSION = 15;

    constexpr static std::uint64_t PRIME_NUMBERS[MAX_DIMENSION] = {
         2,     3,      5,      7,      11,
         13,   17,     19,     23,      29,
         31,   37,     41,     43,      47
    };

    void nextPoint();

    std::uint32_t dimension_;
    std::uint32_t currentDimension_;
    std::uint64_t n_[MAX_DIMENSION]{};
    std::uint64_t d_[MAX_DIMENSION]{};
    std::string outputFile_;
};

#endif // HALTON_HPP_
