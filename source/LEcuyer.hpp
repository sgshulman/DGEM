#ifndef LECUYER_HPP_
#define LECUYER_HPP_

#include "IRandomGenerator.hpp"
#include <cstdint>
#include <string>

// Random number generator
// L’Ecuyer, P. 1988, Communications of the ACM, vol. 31, pp. 742,774.
class LEcuyer final: public IRandomGenerator
{
public:
    explicit LEcuyer(std::int32_t iseed);

    void save() const override;
    void load(std::string const& filename) override;
    void setOutputFile(std::string const& filename) override;

    void save(std::ostream& stream) const;
    void load(std::istream& stream);

    double Get() override;

private:
    std::int32_t iseed_;
    std::int32_t idum2_{ 123456789 };
    std::int32_t iy_{ 0 };
    std::int32_t iv_[32] = {0};
    std::string outputFile_;
};

#endif
