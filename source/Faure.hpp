#ifndef FAURE_HPP_
#define FAURE_HPP_

#include <string>
#include <vector>

#include "IRandomGenerator.hpp"

class Faure final: public IRandomGenerator
{
public:
    Faure(unsigned dimension);
    ~Faure() = default;

    void save() const override;
    void load(std::string const& filename) override;
    void setOutputFile(std::string const& filename) override;
    std::string GetConfiguration() const override;

    void save(std::ostream& stream) const;
    void load(std::istream& stream);

    double Get() override;
    void Skip() override;
private:
    constexpr static std::uint64_t MAX_DIMENSION{ 17 };
    constexpr static std::uint64_t PRIMES_NUMBER{ 7 };
    constexpr static std::uint64_t PRIME_TABLE[PRIMES_NUMBER] = {
         2, 3, 5, 7, 11, 13, 17};

    std::uint64_t selectPrimeNumber() const;
    void initCTable();
    void nextPoint();

    inline std::uint64_t& c(std::uint64_t i, std::uint64_t j)
    {
        return c_[(logqPoint_ + 1) * i + j];
    }

    std::uint32_t dimension_;
    std::uint32_t currentDimension_;
    std::uint64_t pointId_{1};
    std::uint64_t q_{2};
    std::uint64_t qs_{2};
    std::uint64_t logqPoint_{0};
    double x_[MAX_DIMENSION]{};
    std::vector<std::uint64_t> c_;
    std::vector<std::uint64_t> tmp_;
    std::string outputFile_;
};

#endif // FAURE_HPP_
