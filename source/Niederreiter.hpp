#ifndef NIEDERREITER_HPP_
#define NIEDERREITER_HPP_

#include <bitset>
#include <string>
#include <vector>

#include "IRandomGenerator.hpp"

// Niederreiter base 2 sequence based on
// Bratley, Fox, Niederreiter, ACM Trans. Model. Comp. Sim. 2, 195 (1992)
class Niederreiter final: public IRandomGenerator
{
public:
    Niederreiter(unsigned dimension);
    ~Niederreiter() = default;

    void save() const override;
    void load(std::string const& filename) override;
    void setOutputFile(std::string const& filename) override;

    void save(std::ostream& stream) const;
    void load(std::istream& stream);

    double Get() override;
    void Skip() override;
private:
    constexpr static std::uint32_t MAX_DIMENSION = 15;
    constexpr static std::uint32_t BIT_COUNT = 64;
    static double const POW2;

    constexpr static std::uint8_t P_TABLE[MAX_DIMENSION] = {
         2,  3,  7, 11, 13,
        19, 25, 31, 37, 41,
        47, 55, 59, 61, 67};

    void initCTable();
    void nextPoint();
    // Compute V(q,r) values based on BFN1988 sections 2.3 and 3.3.
    static inline std::bitset<128> compute_v(std::bitset<128> const& b, std::uint32_t e, std::uint32_t m, std::uint32_t mOld);

    inline uint64_t& c(unsigned int bit, unsigned int dimension)
    {
        return c_[dimension_*bit + dimension];
    }

    std::uint32_t dimension_;
    std::uint32_t currentDimension_;
    std::uint64_t pointId_{0};
    std::uint64_t x_[MAX_DIMENSION]{};
    std::vector<std::uint64_t> c_;
    std::string outputFile_;
};

#endif // NIEDERREITER_HPP_
