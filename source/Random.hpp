#ifndef RANDOM_HPP_
#define RANDOM_HPP_

#include <cstdint>
#include <string>

// Random number generator
// L’Ecuyer, P. 1988, Communications of the ACM, vol. 31, pp. 742,774.
class Random
{
public:
    explicit Random(std::int32_t iseed);

    Random(Random const&) = delete;
    Random& operator=(Random const&) = delete;

    Random(Random&&) noexcept = default;
    Random& operator=(Random &&) noexcept = default;

    void save() const;
    void save(std::ostream& stream) const;
    void load(std::string const& filename);
    void load(std::istream& stream);
    void setOutputFile(std::string const& filename);

    double Get();

private:
    std::int32_t iseed_;
    std::int32_t idum2_{ 123456789 };
    std::int32_t iy_{ 0 };
    std::int32_t iv_[32] = {0};
    std::string outputFile_;
};

#endif
