#ifndef RANDOM_HPP_
#define RANDOM_HPP_

#include <cstdint>
#include <fstream>
#include <string>

// Random number generator
// L’Ecuyer, P. 1988, Communications of the ACM, vol. 31, pp. 742,774.
class Random
{
public:
    explicit Random(std::int32_t iseed)
        : iseed_(iseed)
    {};

    Random(Random const&) = delete;
    Random& operator=(Random const&) = delete;

    Random(Random&&) noexcept = default;
    Random& operator=(Random &&) noexcept = default;

    void save() const
    {
        if (!outputFile_.empty())
        {
            std::ofstream f(outputFile_);
            f << iseed_ << "\n";
            f << idum2_ << "\n";
            f << iy_ << "\n";
            for (int i=0; i!=32; ++i)
            {
                f << iv_[i] << "\t";
            }
        }
    }

    void load(std::string const& filename)
    {
        std::ifstream f(filename);
        f >> iseed_;
        f >> idum2_;
        f >> iy_;
        for (int i = 0; i != 32; ++i)
        {
            f >> iv_[i];
        }
    }

    void setOutputFile(std::string const& filename)
    {
        outputFile_ = filename;
    }

    double Get()
    {
        std::int32_t const IM1{ 2147483563 };
        std::int32_t const IM2{ 2147483399 };
        std::int32_t const IMM1{ IM1 - 1 };
        std::int32_t const IA1{ 40014 };
        std::int32_t const IA2{ 40692 };
        std::int32_t const IQ1{ 53668 };
        std::int32_t const IQ2{ 52774 };
        std::int32_t const IR1{ 12211 };
        std::int32_t const IR2{ 3791 };
        std::int32_t const NTAB{ 32 };
        std::int32_t const NDIV{ 1 + IMM1 / NTAB };
        double const EPS{ 1.2e-7 };
        double const AM{ 1. / IM1 };
        double const RNMX{ 1. - EPS};

        std::int32_t j, k;

        if (iseed_ <= 0)
        {
            iseed_ = (-iseed_ > 1) ? -iseed_ : 1;
            idum2_ = iseed_;
            for (j=NTAB+7; j>=0; --j)
            {
                k = iseed_ / IQ1;
                iseed_ = IA1 * (iseed_ - k * IQ1) - k * IR1;
                if (iseed_ < 0) iseed_ = iseed_ + IM1;
                if (j < NTAB) iv_[j] = iseed_;
            }
            iy_ = iv_[0];
        }
        k = iseed_ / IQ1;
        iseed_ = IA1 * (iseed_ - k * IQ1) - k * IR1;
        if (iseed_ < 0) iseed_ = iseed_ + IM1;
        k = idum2_ / IQ2;
        idum2_ = IA2 * (idum2_ - k * IQ2) - k * IR2;
        if (idum2_ < 0) idum2_ = idum2_ + IM2;
        j = iy_ / NDIV;
        iy_ = iv_[j] - idum2_;
        iv_[j] = iseed_;
        if(iy_ < 1) iy_ = iy_ + IMM1;

        return (AM * iy_ < RNMX) ? AM * iy_ : RNMX;
    }
private:
    std::int32_t iseed_;
    std::int32_t idum2_{ 123456789 };
    std::int32_t iy_{ 0 };
    std::int32_t iv_[32] = {0};
    std::string outputFile_;
};

#endif
