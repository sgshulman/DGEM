#include <cstdint>
#include <fstream>
#include <string>

#include "DebugUtils.hpp"
#include "Random.hpp"

Random::Random(std::int32_t const iseed)
    : iseed_(iseed)
{}


void Random::save() const
{
    if (!outputFile_.empty())
    {
        std::ofstream stream(outputFile_);
        save(stream);
    }
}

void Random::save(std::ostream& stream) const
{
    stream << iseed_ << "\n";
    stream << idum2_ << "\n";
    stream << iy_ << "\n";
    for (int i=0; i!=32; ++i)
    {
        stream << iv_[i] << "\t";
    }
}

void Random::load(std::string const& filename)
{
    std::ifstream stream(filename);
    DATA_ASSERT(stream.is_open(), filename + " should exist.");
    load(stream);
}

void Random::load(std::istream& stream)
{
    stream >> iseed_;
    stream >> idum2_;
    stream >> iy_;
    for (int i = 0; i != 32; ++i)
    {
        stream >> iv_[i];
    }
}

void Random::setOutputFile(std::string const& filename)
{
    outputFile_ = filename;
}

double Random::Get()
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
