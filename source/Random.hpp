#ifndef RANDOM_HPP_
#define RANDOM_HPP_

#include <cstdint>

// Random number generator
// L’Ecuyer, P. 1988, Communications of the ACM, vol. 31, pp. 742,774.
class Random
{
public:
    explicit Random (std::int32_t iseed)
        : iseed_(iseed)
    {};

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
            idum2 = iseed_;
            for (j=NTAB+7; j>=0;--j)
            {
                k = iseed_/IQ1;
                iseed_ = IA1*(iseed_-k*IQ1)-k*IR1;
                if (iseed_ < 0) iseed_ = iseed_+IM1;
                if (j < NTAB) iv[j] = iseed_;
            }
            iy=iv[0];
        }
        k = iseed_/IQ1;
        iseed_ = IA1*(iseed_-k*IQ1)-k*IR1;
        if (iseed_ < 0) iseed_=iseed_+IM1;
        k = idum2/IQ2;
        idum2 = IA2*(idum2-k*IQ2)-k*IR2;
        if (idum2 < 0) idum2 = idum2+IM2;
        j = iy/NDIV;
        iy = iv[j]-idum2;
        iv[j] = iseed_;
        if(iy < 1) iy = iy+IMM1;

        return (AM*iy < RNMX) ? AM*iy : RNMX;
    }
private:
    std::int32_t iseed_;
    std::int32_t idum2{ 123456789 };
    std::int32_t iy{ 0 };
    std::int32_t iv[32] = {0};
};

#endif
