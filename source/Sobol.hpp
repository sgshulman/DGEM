#ifndef SOBOL_HPP_
#define SOBOL_HPP_

#include <cmath>
#include <string>
#include <vector>

#include "IRandomGenerator.hpp"

class Sobol final: public IRandomGenerator
{
public:
    Sobol(unsigned dimension);
    ~Sobol() = default;

    void save() const override;
    void load(std::string const& filename) override;
    void setOutputFile(std::string const& filename) override;
    std::string GetConfiguration() const override;

    void save(std::ostream& stream) const;
    void load(std::istream& stream);

    double Get() override;
    void Skip() override;
private:
    constexpr static unsigned int MAX_DIMENSION = 15;
    constexpr static unsigned int POLINOMIAL_NUMBER = MAX_DIMENSION - 1;
    constexpr static unsigned int MAX_DEGREE = 6;
    constexpr static unsigned int BIT_COUNT = 64;
    static double const POW2;

    using table_t = std::uint8_t;
    using point_t = std::uint64_t;

    // Joe and Kuo, SIAM J. Sci. Comput. 30, 2635 (2008).

    // p(z) = a_0 + a_1 z + a_2 z^2 + ... a_BIT_COUNT z^BIT_COUNT,
    // here a_i is the i-th bit of P_TABLE[j], j is the polynomial index.
    // our p_table element in binary form is 1(a_j from Joe and Kuo)1
    constexpr static table_t P_TABLE[POLINOMIAL_NUMBER] = {
         3,  7, 11, 13, 19,
        25, 37, 41, 47, 55,
        59, 61, 67, 91};

    // initial m_i values for i=0..degree of the POLINOMIAL_NUMBER primitive polynomials
    constexpr static table_t MI_TABLE[MAX_DEGREE * POLINOMIAL_NUMBER] = {
        1,  0,  0,  0,  0,  0,
        1,  3,  0,  0,  0,  0,

        1,  3,  1,  0,  0,  0,
        1,  1,  1,  0,  0,  0,

        1,  1,  3,  3,  0,  0,
        1,  3,  5, 13,  0,  0,

        1,  1,  5,  5, 17,  0,
        1,  1,  5,  5,  5,  0,
        1,  1,  7, 11, 19,  0,
        1,  1,  5,  1,  1,  0,
        1,  1,  1,  3, 11,  0,
        1,  3,  5,  5, 31,  0,

        1,  3,  3,  9,  7, 49,
        1,  1,  1, 15, 21, 21
    };

    static inline table_t m_initial(unsigned int dimension, unsigned int degree)
    {
      return MI_TABLE[MAX_DEGREE * dimension + degree];
    }

    inline point_t& m(unsigned int i, unsigned int dimension)
    {
        return m_[dimension_*i + dimension];
    }

    void initMiTable();
    void nextPoint();

    unsigned int dimension_;
    unsigned int currentDimension_;
    std::uint64_t pointId_{0};
    point_t x_[MAX_DIMENSION]{};
    std::vector<point_t> m_;
    std::string outputFile_;
};

#endif // SOBOL_HPP_
