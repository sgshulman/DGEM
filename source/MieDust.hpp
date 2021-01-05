#ifndef MIE_DUST_HPP_
#define MIE_DUST_HPP_

#include "IDust.hpp"
#include <string>
#include <vector>

// The dust with Mie theory scatterings.
// P1, P2, P3, P4, I, I/I0, and Pol are tabulated for different scattering angles
class MieDust: public IDust
{
public:
    MieDust(double albedo, std::string const& tableFile);

    void scatteringMatrixElements(
        double &p1,
        double &p2,
        double &p3,
        double &p4,
        double cosTheta) const override;

    double fraction(double cosTheta) const override;
    double cosRandomTheta(double v) const override;

    double albedo() const override
    {	return albedo_;	}

private:
    struct Data
    {
        double cosTheta;
        double p1;
        double p2;
        double p3;
        double p4;

        bool operator<(double cosineTheta) const
        {
            return cosTheta < cosineTheta;
        }
    };

    double const albedo_;
    std::vector<Data> table_;
    std::vector<double> accumulated_;
};

#endif