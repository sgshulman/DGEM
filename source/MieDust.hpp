#ifndef MIE_DUST_HPP_
#define MIE_DUST_HPP_

#include "IDust.hpp"
#include <string>

// The dust with Mie theory scatterings.
// P1, P2, P3, P4, I, I/I0, and Pol are tabulated for different scattering angles
class MieDust: public IDust
{
public:
    MieDust(std::string const& filename)
    {
        (void)filename;
    }

    void scatteringMatrixElements(
        double &p1,
        double &p2,
        double &p3,
        double &p4,
        double cosTheta) const override;

    double fraction(double cosTheta) const override;
    double cosRandomTheta(double v) const override;

    double albedo() const override
    {	return 0.0;	}

private:

};

#endif