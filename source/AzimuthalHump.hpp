#ifndef AZIMUTHAL_HUMP_HPP_
#define AZIMUTHAL_HUMP_HPP_

#include "IDiskHump.hpp"

class AzimuthalHump final: public IDiskHump
{
public:
    AzimuthalHump(double h, double r, double sigma2, double sigma2azimuthalForward, double sigma2azimuthalBackward);

    ~AzimuthalHump() override = default;

    double hump(double value, Vector3d const& position) const override;

private:
    double const h_;
    double const r_;
    double const sigma2_;
    double const sigma2azForward_;
    double const sigma2azBackward_;
};

#endif
