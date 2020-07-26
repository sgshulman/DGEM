#ifndef SAFIER_WIND_HPP_
#define SAFIER_WIND_HPP_

#include "IMatter.hpp"

struct Model;

class SafierWind : public IMatter
{
public:
    SafierWind(char model, double mOut, double mStar, double h0, double rMin, double rMax);

    ~SafierWind() override = default;

    double density(double x, double y, double z) const override;

private:
    Model const * const model_;
    double rho0_;
    double h0_;
};

#endif
