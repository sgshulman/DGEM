#ifndef SAFIER_WIND_HPP_
#define SAFIER_WIND_HPP_

#include "IMatter.hpp"
#include "Predefines.hpp"

struct SafierWindModel;

class SafierWind final: public IMatter
{
public:
    SafierWind(char model, double mOut, double mStar, double h0, double rMin, double rMax, IDiskHumpCPtr hump);

    ~SafierWind() override = default;

    double density(Vector3d const& position) const override;

private:
    SafierWindModel const * const model_;
    double rho0_;
    double h0_;
    IDiskHumpCPtr hump_;
};

#endif
