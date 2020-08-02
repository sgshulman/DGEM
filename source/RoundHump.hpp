#ifndef ROUND_HUMP_HPP_
#define ROUND_HUMP_HPP_

#include "IDiskHump.hpp"

class RoundHump : public IDiskHump
{
public:
    RoundHump(double h, double r, double sigma2);

    ~RoundHump() override = default;

    double hump(double value, Vector3d const& position) const override;

private:
    double const h_;
    double const r_;
    double const sigma2_;
};

#endif
