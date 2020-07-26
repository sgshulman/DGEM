#ifndef FLARED_DISK_HPP_
#define FLARED_DISK_HPP_

#include "IMatter.hpp"

class FlaredDisk : public IMatter
{
    public:
        FlaredDisk(
                double rInner,
                double rOuter,
                double rho0,
                double h0,
                double r0,
                double alpha,
                double beta)
            : rInner_{ rInner }
            , rOuter_{ rOuter }
            , rho0_{ rho0 }
            , h0_{ h0 }
            , r0_{ r0 }
            , alpha_{ alpha }
            , beta_{ beta }
        {}

        ~FlaredDisk() override = default;

        double density(double x, double y, double z) const override;

    private:
        double const rInner_;
        double const rOuter_;
        double const rho0_;
        double const h0_;
        double const r0_;
        double const alpha_;
        double const beta_;
};

#endif
