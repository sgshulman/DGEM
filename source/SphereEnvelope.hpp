#ifndef SPHERE_ENVELOPE_HPP_
#define SPHERE_ENVELOPE_HPP_

#include "IMatter.hpp"

class SphereEnvelope : public IMatter
{
    public:
        SphereEnvelope(
                double rInner,
                double rOuter,
                double rho0,
                double r0,
                double alpha)
            : rInner_{ rInner }
            , rOuter_{ rOuter }
            , rho0_{ rho0 }
            , r0_{ r0 }
            , alpha_{ alpha }
        {}

        ~SphereEnvelope() override = default;

        double density(double x, double y, double z) const override;

    private:
        double const rInner_;
        double const rOuter_;
        double const rho0_;
        double const r0_;
        double const alpha_;
};

#endif
