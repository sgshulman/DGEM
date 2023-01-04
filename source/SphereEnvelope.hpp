#ifndef SPHERE_ENVELOPE_HPP_
#define SPHERE_ENVELOPE_HPP_

#include "IMatter.hpp"
#include "Predefines.hpp"

class SphereEnvelope final: public IMatter
{
    public:
        SphereEnvelope(
                double rInner,
                double rOuter,
                double rho0,
                double r0,
                double alpha,
                MatterTranslationCPtr translation);

        ~SphereEnvelope() override = default;

        double density(Vector3d const& position) const override;

    private:
        double const rInner_;
        double const rOuter_;
        double const rho0_;
        double const r0_;
        double const alpha_;
        MatterTranslationCPtr translation_;
};

#endif
