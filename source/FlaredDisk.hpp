#ifndef FLARED_DISK_HPP_
#define FLARED_DISK_HPP_

#include "IMatter.hpp"
#include "Predefines.hpp"

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
                double beta,
                IMatterCPtr wind,
                MatterTranslationCPtr translation,
                IDiskHumpCPtr hump);

        ~FlaredDisk() override = default;

        double density(Vector3d const& position) const override;

    private:
        double const rInner_;
        double const rOuter_;
        double const rho0_;
        double const h0_;
        double const r0_;
        double const alpha_;
        double const beta_;
        IMatterCPtr wind_;
        MatterTranslationCPtr translation_;
        IDiskHumpCPtr hump_;
};

#endif
