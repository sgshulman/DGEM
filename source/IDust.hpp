#ifndef I_DUST_H
#define I_DUST_H

// Common interface for different dust properties
class IDust
{
    public:
        virtual ~IDust() = default;

        // TODO: use real matrix for future extensions and
        // refactor Photon::Stokes
        virtual void scatteringMatrixElements(
            double &p1,
            double &p2,
            double &p3,
            double &p4,
            double cosTheta) const = 0;

        virtual double fraction(double cosTheta) const = 0;
        virtual double cosRandomTheta(double v) const = 0;
        virtual double albedo() const = 0;
};

#endif //I_DUST_H
