#ifndef MATTER_TRANSLATION_HPP_
#define MATTER_TRANSLATION_HPP_

#include "Vector3d.hpp"

// Translation to move IMatter elements from the coordinate system origin
// operator() converts position in global coordinate system to local position in IMatter system
// Thus operator() provides the reverse conversation
// The direct transformation order:
// 1. precession
// 2. nutation
// 3. intrinsicRotation
// 4. translation
class MatterTranslation final
{
    public:
        MatterTranslation(double precession, double nutation, double intrinsicRotation, Vector3d translation);

        Vector3d operator()(Vector3d const& position) const;

    private:
        double const cosPrecession_;
        double const sinPrecession_;
        double const cosNutation_;
        double const sinNutation_;
        double const cosIntrinsicRotation_;
        double const sinIntrinsicRotation_;
        Vector3d const translation_;
};

#endif
