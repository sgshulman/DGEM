#ifndef I_MATTER_HPP_
#define I_MATTER_HPP_

#include "Vector3d.hpp"

class IMatter
{
    public:
        virtual ~IMatter() = default;
        virtual double density(Vector3d const& position) const = 0;
};

#endif
