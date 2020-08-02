#ifndef I_DISK_HUMP_HPP_
#define I_DISK_HUMP_HPP_

#include "Vector3d.hpp"

class IDiskHump
{
public:
    virtual ~IDiskHump() = default;
    virtual double hump(double value, Vector3d const& position) const = 0;
};

#endif
