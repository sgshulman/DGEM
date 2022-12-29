#ifndef I_MATTER_HPP_
#define I_MATTER_HPP_

#include "Vector3d.hpp"

class IMatter
{
    public:
        IMatter() = default;
        virtual ~IMatter() = default;

        IMatter(IMatter const& other) = delete;
        IMatter(IMatter&& other) = delete;
        IMatter& operator=(IMatter const& other) = delete;
        IMatter& operator=(IMatter&& other) = delete;

        virtual double density(Vector3d const& position) const = 0;
};

#endif
