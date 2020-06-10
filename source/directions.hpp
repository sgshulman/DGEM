#ifndef DIRECTIONS_HPP_
#define DIRECTIONS_HPP_

#include <cstdint>
#include "Vector3d.hpp"

class Directions {
    public:
        Directions(uint32_t NumOfDirectionsLevels);

        ~Directions()
        {
            delete[] dots_;
            delete[] w_;
        }

        Directions(Directions const &) = delete;
        Directions& operator=(Directions const &) = delete;

        inline uint64_t NumOfDirections() const
        {	return ( NumOfDirections_ );	}

        inline Vector3d direction(uint64_t index) const
        {
            return dots_[index];
        }

        inline double w(uint64_t index) const
        {	return w_[index];	}

    private:
        Vector3d *dots_{ nullptr };
        double   *w_{ nullptr };
        uint64_t NumOfDirections_{ 0 };
};

#endif
