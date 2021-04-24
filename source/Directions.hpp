#ifndef DIRECTIONS_HPP_
#define DIRECTIONS_HPP_

#include <cstdint>
#include "Vector3d.hpp"

class Directions {
    public:
        explicit Directions(std::uint32_t NumOfDirectionsLevels, bool useHEALPixGrid);

        ~Directions()
        {
            delete[] points_;
            delete[] w_;
        }

        Directions(Directions const &) = delete;
        Directions& operator=(Directions const &) = delete;

        inline std::uint64_t number() const
        {	return directionsNumber_;	}

        inline Vector3d direction(std::uint64_t index) const
        {
            return points_[index];
        }

        inline double w(std::uint64_t index) const
        {	return w_[index];	}

    private:
        Vector3d *points_{ nullptr };
        double   *w_{ nullptr };
        std::uint64_t directionsNumber_{ 0 };
};

#endif
