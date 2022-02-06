#ifndef DIRECTIONS_HPP_
#define DIRECTIONS_HPP_

#include <cstdint>
#include "Vector3d.hpp"

class Directions {
    public:
        // Directions grid for photon directions
        explicit Directions(std::uint32_t NumOfDirectionsLevels, bool useHEALPixGrid);

        // Direction grid for sphere source
        explicit Directions(std::uint32_t Nphi, std::uint32_t NumOfDirectionsLevels, bool ringScheme);

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
        // Isolatitude sphere tessalation with N_theta = 3, Nphi must be even
        void isolatitudeGrid(std::uint32_t Nphi, std::uint32_t Nside);
        void isolatitudeGridNested(std::uint32_t Nphi, std::uint32_t const Nside);

        Vector3d *points_{ nullptr };
        double   *w_{ nullptr };
        std::uint64_t directionsNumber_{ 0 };
};

#endif
