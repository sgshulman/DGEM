#ifndef OBSERVERS_HPP_
#define OBSERVERS_HPP_

#include <fstream>
#include <cstdio>
#include "Photon.hpp"

// pictures 
class Pictures
{
    public:
        Pictures(std::uint32_t nx, std::uint32_t ny);
        ~Pictures();

        Pictures (Pictures const& other) = delete;
        Pictures& operator=(Pictures const& other) = delete;

        Pictures(Pictures&& other) noexcept;
        Pictures& operator=(Pictures&& other) noexcept;

        // place photon on the images
        void bin(Photon const& ph, int64_t xl, int64_t yl, int64_t id);
        void normalize(std::uint64_t numPhotons);

        void write(double phi, double theta, int key) const;
        void sum(std::ofstream& file);

    private:
        std::uint32_t nx_, ny_;
        double *f_;
        double *q_;
        double *u_;
};


class Observer
{
public:
    Observer(double phi, double theta, double rimage, std::uint32_t Nx=200, std::uint32_t Ny=200);

    void normalize(std::uint64_t numPhotons);
    void writeToMapFiles(bool fWriteSingleAndDoubleScatterings);
    void write(std::ofstream& file);
    void bin(Photon const& photon);

    double phi() const
    {	return direction_.phi();	}
    double theta() const
    {	return theta_;	}
    const Direction3d& direction() const
    {	return direction_;}

    Observer(Observer const& other) = delete;
    Observer& operator=(Observer const& other) = delete;

    Observer(Observer&& other) = default;
    Observer& operator=(Observer&& other) = default;

private:
    Pictures result_, result0_, result1_, result2_;
    Direction3d direction_;
    std::uint32_t nx_, ny_;
    double rimage_;
    double theta_;
    double cosp_, sinp_;
};

#endif
