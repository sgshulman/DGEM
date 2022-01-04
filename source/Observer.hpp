#ifndef OBSERVERS_HPP_
#define OBSERVERS_HPP_

#include <ostream>
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
        void bin(Photon const& ph, int64_t xl, int64_t yl, int64_t id, double weight);
        void normalize(std::uint64_t numPhotons);

        void write(double phi, double theta, int key) const;
        void sum(std::ostream& stream);

        double f(int64_t xl, int64_t yl) const
        {   return f_[xl+yl*nx_]; }

    private:
        std::uint32_t nx_, ny_;
        double *f_;
        double *q_;
        double *u_;
};


class Observer
{
public:
    Observer(double phi, double theta, double rImage, double rMask=0.0, std::uint32_t Nx=200, std::uint32_t Ny=200);

    void normalize(std::uint64_t numPhotons);
    void writeToMapFiles(bool fWriteSingleAndDoubleScatterings, std::uint32_t numberOfScatterings);
    void write(std::ostream& file);
    bool inFov(Photon const& photon) const;
    void bin(Photon const& photon);
    void bin(Photon const& photon, const Vector3d &pos1, const Vector3d &pos2);

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

    double totalLuminosity(int64_t xl, int64_t yl) const
    {   return result_.f(xl, yl); }

private:
    inline double imageX(Vector3d const &position) const;
    inline double imageY(Vector3d const &position) const;
    inline void bin(Photon const& photon, int64_t x, int64_t y, double weight);
    inline void bin(Photon const& photon, int64_t x, int64_t y, bool onXBorder, bool onYBorder, double weight);

    Pictures result_, result0_, result1_, result2_;
    Direction3d direction_;
    std::uint32_t const nx_, ny_;
    double const rImage_;
    double const rMask2_;
    double const theta_;
    double const cosp_, sinp_;
};

#endif
