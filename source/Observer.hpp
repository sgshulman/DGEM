#ifndef OBSERVERS_HPP_
#define OBSERVERS_HPP_

#include <fstream>
#include <cstdio>
#include "Photon.hpp"
#include "Vector2d.hpp"

class Observer
{
public:
    Observer(double phi, double theta, double rImage, double rMask=0.0, std::uint32_t Nx=200, std::uint32_t Ny=200);
    ~Observer();

    void normalize(std::uint64_t numPhotons);
    void writeToMapFiles(bool fWriteSingleAndDoubleScatterings, std::uint32_t numberOfScatterings);
    void write(std::ostream& file);
    bool inFov(Vector3d const& pos) const;
    void bin(Photon const& photon);
    void bin(Photon const& photon, Vector3d const& pos1, Vector3d const& pos2);
    void binHex(Photon const& photon, Vector3d const& pos1, Vector3d const& pos2);

    double phi() const
    {	return direction_.phi();	}
    double theta() const
    {	return theta_;	}
    const Direction3d& direction() const
    {	return direction_;}

    Observer(Observer const& other) = delete;
    Observer& operator=(Observer const& other) = delete;

    Observer(Observer&& other) noexcept;
    Observer& operator=(Observer&& other) = delete;

    double totalLuminosity(std::int64_t x, std::int64_t y) const
    {   return result_[3*(x + y * nx_)]; }

private:
    inline Vector2d project(Vector3d const &position) const;

    inline void binToPixel(Photon const& photon, int64_t x, int64_t y, double weight);
    inline void binPoint(Photon const& photon, int64_t x, int64_t y, bool onXBorder, bool onYBorder, double weight);
    inline void binLine(Photon const& photon, const Vector2d &pos1, const Vector2d &pos2, double lineWeight);
    inline void binMaskedLine(Photon const& photon, const Vector3d &pos1, const Vector3d &pos2, double lineWeight);

    constexpr static int NUM_OF_RESULTS{ 3 };

    double* result_{ nullptr };
    double* results_[NUM_OF_RESULTS];
    Direction3d direction_;
    std::uint32_t nx_, ny_;
    double rImage_, rImageRev_;
    Vector2d pixelSize_;
    double rMask2_;
    double theta_;
    double cosp_, sinp_;
};

#endif
