#ifndef OBSERVERS_HPP_
#define OBSERVERS_HPP_

#include <fstream>
#include <iostream>
#include <cstdio>
#include "Photon.hpp"

inline void formatAngle(char* str, int const length, double const angleDeg)
{
    long const angleMilliDeg = std::lround(degrees(angleDeg) * 1000);
    long const degrees = angleMilliDeg / 1000 % 360;
    long const milliDegrees = angleMilliDeg % 1000;

    if (milliDegrees == 0)
    {
        snprintf(str, length, "%2.2li", degrees);
    } else {
        snprintf(str, length, "%2.2li-%3.3li", degrees, milliDegrees);
    }
}


// pictures 
class Pictures
{
    public:
        Pictures(std::uint32_t nx, std::uint32_t ny)
            : nx_{ nx }
            , ny_{ ny }
        {
            f_ = new double[ nx*ny ]();
            q_ = new double[ nx*ny ]();
            u_ = new double[ nx*ny ]();
        }

        ~Pictures()
        {
            delete[] f_;
            delete[] q_;
            delete[] u_;
        }

        Pictures (Pictures const& other) = delete;
        Pictures& operator=(Pictures const& other) = delete;

        Pictures(Pictures&& other) noexcept
            : nx_{ other.nx_ }
            , ny_{ other.ny_ }
            , f_{ other.f_ }
            , q_{ other.q_ }
            , u_{ other.u_ }
        {
            other.f_ = nullptr;
            other.q_ = nullptr;
            other.u_ = nullptr;
        }

        Pictures& operator=(Pictures&& other) noexcept
        {
            std::swap(f_, other.f_);
            std::swap(q_, other.q_);
            std::swap(u_, other.u_);
            nx_ = other.nx_;
            ny_ = other.ny_;
            return *this;
        }

        // place photon on the images
        void bin(Photon const& ph, int64_t const xl, int64_t const yl, int64_t const id)
        {
            // place weighted photon into image location
            if((id == -1 || static_cast<std::uint64_t>(id) == ph.nscat()) && (xl>=0) && (yl >= 0) && ((std::uint32_t)xl < nx_) && ((std::uint32_t)yl < ny_) )
            {
                f_[xl+yl*nx_] += ph.weight()*ph.fi();
                q_[xl+yl*nx_] += ph.weight()*ph.fq();
                u_[xl+yl*nx_] += ph.weight()*ph.fu();
            }
        }
            
        void normalize(std::uint64_t const numPhotons)
        {
            for (std::uint64_t i=0; i!=nx_*ny_; ++i)
            {
               f_[i] /= numPhotons;
               q_[i] /= numPhotons;
               u_[i] /= numPhotons;
            }
        }

        void write(double const phi, double const theta, int const key) const
        {
            int const FILENAME_LENGTH{ 40 };
            int const ANGLE_LENGTH{ 10 };

            char fname[FILENAME_LENGTH];
            char phiStr[ANGLE_LENGTH];
            char thetaStr[ANGLE_LENGTH];

            formatAngle(phiStr, ANGLE_LENGTH, phi);
            formatAngle(thetaStr, ANGLE_LENGTH, theta);

            snprintf(fname, FILENAME_LENGTH, "fimage%s_%s_%2.2i.dat", phiStr, thetaStr, key);
            std::ofstream f(fname);
            for (std::uint64_t y=0; y != ny_; ++y)
            {
                for (std::uint64_t x=0; x!=nx_; ++x)
                {
                    f << f_[x + y * nx_] << "\t";
                }
                f << "\n";
            }
            f.close();

            snprintf(fname, FILENAME_LENGTH, "qimage%s_%s_%2.2i.dat", phiStr, thetaStr, key);
            std::ofstream q(fname);
            for (std::uint64_t y=0; y != ny_; ++y)
            {
                for (std::uint64_t x=0; x!=nx_; ++x)
                {
                    q << q_[x + y * nx_] << "\t";
                }
                q << "\n";
            }
            q.close();

            snprintf(fname, FILENAME_LENGTH, "uimage%s_%s_%2.2i.dat", phiStr, thetaStr, key);
            std::ofstream u(fname);
            for (std::uint64_t y=0; y != ny_; ++y)
            {
                for (std::uint64_t x=0; x!=nx_; ++x)
                {
                    u << u_[x + y * nx_] << "\t";
                }
                u << "\n";
            }
            u.close();
        }

        void sum(std::ofstream& file)
        {
            double fsum=0, qsum=0, usum=0;
            for (std::uint64_t x=0; x!=nx_; ++x)
            {
                for (std::uint64_t y=0; y != ny_; ++y)
                {
                    fsum += f_[y + x * ny_];
                    qsum += q_[y + x * ny_];
                    usum += u_[y + x * ny_];
                }
            }
            file << "\tF= " << fsum << "\tQ= " << qsum << "\tU= " << usum
                    << "\tp= " << std::sqrt(qsum*qsum + usum*usum)/fsum
                    << "\tphi= " << 90 * std::atan2(usum, qsum)/PI;
        }

    private:
        std::uint32_t nx_, ny_;
        double *f_;
        double *q_;
        double *u_;
};


class Observer
{
public:
    Observer(double const phi, double const theta, double const rimage, std::uint32_t const Nx=200, std::uint32_t const Ny=200)
            : result_(Nx, Ny)
            , result0_(Nx, Ny)
            , result1_(Nx, Ny)
            , result2_(Nx, Ny)
            , direction_{phi, theta}
            , nx_(Nx)
            , ny_(Ny)
    {
        rimage_ = rimage;
        theta_ = theta;

        cosp_ = std::cos(phi);
        sinp_ = std::sin(phi);
    }

    void normalize(std::uint64_t const numPhotons)
    {
        result_.normalize(numPhotons);
        result0_.normalize(numPhotons);
        result1_.normalize(numPhotons);
        result2_.normalize(numPhotons);
    }

    void writeToMapFiles(bool const fWriteSingleAndDoubleScatterings)
    {
        result_.write(direction_.phi(), theta_, 0);

        if (fWriteSingleAndDoubleScatterings)
        {
            result1_.write(direction_.phi(), theta_, 1);
            result2_.write(direction_.phi(), theta_, 2);
        }
    }

    void write(std::ofstream& file)
    {
        file << "phi= " << degrees(direction_.phi()) << "\ttheta= " << degrees(theta_);
        result_.sum(file);
        result0_.sum(file);
        result1_.sum(file);
        result2_.sum(file);
        file << "\n";
    }

    void bin(Photon const& photon)
    {
        double const yimage = rimage_ + photon.pos().z() * direction_.sinTheta() -
            direction_.cosTheta() * (photon.pos().y()*sinp_ + photon.pos().x()*cosp_);

        double const ximage = rimage_ + photon.pos().y()*cosp_ - photon.pos().x()*sinp_;
        
        auto const xl = static_cast<int64_t>(nx_ * ximage / (2.0 * rimage_));
        auto const yl = static_cast<int64_t>(ny_ * yimage / (2.0 * rimage_));
        result_.bin(photon, xl, yl, -1);
        result0_.bin(photon, xl, yl, 0);
        result1_.bin(photon, xl, yl, 1);
        result2_.bin(photon, xl, yl, 2);
    }

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
