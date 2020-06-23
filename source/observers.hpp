#ifndef OBSERVERS_HPP_
#define OBSERVERS_HPP_

#include <fstream>
#include <iostream>
#include "model.hpp"
#include "Photon.hpp"

// pictures 
class Pictures
{
    public:
        Pictures(uint32_t nx, uint32_t ny)
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

        Pictures (Pictures&& other) = default;
        Pictures& operator=(Pictures&& other) = default;

        // place photon on the images
        void bin(Photon const& ph, int64_t const xl, int64_t const yl, int64_t const id)
        {
            // place weighted photon into image location
            if((id == -1 || static_cast<size_t>(id) == ph.nscat()) && (xl>=0) && (yl >= 0) && ((uint32_t)xl < nx_) && ((uint32_t)yl < ny_) )
            {
                f_[xl+yl*nx_] += ph.weight()*ph.fi();
                q_[xl+yl*nx_] += ph.weight()*ph.fq();
                u_[xl+yl*nx_] += ph.weight()*ph.fu();
            }
        }
            
        void normalize(size_t const numPhotons)
        {
            for (size_t i=0; i!=nx_*ny_; ++i)
            {
               f_[i] /= numPhotons;
               q_[i] /= numPhotons;
               u_[i] /= numPhotons;
            }
        }

        void write(double const phi, double const theta, int const key) const
        {
            char fname[30];
            int const iPhi = int(phi*180/3.1415926+0.5);
            int const iTheta = int(theta*180/3.1415926+0.5);

            sprintf(fname, "fimage%2.2i_%2.2i_%2.2i.dat", iPhi, iTheta, key);
            std::ofstream f(fname);
            for (size_t y=0; y != ny_; ++y)
            {
                for (size_t x=0; x!=nx_; ++x)
                    f << f_[x+y*nx_] << "\t";
                f << "\n";
            }
            f.close();
            sprintf(fname, "qimage%2.2i_%2.2i_%2.2i.dat", iPhi, iTheta, key);
            std::ofstream q(fname);
            for (size_t y=0; y != ny_; ++y)
            {
                for (size_t x=0; x!=nx_; ++x)
                    q << q_[x+y*nx_] << "\t";
                q << "\n";
            }
            q.close();

            sprintf(fname, "uimage%2.2i_%2.2i_%2.2i.dat", iPhi, iTheta, key);
            std::ofstream u(fname);
            for (size_t y=0; y != ny_; ++y)
            {
                for (size_t x=0; x!=nx_; ++x)
                    u << u_[x+y*nx_] << "\t";
                u << "\n";
            }
            u.close();
        }

        void sum(std::ofstream& file)
        {
            double fsum=0, qsum=0, usum=0;
            for (size_t x=0; x!=nx_; ++x)
            {
                for (size_t y=0; y != ny_; ++y)
                {
                    fsum += f_[y + x * ny_];
                    qsum += q_[y + x * ny_];
                    usum += u_[y + x * ny_];
                }
            }
            file << "\tF= " << fsum << "\tQ= " << qsum << "\tU= " << usum
                    << "\tp= " << std::sqrt(qsum*qsum + usum*usum)/fsum
                    << "\tphi= " << 90 * std::atan2(usum, qsum)/3.1415926 << "\n" ;
        }

    private:
        uint32_t nx_, ny_;
        double *f_;
        double *q_;
        double *u_;
};


class Observer
{
public:
    Observer(double phi, double theta, double rimage, uint32_t Nx=200, uint32_t Ny=200)
            : result_(Nx, Ny)
            , result0_(Nx, Ny)
            , result1_(Nx, Ny)
            , result2_(Nx, Ny)
            , nx_(Nx)
            , ny_(Ny)
    {
        rimage_ = rimage;

        phi_ = phi;
        if(phi_ > 2*3.1415926 ) phi_=phi_-2*3.1415926;
        if(phi_ < 0.0)          phi_=phi_+2*3.1415926;
        theta_ = theta;

        cosp_ = cos(phi_);
        sinp_ = sin(phi_);
        sint_ = sin(theta_);
        cost_ = cos(theta_);
    }

    void normalize(size_t numPhotons)
    {
        result_.normalize(numPhotons);
        result0_.normalize(numPhotons);
        result1_.normalize(numPhotons);
        result2_.normalize(numPhotons);
    }
    void writeToMapFiles(bool fWriteSingleAndDoubleScatterings)
    {
        result_.write(phi_, theta_, 0);

        if (fWriteSingleAndDoubleScatterings)
        {
            result1_.write(phi_, theta_, 1);
            result2_.write(phi_, theta_, 2);
        }
    }

    void write(std::ofstream & file)
    {
        file << "phi = " << (phi_*180/3.1415926) << "\ttheta = " << (theta_*180/3.1415926);
        result_.sum(file);
        result0_.sum(file);
    }

    void bin(Photon const& photon)
    {
        double yimage=rimage_+photon.pos().z()*sint_-photon.pos().y()*cost_*sinp_-photon.pos().x()*cost_*cosp_;
        double ximage=rimage_+photon.pos().y()*cosp_-photon.pos().x()*sinp_;
        int64_t xl=int(nx_ * ximage / (2.0 * rimage_));
        int64_t yl=int(ny_ * yimage / (2.0 * rimage_));
        result_.bin(photon, xl, yl, -1);
        result0_.bin(photon, xl, yl, 0);
        result1_.bin(photon, xl, yl, 1);
        result2_.bin(photon, xl, yl, 2);
    }

    double phi() const
    {	return phi_;	}
    double theta() const
    {	return theta_;	}
    Vector3d pos() const
    {	return {sin(theta_)*cos(phi_), sin(theta_)*sin(phi_), cos(theta_)};}

    Observer(Observer const& other) = delete;
    Observer& operator=(Observer const& other) = delete;

    Observer(Observer&& other) = default;
    Observer& operator=(Observer&& other) = default;

private:
    Pictures result_, result0_, result1_, result2_;
    uint32_t nx_, ny_;
    double rimage_;
    double phi_, theta_;
    double cosp_, sinp_;
    double cost_, sint_;
};

#endif
