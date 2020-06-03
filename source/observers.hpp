#ifndef OBSERVERS_HPP_
#define OBSERVERS_HPP_

#include <fstream>
#include <iostream>
#include "model.hpp"
#include "photons.hpp"

// pictures 
class Pictures
{
    public:
        Pictures(uint32_t Nx, uint32_t Ny)
        {
            Nx_=Nx;
            Ny_=Ny;
            f_ = new double[ Nx*Ny ]();
            q_ = new double[ Nx*Ny ]();
            u_ = new double[ Nx*Ny ]();
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
        void Bin(Photon ph, int64_t xl, int64_t yl, int64_t id)
        {
            // place weighted photon into image location
            if((id == -1 || static_cast<size_t>(id) == ph.nscat()) && (xl>=0) && (yl >= 0) && ((uint32_t)xl < Nx_) && ((uint32_t)yl < Ny_) )
            {
                f_[xl+yl*Nx_] += ph.weight()*ph.fi();
                q_[xl+yl*Nx_] += ph.weight()*ph.fq();
                u_[xl+yl*Nx_] += ph.weight()*ph.fu();
            }
        }
            
        void Normalize(size_t numPhotons)
        {
            double imtot = 0.0;
            for (size_t i=0; i!=Nx_*Ny_; ++i)
            {
               f_[i]=f_[i]/numPhotons;
               q_[i]=q_[i]/numPhotons;
               u_[i]=u_[i]/numPhotons;
               imtot += f_[i];
            }
        }

        void Write(double phi, double theta, int key) const
        {
            char fname[30];
            sprintf(fname, "fimage%2.2i_%2.2i_%2.2i.dat", int(phi*180/3.1415926+0.5), int(theta*180/3.1415926+0.5), key);
            std::ofstream f(fname);
            for (size_t y=0; y!=Ny_; ++y)
            {
                for (size_t x=0; x!=Nx_; ++x)
                    f << f_[x+y*Nx_] << "\t";
                f << "\n";
            }
            f.close();
            sprintf(fname, "qimage%2.2i_%2.2i_%2.2i.dat", int(phi*180/3.1415926+0.5), int(theta*180/3.1415926+0.5), key);
            std::ofstream q(fname);
            for (size_t y=0; y!=Ny_; ++y)
            {
                for (size_t x=0; x!=Nx_; ++x)
                    q << q_[x+y*Nx_] << "\t";
                q << "\n";
            }
            q.close();

            sprintf(fname, "uimage%2.2i_%2.2i_%2.2i.dat", int(phi*180/3.1415926+0.5), int(theta*180/3.1415926+0.5), key);
            std::ofstream u(fname);
            for (size_t y=0; y!=Ny_; ++y)
            {
                for (size_t x=0; x!=Nx_; ++x)
                    u << u_[x+y*Nx_] << "\t";
                u << "\n";
            }
            u.close();
        }

        void Sum(std::ofstream& file)
        {
            double fsum=0, qsum=0, usum=0;
            for (size_t x=0; x!=Nx_; ++x)
            {
                for (size_t y=0; y!=Ny_; ++y)
                {
                    fsum += f_[y+x*Ny_];
                    qsum += q_[y+x*Ny_];
                    usum += u_[y+x*Ny_];
                }
            }
            file << "\tF= " << fsum << "\tQ= " << qsum << "\tU= " << usum << "\tp= " << sqrt(qsum*qsum+usum*usum)/fsum << "\tphi= " << 90*atan2(usum, qsum)/3.1415926 <<"\n" ;
        }

    private:
        double rimage_;
        uint32_t Nx_, Ny_;
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
            , Nx_(Nx)
            , Ny_(Ny)
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

    void Normalize(size_t numPhotons)
    {
        result_.Normalize(numPhotons);
        result0_.Normalize(numPhotons);
        result1_.Normalize(numPhotons);
        result2_.Normalize(numPhotons);
    }
    void WriteToMapFiles(bool fWriteSingleAndDoubleScatterings)
    {
        result_.Write(phi_, theta_, 0);
        if (fWriteSingleAndDoubleScatterings)
        {
            result1_.Write(phi_, theta_, 1);
            result2_.Write(phi_, theta_, 2);
        }
    }
    void Write(std::ofstream & file)
    {
        file << "phi = " << (phi_*180/3.1415926) << "\ttheta = " << (theta_*180/3.1415926);
        result_.Sum(file);
        result0_.Sum(file);
    }
    void Bin(Photon photon)
    {
        double yimage=rimage_+photon.pos().z()*sint_-photon.pos().y()*cost_*sinp_-photon.pos().x()*cost_*cosp_;
        double ximage=rimage_+photon.pos().y()*cosp_-photon.pos().x()*sinp_;
        int64_t xl=int(Nx_*ximage/(2.0*rimage_));
        int64_t yl=int(Ny_*yimage/(2.0*rimage_));
        result_.Bin(photon, xl, yl, -1);
        result0_.Bin(photon, xl, yl, 0);
        result1_.Bin(photon, xl, yl, 1);
        result2_.Bin(photon, xl, yl, 2);
    }
    double phi() const
    {	return phi_;	}
    double theta() const
    {	return theta_;	}
    Position pos() const
    {	return {sin(theta_)*cos(phi_), sin(theta_)*sin(phi_), cos(theta_)};}

    Observer(Observer const& other) = delete;
    Observer& operator=(Observer const& other) = delete;

    Observer(Observer&& other) = default;
    Observer& operator=(Observer&& other) = default;

private:
    Pictures result_, result0_, result1_, result2_;
    uint32_t Nx_, Ny_;
    double rimage_;
    double phi_, theta_;
    double cosp_, sinp_;
    double cost_, sint_;
};

#endif
