#ifndef PHOTONS_HPP_
#define PHOTONS_HPP_

#include <limits>
#include "model.hpp"
#include "Vector3d.hpp"
#include "Direction3d.hpp"
#include "Directions.hpp"

// photon. It's properties and methods for the work with it
class Photon
{
    public:
        Photon(Vector3d const &position, Direction3d const &dir, double weight, uint32_t nscat, double fi=1.0, double fq=0.0, double fu=0.0, double fv=0.0);
        // photon scattering
        double Scatt(DustCRef dust, Direction3d const & dir );
        void Scatt( Model const &m, Directions const &dirs, GridCRef grid, std::vector<Observer>& observers);
        void Stokes(DustCRef dust, Direction3d const &dir, double calpha, bool fDir );
        void Move(double t)
        {
            pos_ = pos_ + t * dir_.vector();
        }
        Vector3d& pos()
        {	return pos_;	}
        Direction3d& dir()
        {	return dir_;	}
        Vector3d const& pos() const
        {	return pos_;	}
        Direction3d const& dir() const
        {	return dir_;	}

        double& x()
        {	return pos_.x();}
        double& y()
        {	return pos_.y();}
        double& z()
        {	return pos_.z();}

        double fi() const
        {	return fi_;		}
        double fq() const
        {	return fq_;		}
        double fu() const
        {	return fu_;		}
        double fv() const
        {	return fv_;		}

        uint32_t& nscat()
        {	return nscat_;	}
        uint32_t nscat() const
        {	return nscat_;	}
        double& weight()
        {	return weight_;	}
        double weight() const
        {	return weight_;	}

        bool termination() const
        {   return nscat_ == std::numeric_limits<uint32_t>::max(); }

    private:
        Vector3d	pos_;	// outpoint
        Direction3d	dir_;	// vector of the direction
        uint32_t 	nscat_;
        double weight_;
        double fi_, fq_, fu_, fv_;
};

#endif
