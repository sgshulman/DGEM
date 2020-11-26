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
        Photon(Vector3d const &position, std::uint64_t cellId, Direction3d const &dir, double weight, std::uint32_t nscat, double fi=1.0, double fq=0.0, double fu=0.0, double fv=0.0);
        // photon scattering
        double Scatt(IDustCRef dust, Direction3d const & dir, Random* ran);
        void Scatt( Model const &m, Directions const &dirs, IGridCRef grid, std::vector<Observer>& observers, Random* ran);
        void Stokes(IDustCRef dust, Direction3d const &dir, double calpha, bool fDir, Random* ran);

        void Move(double t, std::uint64_t cellId)
        {
            cellId_ = cellId;
            pos_ = pos_ + t * dir_.vector();
        }
        Vector3d& pos()
        {	return pos_;	}
        Vector3d const& pos() const
        {	return pos_;	}
        Direction3d const& dir() const
        {	return dir_;	}

        double fi() const
        {	return fi_;		}
        double fq() const
        {	return fq_;		}
        double fu() const
        {	return fu_;		}
        double fv() const
        {	return fv_;		}

        std::uint32_t& nscat()
        {	return nscat_;	}
        std::uint32_t nscat() const
        {	return nscat_;	}
        double& weight()
        {	return weight_;	}
        double weight() const
        {	return weight_;	}
        std::uint64_t& cellId()
        {	return cellId_;	}
        std::uint64_t cellId() const
        {	return cellId_;	}

        bool termination() const
        {   return nscat_ == std::numeric_limits<std::uint32_t>::max(); }

    private:
        Vector3d	pos_;	// outpoint
        Direction3d	dir_;	// vector of the direction
        std::uint32_t 	nscat_;
        std::uint64_t    cellId_;
        double weight_;
        double fi_, fq_, fu_, fv_;
};

#endif
