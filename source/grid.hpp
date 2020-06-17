#ifndef GRID_HPP_
#define GRID_HPP_

#include <math.h>
#include <cstdint>
#include <memory>

#include "Predefines.hpp"

class FlaredDisk
{
public:
    FlaredDisk(
            double rInner,
            double rOuter,
            double rho0,
            double h0,
            double r0,
            double alpha,
            double beta)
        : rInner_{ rInner }
        , rOuter_{ rOuter }
        , rho0_{ rho0 }
        , h0_{ h0 }
        , r0_{ r0 }
        , alpha_{ alpha }
        , beta_{ beta }
    {}

    double density(double x, double y, double z) const;

private:
    double const rInner_;
    double const rOuter_;
    double const rho0_;
    double const h0_;
    double const r0_;
    double const alpha_;
    double const beta_;
};

// cartezian grid
class Grid
{
    public:
        Grid(double xmax, double ymax, double zmax, double kappa,
            uint32_t Nx, uint32_t Ny, uint32_t Nz, FlaredDiskCPtr disk);

        ~Grid()
        {
            delete[] rhokappa_;
        }

        Grid(Grid const &) = delete;
        Grid& operator=(Grid const &) = delete;

        double PhotonSMax( Photon &ph ) const;
        double PhotonCWall( Photon &ph, double delta ) const;
        double TauFind( Photon ph, double delta=-0.001 ) const;
        int TauInt( Photon & ph, double tau, double tauold=0.0, double delta=-0.001 ) const;
        int TauInt2( Photon & ph, double delta=-0.001 ) const;
        void Peeloff( Photon ph, Observer &obs, DustCRef dust) const;
    private:
        uint32_t Nx_, Ny_, Nz_;
        double *rhokappa_{ nullptr };
        double xmax_, ymax_, zmax_;
        double minrho_;
        FlaredDiskCPtr disk_;
};
#endif

