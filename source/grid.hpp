#ifndef GRID_HPP_
#define GRID_HPP_

#include <math.h>
#include "model.hpp"

// cartezian grid
class Grid
{
    public:
        Grid() : Nx_(0), Ny_(0), Nz_(0), rhokappa_(nullptr) {};
        ~Grid()
        {
            if (rhokappa_ != nullptr) delete[] rhokappa_;
        }
        void Init(const Model &m, double R_i, double R_d, double rho_0, double h_0, double R_0,
                                    double alpha, double beta, uint32_t Nx, uint32_t Ny, uint32_t Nz );
        double PhotonSMax( Photon &ph ) const;
        double PhotonCWall( Photon &ph, double delta ) const;
        double TauFind( Photon ph, double delta=-0.001 ) const;
        int TauInt( Photon & ph, double tau, double tauold=0.0, double delta=-0.001 ) const;
        int TauInt2( Photon & ph, double delta=-0.001 ) const;
        void Peeloff( Photon ph, Observer &obs, std::shared_ptr<Dust const> const& dust) const;
    private:
        uint32_t Nx_, Ny_, Nz_;
        double *rhokappa_;
        double xmax_, ymax_, zmax_;
        double minrho_;
        Grid ( Grid const &);
        Grid & operator =( Grid const &);
};
#endif

