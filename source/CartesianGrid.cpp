#include <algorithm>
#include <cmath>
#include <limits>

#include "CartesianGrid.hpp"
#include "IMatter.hpp"
#include "IRandomGenerator.hpp"
#include "Observer.hpp"
#include "Photon.hpp"
#include "Sources.hpp"
#include "Units.hpp"

CartesianGrid::CartesianGrid(
        double const xmax,
        double const ymax,
        double const zmax,
        double const kappa,
        std::uint32_t const nx,
        std::uint32_t const ny,
        std::uint32_t const nz,
        IMatterCPtr matter)
    : nx_{ nx }
    , ny_{ ny }
    , nz_{ nz }
    , maxCellId_{ (static_cast<std::uint64_t>(nx_) | (static_cast<std::uint64_t>(ny_) << 20U) | (static_cast<std::uint64_t>(nz_) << 40U)) }
    , rhokappa_{ new double[nx_ * ny_ * nz_] }
    , xmax_{ xmax }
    , ymax_{ ymax }
    , zmax_{ zmax }
    , xCellSize_ { 2. * xmax_ / nx_ }
    , yCellSize_ { 2. * ymax_ / ny_ }
    , zCellSize_ { 2. * zmax_ / nz_ }
    , xCellSizeInv_{ 1. / xCellSize_ }
    , yCellSizeInv_{ 1. / yCellSize_ }
    , zCellSizeInv_{ 1. / zCellSize_ }
    , minrho_{ 1e+38 }
    , kappa_{ kappa }
    , matter_{ std::move(matter) }
{
    for (std::uint32_t cntx=0; cntx!=nx_; ++cntx)
    {
        double const x = (cntx*2.0+1) * xmax_/nx_ - xmax_;
        for (std::uint32_t cnty=0; cnty!=ny_; ++cnty)
        {
            double const y=(cnty*2.0+1) * ymax_/ny_ - ymax_;
            for (std::uint32_t cntz=0; cntz!=nz_; ++cntz)
            {
                double const z=(cntz*2.0+1) * zmax_/nz_ - zmax_;
                std::uint64_t const idx = cntx+cnty*nx+cntz*ny*nx;
                rhokappa_[idx] = matter_->density({x, y, z}) * kappa * AU_Cm; // rho*kappa*R,
                if (minrho_ > rhokappa_[idx] && rhokappa_[idx] > 0)
                {
                    minrho_ = rhokappa_[idx];
                }
            }
        }
    }
}


double CartesianGrid::findRealOpticalDepth(Vector3d const& position, Vector3d const& direction) const
{
    double tau{ 0.0 };
    Vector3d pos{ position };
    Vector3d dirNormalized{ direction.normalized() };

    double const step{0.001 * std::min({1., xmax_, ymax_, zmax_})};

    while (std::abs(pos.x()) <= xmax_ && std::abs(pos.y()) <= ymax_ && std::abs(pos.z()) <= zmax_)
    {
        tau += step * matter_->density(pos) * kappa_ * AU_Cm;
        pos = pos + step * dirNormalized;
    }

    return tau;
}


double CartesianGrid::findOpticalDepth(Photon ph) const
{
    double taurun=0.0;

    Vector3d const phDirInv = ph.dir().vector().inverse();

    Vector3d const phDirPos(
        ph.dir().x() > 0.0 ? 1. : ph.dir().x() < 0.0 ? 0. : std::numeric_limits<double>::infinity(),
        ph.dir().y() > 0.0 ? 1. : ph.dir().y() < 0.0 ? 0. : std::numeric_limits<double>::infinity(),
        ph.dir().z() > 0.0 ? 1. : ph.dir().z() < 0.0 ? 0. : std::numeric_limits<double>::infinity());

    std::int64_t const dCellX = ph.dir().x() > 0.0 ? 1 : ph.dir().x() < 0.0 ? -1 : 0;
    std::int64_t const dCellY = ph.dir().y() > 0.0 ? 0x000000000100000 : ph.dir().y() < 0.0 ? -0x000000000100000 : 0;
    std::int64_t const dCellZ = ph.dir().z() > 0.0 ? 0x000010000000000 : ph.dir().z() < 0.0 ? -0x000010000000000 : 0;

    auto x = static_cast<std::uint32_t>( ph.cellId() & 0x0000000000FFFFU);
    auto y = static_cast<std::uint32_t>((ph.cellId() & 0x00000FFFF00000U) >> 20U);
    auto z = static_cast<std::uint32_t>((ph.cellId() & 0xFFFF0000000000U) >> 40U);

    Vector3d border(
        (x + phDirPos.x()) * xCellSize_ - xmax_,
        (y + phDirPos.y()) * yCellSize_ - ymax_,
        (z + phDirPos.z()) * zCellSize_ - zmax_);

    // integrate through grid
    while (inside_inner(ph.cellId()))
    {
        double const dx = (border.x() - ph.pos().x()) * phDirInv.x();
        double const dy = (border.y() - ph.pos().y()) * phDirInv.y();
        double const dz = (border.z() - ph.pos().z()) * phDirInv.z();

        double const rhocell = rhokappa_[ x+y*nx_+z*ny_*nx_ ];
        std::uint64_t newCellId = ph.cellId();
        double dcell = dx;

        if (dx < dy && dx < dz)
        {
            newCellId += dCellX;
            x += static_cast<std::uint32_t>(dCellX);
            border.x() = (x + phDirPos.x()) * xCellSize_ - xmax_;
        } else if (dy < dz) {
            dcell = dy;
            newCellId += dCellY;
            y = static_cast<std::uint32_t>((newCellId & 0x00000FFFF00000U) >> 20U);
            border.y() = (y + phDirPos.y()) * yCellSize_ - ymax_;
        } else {
            dcell = dz;
            newCellId += dCellZ;
            z = static_cast<std::uint32_t>((newCellId & 0xFFFF0000000000U) >> 40U);
            border.z() = (z + phDirPos.z()) * zCellSize_ - zmax_;
        }

        taurun += dcell * rhocell;
        ph.Move(dcell, newCellId);
    }

    return taurun;
}


double CartesianGrid::movePhotonAtDistance(Photon &ph, double distance) const
{
    double taurun=0.0, d=0.0;

    Vector3d const phDirInv = ph.dir().vector().inverse();

    Vector3d const phDirPos(
        ph.dir().x() > 0.0 ? 1. : ph.dir().x() < 0.0 ? 0. : std::numeric_limits<double>::infinity(),
        ph.dir().y() > 0.0 ? 1. : ph.dir().y() < 0.0 ? 0. : std::numeric_limits<double>::infinity(),
        ph.dir().z() > 0.0 ? 1. : ph.dir().z() < 0.0 ? 0. : std::numeric_limits<double>::infinity());

    std::int64_t const dCellX = ph.dir().x() > 0.0 ? 1 : ph.dir().x() < 0.0 ? -1 : 0;
    std::int64_t const dCellY = ph.dir().y() > 0.0 ? 0x000000000100000 : ph.dir().y() < 0.0 ? -0x000000000100000 : 0;
    std::int64_t const dCellZ = ph.dir().z() > 0.0 ? 0x000010000000000 : ph.dir().z() < 0.0 ? -0x000010000000000 : 0;

    auto x = static_cast<std::uint32_t>( ph.cellId() & 0x0000000000FFFFU);
    auto y = static_cast<std::uint32_t>((ph.cellId() & 0x00000FFFF00000U) >> 20U);
    auto z = static_cast<std::uint32_t>((ph.cellId() & 0xFFFF0000000000U) >> 40U);

    Vector3d border(
        (x + phDirPos.x()) * xCellSize_ - xmax_,
        (y + phDirPos.y()) * yCellSize_ - ymax_,
        (z + phDirPos.z()) * zCellSize_ - zmax_);

    // integrate through grid
    while (d < distance && inside_inner(ph.cellId()))
    {
        double const dx = (border.x() - ph.pos().x()) * phDirInv.x();
        double const dy = (border.y() - ph.pos().y()) * phDirInv.y();
        double const dz = (border.z() - ph.pos().z()) * phDirInv.z();

        double const rhocell = rhokappa_[ x+y*nx_+z*ny_*nx_ ];
        std::uint64_t newCellId = ph.cellId();
        double dcell = dx;

        if (dx < dy && dx < dz)
        {
             newCellId += dCellX;
             x += static_cast<std::uint32_t>(dCellX);
             border.x() = (x + phDirPos.x()) * xCellSize_ - xmax_;
        } else if (dy < dz) {
             dcell = dy;
             newCellId += dCellY;
             y = static_cast<std::uint32_t>((newCellId & 0x00000FFFF00000U) >> 20U);
             border.y() = (y + phDirPos.y()) * yCellSize_ - ymax_;
        } else {
             dcell = dz;
             newCellId += dCellZ;
             z = static_cast<std::uint32_t>((newCellId & 0xFFFF0000000000U) >> 40U);
             border.z() = (z + phDirPos.z()) * zCellSize_ - zmax_;
        }

        if(d + dcell >= distance)
        {
            double const d1 = distance - d;
            taurun += d1 * rhocell;
            ph.Move(d1, ph.cellId());
            break;
        }
        d += dcell;
        taurun += dcell * rhocell;
        ph.Move(dcell, newCellId);
    }

    return taurun;
}


bool CartesianGrid::movePhotonAtDepth(Photon & ph, double tau, double tauold) const
{
    double taurun=tauold;

    Vector3d const phDirInv = ph.dir().vector().inverse();

    Vector3d const phDirPos(
        ph.dir().x() > 0.0 ? 1. : ph.dir().x() < 0.0 ? 0. : std::numeric_limits<double>::infinity(),
        ph.dir().y() > 0.0 ? 1. : ph.dir().y() < 0.0 ? 0. : std::numeric_limits<double>::infinity(),
        ph.dir().z() > 0.0 ? 1. : ph.dir().z() < 0.0 ? 0. : std::numeric_limits<double>::infinity());

    std::int64_t const dCellX = ph.dir().x() > 0.0 ? 1 : ph.dir().x() < 0.0 ? -1 : 0;
    std::int64_t const dCellY = ph.dir().y() > 0.0 ? 0x000000000100000 : ph.dir().y() < 0.0 ? -0x000000000100000 : 0;
    std::int64_t const dCellZ = ph.dir().z() > 0.0 ? 0x000010000000000 : ph.dir().z() < 0.0 ? -0x000010000000000 : 0;

    auto x = static_cast<std::uint32_t>( ph.cellId() & 0x0000000000FFFFU);
    auto y = static_cast<std::uint32_t>((ph.cellId() & 0x00000FFFF00000U) >> 20U);
    auto z = static_cast<std::uint32_t>((ph.cellId() & 0xFFFF0000000000U) >> 40U);

    Vector3d border(
        (x + phDirPos.x()) * xCellSize_ - xmax_,
        (y + phDirPos.y()) * yCellSize_ - ymax_,
        (z + phDirPos.z()) * zCellSize_ - zmax_);

    // integrate through grid
    while (taurun < tau && inside_inner(ph.cellId()))
    {
        double const dx = (border.x() - ph.pos().x()) * phDirInv.x();
        double const dy = (border.y() - ph.pos().y()) * phDirInv.y();
        double const dz = (border.z() - ph.pos().z()) * phDirInv.z();

        double const rhocell = rhokappa_[ x+y*nx_+z*ny_*nx_ ];
        std::uint64_t newCellId = ph.cellId();
        double dcell = dx;

        if (dx < dy && dx < dz)
        {
             newCellId += dCellX;
             x += static_cast<std::uint32_t>(dCellX);
             border.x() = (x + phDirPos.x()) * xCellSize_ - xmax_;
        } else if (dy < dz) {
             dcell = dy;
             newCellId += dCellY;
             y = static_cast<std::uint32_t>((newCellId & 0x00000FFFF00000U) >> 20U);
             border.y() = (y + phDirPos.y()) * yCellSize_ - ymax_;
        } else {
             dcell = dz;
             newCellId += dCellZ;
             z = static_cast<std::uint32_t>((newCellId & 0xFFFF0000000000U) >> 40U);
             border.z() = (z + phDirPos.z()) * zCellSize_ - zmax_;
        }

        double const taucell = dcell * rhocell;
        if(taurun + taucell >= tau)
        {
            double const d1 = (tau-taurun) / rhocell;
            ph.Move(d1, ph.cellId());
        } else {
            ph.Move(dcell, newCellId);
        }

        taurun += taucell;
    }

    return inside_inner(ph.cellId());
}


bool CartesianGrid::movePhotonAtRandomDepth(Photon &ph, IRandomGenerator *ran) const
{
    double const tau = -std::log(ran->Get());
    return movePhotonAtDepth(ph, tau, 0.0);
}


void CartesianGrid::peeloff(Photon ph, Observer& observer, IDustCRef dust) const
{
    if (!observer.inFov(ph.pos()) || (sources_ && sources_->intersectSphereSource(ph.pos(), ph.dir().vector())))
    {
        return;
    }

    double const hgfac = ph.Scatt(dust, observer.direction(), nullptr);
    double const tau = findOpticalDepth(ph);

    if (tau == 0.0)
    {
        return;
    }

    ph.weight() *= hgfac * exp(-tau);
    // Bin the photon into the image according to its position and direction of travel.
    observer.bin(ph);
}



void CartesianGrid::peeloff(Photon ph, Observer &observer, const IDustCPtr &dust, Vector3d const& pos1, Vector3d const& pos2) const
{
    if (!observer.inFov(ph.pos()) || (sources_ && sources_->intersectSphereSource(ph.pos(), ph.dir().vector())))
    {
        return;
    }

    double const hgfac = ph.Scatt(dust, observer.direction(), nullptr);
    double const tau = findOpticalDepth(ph);

    if (tau == 0.0)
    {
        return;
    }

    ph.weight() *= hgfac * exp(-tau);
    // Bin the photon into the image according to its position and direction of travel.
    observer.bin(ph, pos1, pos2);
}


void CartesianGrid::peeloffHex(Photon ph, Observer &observer, const IDustCPtr &dust, Vector3d const& pos1, Vector3d const& pos2) const
{
    if (!observer.inFov(ph.pos()) || (sources_ && sources_->intersectSphereSource(ph.pos(), ph.dir().vector())))
    {
        return;
    }

    double const hgfac = ph.Scatt(dust, observer.direction(), nullptr);
    double const tau = findOpticalDepth(ph);

    if (tau == 0.0)
    {
        return;
    }

    ph.weight() *= hgfac * exp(-tau);
    // Bin the photon into the image according to its position and direction of travel.
    observer.binHex(ph, pos1, pos2);
}


double CartesianGrid::computeMatterMass() const
{
    double density = 0.0;

    for (std::uint32_t cntx=0; cntx!=nx_; ++cntx)
    {
        double const x = (cntx*2.0+1) * xmax_/nx_ - xmax_;
        for (std::uint32_t cnty=0; cnty!=ny_; ++cnty)
        {
            double const y=(cnty*2.0+1) * ymax_/ny_ - ymax_;
            for (std::uint32_t cntz=0; cntz!=nz_; ++cntz)
            {
                double const z=(cntz*2.0+1) * zmax_/nz_ - zmax_;
                density += matter_->density({x, y, z});
            }
        }
    }

    return density * GPerCm3_MSunPerAU3 * 8. * xmax_ / nx_ * ymax_ / ny_ * zmax_ / nz_;
}


double CartesianGrid::max() const
{
    return std::max({xmax_, ymax_, zmax_});
}


std::uint64_t CartesianGrid::cellId(const Vector3d& position) const
{
    double const x = (position.x()+xmax_)*xCellSizeInv_;
    double const y = (position.y()+ymax_)*yCellSizeInv_;
    double const z = (position.z()+zmax_)*zCellSizeInv_;

    std::uint64_t const xCell = (x >= 0 ? 0x10000 : 0xFFFF) + static_cast<std::int32_t>(x);
    std::uint64_t const yCell = (y >= 0 ? 0x10000 : 0xFFFF) + static_cast<std::int32_t>(y);
    std::uint64_t const zCell = (z >= 0 ? 0x10000 : 0xFFFF) + static_cast<std::int32_t>(z);

    return (xCell | (yCell << 20U) | (zCell << 40U));
}


bool CartesianGrid::inside(const Photon& ph) const
{
    return inside_inner(ph.cellId());
}


bool CartesianGrid::inside_inner(std::uint64_t const cellId) const
{
    return ((cellId ^ (cellId - maxCellId_)) & 0x100001000010000U) == 0x100001000010000U;
}


void CartesianGrid::registerSources(SourcesCPtr sources)
{
    sources_ = std::move(sources);
}
