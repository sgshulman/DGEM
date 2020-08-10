#include "TetrahedralGrid.hpp"
#include "IMatter.hpp"
#include "Photon.hpp"

#include <algorithm>
#include <fstream>
#include <map>
#include <iostream>

class TetrahedralGrid::Node
{
public:
    explicit Node(Vector3d v)
        : pos_(v)
        , rhokappa_(0.0)
    {}

    void setRhoKappa(double rhokappa)
    {
        rhokappa_ = rhokappa;
    }

    double rhoKappa() const
    {
        return rhokappa_;
    }

    operator Vector3d() const // NOLINT
    {
        return pos_;
    }

    double x() const
    {
        return pos_.x();
    }

    double y() const
    {
        return pos_.y();
    }

    double z() const
    {
        return pos_.z();
    }

private:
    Vector3d pos_;
    double rhokappa_;
};


class TetrahedralGrid::Tetrahedron
{
public:
    Tetrahedron(
        std::uint32_t dot_1,
        std::uint32_t dot_2,
        std::uint32_t dot_3,
        std::uint32_t dot_4)
        : dot1(dot_1), dot2(dot_2), dot3(dot_3), dot4(dot_4)
        , neighbor1(0), neighbor2(0), neighbor3(0), neighbor4(0)
        , d1(0.0), d2(0.0), d3(0.0)
        , fEmpty(false), size_(0)
    {}

    void setSize(double size)
    {   size_ = size;   }

    double size() const
    {   return size_;        }

    std::uint32_t dot1, dot2, dot3, dot4; // must be sorted
    std::uint32_t neighbor1, neighbor2, neighbor3, neighbor4;
    double d1, d2, d3;
    bool fEmpty;
private:
    double size_;
};


double TetrahedralGrid::calculateRealTau(const Vector3d& v, double const kappa) const
{
    double tau=0.0;
    for (double x = 0.0, y = 0.0, z=0.0;
         x <= max_ && y <= max_ && z <= max_;
         x+=0.001*v.x(), y+=0.001*v.y(), z+=0.001*v.z())
    {
        tau += 0.001*matter_->density({x,y,z})*kappa*1.5e13;
    }
    return tau;
}

// Nodes for the grid
void TetrahedralGrid::readNodes(std::string const& file)
{
    uint32_t NumberOfDots;
    std::ifstream nodes(file);
    nodes >> NumberOfDots;
    dots_.reserve(NumberOfDots);
    for (uint32_t cnt=0; cnt != NumberOfDots; ++cnt)
    {
        std::uint32_t id;
        double x, y, z;
        nodes >> id >> x >> y >> z;
        dots_.emplace_back(Vector3d{x,y,z});
    }
}

void TetrahedralGrid::calculateNodesRhoKappa(double const kappa)
{
    for (uint32_t cnt=0; cnt != dots_.size(); ++cnt)
    {
        double rhoKappa = matter_->density(dots_[cnt])*kappa*1.5e13;
        dots_[cnt].setRhoKappa(rhoKappa);
    }
}

// Elements of the grid
void TetrahedralGrid::readElements(std::string const & file)
{
    uint32_t NumberOfElements;
    std::ifstream nodes(file);
    nodes >> NumberOfElements;
    elements_.reserve(NumberOfElements);

    for (uint32_t cnt=0; cnt != NumberOfElements; ++cnt)
    {
        std::uint32_t id, rubbish;
        std::uint32_t dot[4];
        nodes >> id >> rubbish >> rubbish >> rubbish >> rubbish >> dot[0] >> dot[1] >> dot[2] >> dot[3];
        std::sort(dot, dot + 4);
        elements_.emplace_back(Tetrahedron(--dot[0], --dot[1], --dot[2], --dot[3]));
    }
}


void TetrahedralGrid::calculateElementSizes()
{
    for (uint32_t i = 0; i != elements_.size(); ++i)
    {
        double size = (dots_[elements_[i].dot1] - dots_[elements_[i].dot2]).norm();
        size = std::min(size, (dots_[elements_[i].dot2] - dots_[elements_[i].dot3]).norm());
        size = std::min(size, (dots_[elements_[i].dot1] - dots_[elements_[i].dot3]).norm());
        size = std::min(size, (dots_[elements_[i].dot1] - dots_[elements_[i].dot4]).norm());
        size = std::min(size, (dots_[elements_[i].dot2] - dots_[elements_[i].dot4]).norm());
        size = std::min(size, (dots_[elements_[i].dot3] - dots_[elements_[i].dot4]).norm());
        elements_[i].setSize(size);

        elements_[i].d1 = distanceToPlane(dots_[elements_[i].dot1], dots_[elements_[i].dot2], dots_[elements_[i].dot3], dots_[elements_[i].dot4]);
        elements_[i].d2 = distanceToPlane(dots_[elements_[i].dot2], dots_[elements_[i].dot1], dots_[elements_[i].dot3], dots_[elements_[i].dot4]);
        elements_[i].d3 = distanceToPlane(dots_[elements_[i].dot3], dots_[elements_[i].dot1], dots_[elements_[i].dot2], dots_[elements_[i].dot4]);

        Vector3d middle = 0.25*(dots_[elements_[i].dot1]+dots_[elements_[i].dot2]+dots_[elements_[i].dot3]+dots_[elements_[i].dot4]);
        elements_[i].fEmpty = matter_->density({middle.x(), middle.y(), middle.z()}) == 0.0;
    }
}

class Plane
{
public:
    Plane(
        std::uint32_t dot1,
        std::uint32_t dot2,
        std::uint32_t dot3)
        : a(dot1)
        , b(dot2)
        , c(dot3)
    {}

    std::uint32_t a, b, c;
};

bool operator == (const Plane& l, const Plane& r)
{
    return l.a == r.a && l.b == r.b && l.c == r.c;
}

bool operator < (const Plane& l, const Plane& r)
{
    return l.a < r.a || (l.a == r.a && l.b < r.b) || (l.a == r.a && l.b == r.b && l.c < r.c);
}

void TetrahedralGrid::findElementsNeighbours()
{
    std::map<Plane, std::uint32_t> hash;
    for (uint32_t cnt = 0; cnt != elements_.size(); ++cnt)
    {
        elements_[cnt].neighbor1 = static_cast<uint32_t>(elements_.size());
        elements_[cnt].neighbor2 = static_cast<uint32_t>(elements_.size());
        elements_[cnt].neighbor3 = static_cast<uint32_t>(elements_.size());
        elements_[cnt].neighbor4 = static_cast<uint32_t>(elements_.size());
    }

    for (uint32_t cnt = 0; cnt != elements_.size(); ++cnt)
    {
        Plane firstPlane(elements_[cnt].dot2, elements_[cnt].dot3, elements_[cnt].dot4);
        auto iter = hash.find(firstPlane);

        if ( iter != hash.end() )
        {
            std::uint32_t cnt2 = iter->second;
            elements_[cnt].neighbor1 = cnt2;

            if (firstPlane == Plane(elements_[cnt2].dot1, elements_[cnt2].dot2, elements_[cnt2].dot3)) {
                elements_[cnt2].neighbor4 = cnt;
            } else if (firstPlane == Plane(elements_[cnt2].dot1, elements_[cnt2].dot2, elements_[cnt2].dot4)) {
                elements_[cnt2].neighbor3 = cnt;
            } else if(firstPlane == Plane(elements_[cnt2].dot1, elements_[cnt2].dot3, elements_[cnt2].dot4)) {
                elements_[cnt2].neighbor2 = cnt;
            } else {
                elements_[cnt2].neighbor1 = cnt;
            }
        } else {
            hash[firstPlane] = cnt;
        }

        Plane secondPlane(elements_[cnt].dot1, elements_[cnt].dot3, elements_[cnt].dot4);
        iter = hash.find(secondPlane);
        if ( iter != hash.end() )
        {
            std::uint32_t cnt2 = iter->second;
            elements_[cnt].neighbor2 = cnt2;

            if (secondPlane == Plane(elements_[cnt2].dot1, elements_[cnt2].dot2, elements_[cnt2].dot3)) {
                elements_[cnt2].neighbor4 = cnt;
            } else if (secondPlane == Plane(elements_[cnt2].dot1, elements_[cnt2].dot2, elements_[cnt2].dot4)) {
                elements_[cnt2].neighbor3 = cnt;
            } else if(secondPlane == Plane(elements_[cnt2].dot1, elements_[cnt2].dot3, elements_[cnt2].dot4)) {
                elements_[cnt2].neighbor2 = cnt;
            } else {
                elements_[cnt2].neighbor1 = cnt;
            }
        } else {
            hash[secondPlane] = cnt;
        }

        Plane thirdPlane(elements_[cnt].dot1, elements_[cnt].dot2, elements_[cnt].dot4);
        iter = hash.find(thirdPlane);
        if ( iter != hash.end() )
        {
            std::uint32_t cnt2 = iter->second;
            elements_[cnt].neighbor3 = cnt2;

            if (thirdPlane == Plane(elements_[cnt2].dot1, elements_[cnt2].dot2, elements_[cnt2].dot3)) {
                elements_[cnt2].neighbor4 = cnt;
            } else if (thirdPlane == Plane(elements_[cnt2].dot1, elements_[cnt2].dot2, elements_[cnt2].dot4)) {
                elements_[cnt2].neighbor3 = cnt;
            } else if(thirdPlane == Plane(elements_[cnt2].dot1, elements_[cnt2].dot3, elements_[cnt2].dot4)) {
                elements_[cnt2].neighbor2 = cnt;
            } else {
                elements_[cnt2].neighbor1 = cnt;
            }
        } else {
            hash[thirdPlane] = cnt;
        }

        Plane forthPlane(elements_[cnt].dot1, elements_[cnt].dot2, elements_[cnt].dot3);
        iter = hash.find(forthPlane);
        if ( iter != hash.end() )
        {
            std::uint32_t cnt2 = iter->second;
            elements_[cnt].neighbor4 = cnt2;

            if (forthPlane == Plane(elements_[cnt2].dot1, elements_[cnt2].dot2, elements_[cnt2].dot3)) {
                elements_[cnt2].neighbor4 = cnt;
            } else if (forthPlane == Plane(elements_[cnt2].dot1, elements_[cnt2].dot2, elements_[cnt2].dot4)) {
                elements_[cnt2].neighbor3 = cnt;
            } else if(forthPlane == Plane(elements_[cnt2].dot1, elements_[cnt2].dot3, elements_[cnt2].dot4)) {
                elements_[cnt2].neighbor2 = cnt;
            } else {
                elements_[cnt2].neighbor1 = cnt;
            }
        } else {
            hash[forthPlane] = cnt;
        }
    }
}


TetrahedralGrid::TetrahedralGrid(
    std::string const& nodes_file,
    std::string const& elements_file,
    double const max,
    double const kappa,
    IMatterCPtr matter)
    : max_{ max }
    , matter_{ std::move(matter) }
{
    readNodes(nodes_file);
    calculateNodesRhoKappa(kappa);
    readElements(elements_file);
    calculateElementSizes();
    findElementsNeighbours();
}


// calculate smax -- maximum distance photon can travel *******
double TetrahedralGrid::maxDistance(Photon const &ph) const
{
    double dsx, dsy, dsz, smax;
    if(ph.dir().x() > 0.0) {
        dsx= ( max_ - ph.pos().x() )/ph.dir().x();
    } else if(ph.dir().x() < 0.0) {
        dsx=-( ph.pos().x() + max_ )/ph.dir().x();
    } else {
        dsx=200.0*max_;
    }
    if(ph.dir().y() > 0.0) {
        dsy= ( max_ - ph.pos().y() )/ph.dir().y();
    } else if(ph.dir().y() < 0.0) {
        dsy=-( ph.pos().y() + max_ )/ph.dir().y();
    } else {
        dsy=200.0*max_;
    }
    if(ph.dir().z() > 0.0) {
        dsz= ( max_ - ph.pos().z() )/ph.dir().z();
    } else if(ph.dir().z() < 0.0) {
        dsz=-( ph.pos().z() + max_ )/ph.dir().z();
    } else {
        dsz=200.0*max_;
    }
    smax = ( dsx < dsy) ? dsx  : dsy;
    smax = (smax < dsz) ? smax : dsz;
    return smax;
}


double TetrahedralGrid::timeToPlane(Vector3d const& dot1, Vector3d const& dot2, Vector3d const& dot3, Photon const& ph) const
{
    Vector3d norm1 = vectorProduct(dot2 - dot1, dot3 - dot1);
    Vector3d v = dot1 - ph.pos();
    double d = norm1 * v;
    double e = norm1 * ph.dir().vector();
    double t = 200 * max_;
    if (e != 0.0) t = d / e;
    if (t <= 0.0) t = 200 * max_;
    return t;
}

double TetrahedralGrid::distanceToPlane(Vector3d const& dot, Vector3d const& dot1, Vector3d const& dot2, Vector3d const& dot3) const
{
    Vector3d norm1 = vectorProduct(dot2 - dot1, dot3 - dot1);
    Vector3d v = dot1 - dot;
    return (1./norm1.norm()) * norm1 * v;
}


double TetrahedralGrid::rhoInDot(const Vector3d& dot, const Tetrahedron& el) const
{
    double L1 = distanceToPlane(dot, dots_[el.dot2], dots_[el.dot3], dots_[el.dot4]) / el.d1;
    double L2 = distanceToPlane(dot, dots_[el.dot1], dots_[el.dot3], dots_[el.dot4]) / el.d2;
    double L3 = distanceToPlane(dot, dots_[el.dot1], dots_[el.dot2], dots_[el.dot4]) / el.d3;
    double L4 = 1.0 - L1 - L2 - L3;

    return L1*dots_[el.dot1].rhoKappa() + L2*dots_[el.dot2].rhoKappa() + L3*dots_[el.dot3].rhoKappa() + L4*dots_[el.dot4].rhoKappa();
}


std::pair<double, std::uint32_t> TetrahedralGrid::cellDistance(Photon& ph) const
{
    double const delta=0.001*(elements_[ph.cellId()].size());
    Vector3d const dot1 = dots_[elements_[ph.cellId()].dot1];
    Vector3d const dot2 = dots_[elements_[ph.cellId()].dot2];
    Vector3d const dot3 = dots_[elements_[ph.cellId()].dot3];
    Vector3d const dot4 = dots_[elements_[ph.cellId()].dot4];

    // сдвиг к центру элемента (изначально должны быть на границе)
    Vector3d v = 0.25*(dot1+dot2+dot3+dot4)-ph.pos();
    v = (1.0/v.norm())*v;
    ph.pos() = ph.pos() + 11*delta*v;

    // поиск расстояния
    double d1, d2, d3, d4, d;

    d1 = timeToPlane(dot2, dot3, dot4, ph);
    d2 = timeToPlane(dot1, dot3, dot4, ph);
    d3 = timeToPlane(dot1, dot2, dot4, ph);
    d4 = timeToPlane(dot1, dot2, dot3, ph);
    d = std::min(std::min(d1, d2), std::min(d3, d4));

    if (std::abs(d - d1) < 0.0000001) return std::make_pair(d, elements_[ph.cellId()].neighbor1);
    if (std::fabs(d - d2) < 0.0000001) return std::make_pair(d, elements_[ph.cellId()].neighbor2);
    if (std::fabs(d - d3) < 0.0000001) return std::make_pair(d, elements_[ph.cellId()].neighbor3);
    return std::make_pair(d, elements_[ph.cellId()].neighbor4);
}


// NEED TO BE REDONE
double TetrahedralGrid::findOpticalDepth(Photon ph) const
{
    if (ph.cellId() >= elements_.size()) return 0.0;
    double taurun=0.0, taucell, d=0.0;
    double const delta=0.001*(elements_[ph.cellId()].size());
    double smax = maxDistance(ph);

    if(smax < delta) return 0.0;
    while (ph.pos().x() < max_ && ph.pos().y() < max_ && ph.pos().z() < max_
           && -max_ < ph.pos().x() && -max_ < ph.pos().y() && -max_ < ph.pos().z()
           && ph.cellId() < elements_.size())
    {
        std::pair<double, std::uint32_t> dcell = cellDistance(ph);
        double rho1=0.0, rho2=0.0;
        if (!elements_[ph.cellId()].fEmpty)
        {
            rho1 = rhoInDot(ph.pos(), elements_[ph.cellId()]);
            rho2 = rhoInDot(ph.pos() + dcell.first * ph.dir().vector(), elements_[ph.cellId()]);
        }
        taucell=dcell.first* (rho1 + rho2)*0.5;
        taurun+=std::max(taucell, 0.0);
        ph.Move( dcell.first, dcell.second );
        d += dcell.first;
    }
    return taurun;
}

// NEED TO BE REDONE
int TetrahedralGrid::movePhotonAtDepth(Photon& ph, double tau, double tauold) const
{
    if (ph.cellId() >= elements_.size()) return 1;
    double taurun=tauold, taucell, d=0.0;
    double const delta=0.0001*(elements_[ph.cellId()].size());
    double smax = maxDistance(ph);
    if(smax < delta) return 1;
    // integrate through grid
    while ( taurun < tau && (ph.pos().x() < max_ && ph.pos().y() < max_ && ph.pos().z() < max_)
            && (-max_ < ph.pos().x() && -max_ < ph.pos().y() && -max_ < ph.pos().z())
            && ph.cellId() < elements_.size())
    {
        std::pair<double, std::uint32_t> dcell = cellDistance(ph);

        double rho1=0.0, rho2=0.0;
        if (!elements_[ph.cellId()].fEmpty)
        {
            rho1 = rhoInDot(ph.pos(), elements_[ph.cellId()]);
            rho2 = rhoInDot(ph.pos() + dcell.first * ph.dir().vector(), elements_[ph.cellId()]);
        }
        taucell=dcell.first* (rho1 + rho2)*0.5;
        if( (taurun+taucell) >= tau)
        {
            double a = (rho2 - rho1)/2/dcell.first;
            double b = rho1;
            double c = taurun-tau;

            double d1=0.5*(-b+sqrt(b*b-4*a*c))/a;
            d=d+d1;
            taurun=taurun+taucell;
            ph.Move( d1, ph.cellId() );
        } else {
            d=d+dcell.first;
            taurun=taurun+taucell;
            ph.Move( dcell.first, dcell.second );
        }
    }
    if((d>=(0.999*smax))) return 1;
    return 0;
}


double TetrahedralGrid::computeMatterMass() const
{
    double m = 0;
    for(uint32_t id=0; id != elements_.size(); ++id)
    {
        if (!elements_[id].fEmpty)
        {
            double V;

            Node const d1 = dots_[elements_[id].dot1];
            Node const d2 = dots_[elements_[id].dot2];
            Node const d3 = dots_[elements_[id].dot3];
            Node const d4 = dots_[elements_[id].dot4];

            V = d2.x() * d3.y() * d4.z() + d4.x() * d2.y() * d3.z() +
                d3.x() * d4.y() * d2.z() - d4.x() * d3.y() * d2.z() -
                d3.x() * d2.y() * d4.z() - d2.x() * d4.y() * d3.z();

            V = V - d1.x() * d3.y() * d4.z() - d4.x() * d1.y() * d3.z() -
                d3.x() * d4.y() * d1.z() + d4.x() * d3.y() * d1.z() +
                d3.x() * d1.y() * d4.z() + d1.x() * d4.y() * d3.z();

            V = V + d1.x() * d2.y() * d4.z() + d4.x() * d1.y() * d2.z() +
                d2.x() * d4.y() * d1.z() - d4.x() * d2.y() * d1.z() -
                d2.x() * d1.y() * d4.z() - d1.x() * d4.y() * d2.z();

            V = V - d1.x() * d2.y() * d3.z() - d3.x() * d1.y() * d2.z() -
                d2.x() * d3.y() * d1.z() + d3.x() * d2.y() * d1.z() +
                d2.x() * d1.y() * d3.z() + d1.x() * d3.y() * d2.z();

            V = fabs(V) / 6;

            Vector3d middle = 0.25*(d1+d2+d3+d4);
            m += V * matter_->density({middle.x(), middle.y(), middle.z()});
        }
    }
    return m*3.34792898e+6/2;
}


std::uint32_t TetrahedralGrid::cellId(const Vector3d &position) const
{
    for (size_t i = 0; i != elements_.size(); ++i)
    {
        Vector3d const d1 = dots_[elements_[i].dot1];
        Vector3d const d2 = dots_[elements_[i].dot2];
        Vector3d const d3 = dots_[elements_[i].dot3];
        Vector3d const d4 = dots_[elements_[i].dot4];
        Vector3d const center = 0.25 * (d1 + d2 + d3 + d4);

        Vector3d r1 = center - d1;
        Vector3d r2 = center - d2;
        Vector3d r3 = center - d3;
        Vector3d r4 = center - d4;
        double const size = std::max( { r1*r1, r2*r2, r3*r3, r4*r4 } );

        Vector3d r{ center - position };
        if (r * r <= size)
        {
            if (tripleProduct(position-d1, d2-d1, d3-d1)*tripleProduct(d4-d1, d2-d1, d3-d1) >= 0
                && tripleProduct(position-d1, d2-d1, d4-d1)*tripleProduct(d3-d1, d2-d1, d4-d1) >= 0
                && tripleProduct(position-d1, d3-d1, d4-d1)*tripleProduct(d2-d1, d3-d1, d4-d1) >= 0
                && tripleProduct(position-d2, d3-d2, d4-d2)*tripleProduct(d1-d2, d3-d2, d4-d2) >= 0)
            {
                return i;
            }
        }
    }

    return elements_.size();
}

