#include "TetrahedralGrid.hpp"
#include "DebugUtils.hpp"
#include "IMatter.hpp"
#include "observers.hpp"
#include "Photon.hpp"
#include "Random.hpp"
#include "Units.hpp"

#include <algorithm>
#include <map>
#include <iostream>

namespace
{
    // distance to plane in plane normal length units
    inline double distanceToPlane(Vector3d const &dot, Vector3d const &dot1, Vector3d const &dot2, Vector3d const &dot3)
    {
        Vector3d const norm = vectorProduct(dot2 - dot1, dot3 - dot1);
        Vector3d const v = dot1 - dot;
        return norm * v;
    }


    inline double timeToPlane(Vector3d const &dot1, Vector3d const &dot2, Vector3d const &dot3, Photon const &ph)
    {
        Vector3d const norm = vectorProduct(dot2 - dot1, dot3 - dot1);
        Vector3d const v = dot1 - ph.pos();
        double const d = norm * v;
        double const e = norm * ph.dir().vector();
        double t = std::numeric_limits<double>::max();
        if (e != 0.0) t = d / e;
        if (t <= 0.0) t = std::numeric_limits<double>::max();
        return t;
    }
}

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

    Vector3d const& pos() const
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

    Tetrahedron(
        std::uint32_t dot_1,
        std::uint32_t dot_2,
        std::uint32_t dot_3,
        std::uint32_t dot_4,
        std::uint32_t neighbor1,
        std::uint32_t neighbor2,
        std::uint32_t neighbor3,
        std::uint32_t neighbor4,
        double d1,
        double d2,
        double d3,
        double size)
        : dot1(dot_1), dot2(dot_2), dot3(dot_3), dot4(dot_4)
        , neighbor1(neighbor1), neighbor2(neighbor2), neighbor3(neighbor3), neighbor4(neighbor4)
        , d1(d1), d2(d2), d3(d3)
        , fEmpty(false), size_(size)
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


class TetrahedralGrid::PlaneTetrahedron
{
public:
    PlaneTetrahedron(
        TetrahedralGrid::Node const& node1,
        TetrahedralGrid::Node const& node2,
        TetrahedralGrid::Node const& node3,
        TetrahedralGrid::Node const& node4,
        TetrahedralGrid::Tetrahedron const& element)
        : dot1(node1), dot2(node2), dot3(node3), dot4(node4)
        , neighbor1(element.neighbor1), neighbor2(element.neighbor2), neighbor3(element.neighbor3), neighbor4(element.neighbor4)
        , d1(element.d1), d2(element.d2), d3(element.d3)
        , size_(element.size()), fEmpty_(element.fEmpty)
    {}

    std::pair<double, std::uint64_t> cellDistance(Photon& ph) const
    {
        double const delta=0.001*size_;

        // сдвиг к центру элемента (изначально должны быть на границе)
        Vector3d v = 0.25*(dot1.pos()+dot2.pos()+dot3.pos()+dot4.pos())-ph.pos();
        v = (1.0/v.norm())*v;
        ph.pos() = ph.pos() + 11*delta*v;

        // поиск расстояния
        std::pair<double, std::uint64_t> d[4];

        d[0] = std::make_pair(timeToPlane(dot2.pos(), dot3.pos(), dot4.pos(), ph), neighbor1);
        d[1] = std::make_pair(timeToPlane(dot1.pos(), dot3.pos(), dot4.pos(), ph), neighbor2);
        d[2] = std::make_pair(timeToPlane(dot1.pos(), dot2.pos(), dot4.pos(), ph), neighbor3);
        d[3] = std::make_pair(timeToPlane(dot1.pos(), dot2.pos(), dot3.pos(), ph), neighbor4);

        return *std::min_element(d, d+4);
    }


    double rhoInDot(const Vector3d& dot) const
    {
        double const L1 = distanceToPlane(dot, dot2.pos(), dot3.pos(), dot4.pos()) * d1;
        double const L2 = distanceToPlane(dot, dot1.pos(), dot3.pos(), dot4.pos()) * d2;
        double const L3 = distanceToPlane(dot, dot1.pos(), dot2.pos(), dot4.pos()) * d3;
        double const L4 = 1.0 - L1 - L2 - L3;

        return L1*dot1.rhoKappa() + L2*dot2.rhoKappa() + L3*dot3.rhoKappa() + L4*dot4.rhoKappa();
    }

    bool fEmpty() const
    { return fEmpty_; }

private:
    TetrahedralGrid::Node dot1, dot2, dot3, dot4;
    std::uint32_t neighbor1, neighbor2, neighbor3, neighbor4;
    double d1, d2, d3;
    double size_;
    bool fEmpty_;
};


double TetrahedralGrid::findRealOpticalDepth(Vector3d const& position, Vector3d const& direction) const
{
    double tau{ 0.0 };
    Vector3d pos{ position };
    Vector3d dirNormalized{ direction.normalized() };
    double const step{0.001 * std::min(1., max_)};

    while (std::abs(pos.x()) <= max_ && std::abs(pos.y()) <= max_ && std::abs(pos.z()) <= max_)
    {
        tau += step * matter_->density(pos) * kappa_ * AU_Cm;
        pos = pos + step * dirNormalized;
    }

    return tau;
}

// Nodes for the grid
void TetrahedralGrid::readNodes(std::istream& nodes)
{
    std::uint32_t NumberOfDots;

    nodes >> NumberOfDots;
    dots_.reserve(NumberOfDots);
    for (std::uint32_t cnt=0; cnt != NumberOfDots; ++cnt)
    {
        std::uint32_t id;
        double x, y, z;
        nodes >> id >> x >> y >> z;
        dots_.emplace_back(Vector3d{x,y,z});
    }
}

void TetrahedralGrid::calculateNodesRhoKappa(double const kappa)
{
    for (std::uint32_t cnt=0; cnt != dots_.size(); ++cnt)
    {
        double rhoKappa = matter_->density(dots_[cnt].pos())*kappa * AU_Cm;
        dots_[cnt].setRhoKappa(rhoKappa);
    }
}

// Elements of the grid
void TetrahedralGrid::readElements(std::istream& elements)
{
    std::uint32_t NumberOfElements;

    elements >> NumberOfElements;
    elements_.reserve(NumberOfElements);

    for (std::uint32_t cnt=0; cnt != NumberOfElements; ++cnt)
    {
        std::uint32_t dot[4];
        elements >> dot[0] >> dot[1] >> dot[2] >> dot[3];
        std::sort(dot, dot + 4);
        elements_.emplace_back(Tetrahedron(--dot[0], --dot[1], --dot[2], --dot[3]));
    }
}


void TetrahedralGrid::calculateElementSizes()
{
    for (std::uint32_t i = 0; i != elements_.size(); ++i)
    {
        double size = (dots_[elements_[i].dot1].pos() - dots_[elements_[i].dot2].pos()).norm();
        size = std::min(size, (dots_[elements_[i].dot2].pos() - dots_[elements_[i].dot3].pos()).norm());
        size = std::min(size, (dots_[elements_[i].dot1].pos() - dots_[elements_[i].dot3].pos()).norm());
        size = std::min(size, (dots_[elements_[i].dot1].pos() - dots_[elements_[i].dot4].pos()).norm());
        size = std::min(size, (dots_[elements_[i].dot2].pos() - dots_[elements_[i].dot4].pos()).norm());
        size = std::min(size, (dots_[elements_[i].dot3].pos() - dots_[elements_[i].dot4].pos()).norm());
        elements_[i].setSize(size);

        elements_[i].d1 = 1. / distanceToPlane(dots_[elements_[i].dot1].pos(), dots_[elements_[i].dot2].pos(), dots_[elements_[i].dot3].pos(), dots_[elements_[i].dot4].pos());
        elements_[i].d2 = 1. / distanceToPlane(dots_[elements_[i].dot2].pos(), dots_[elements_[i].dot1].pos(), dots_[elements_[i].dot3].pos(), dots_[elements_[i].dot4].pos());
        elements_[i].d3 = 1. / distanceToPlane(dots_[elements_[i].dot3].pos(), dots_[elements_[i].dot1].pos(), dots_[elements_[i].dot2].pos(), dots_[elements_[i].dot4].pos());

        Vector3d middle = 0.25*(dots_[elements_[i].dot1].pos()+dots_[elements_[i].dot2].pos()+dots_[elements_[i].dot3].pos()+dots_[elements_[i].dot4].pos());
        elements_[i].fEmpty = matter_->density(middle) == 0.0;
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
    for (std::uint32_t cnt = 0; cnt != elements_.size(); ++cnt)
    {
        elements_[cnt].neighbor1 = static_cast<std::uint32_t>(elements_.size());
        elements_[cnt].neighbor2 = static_cast<std::uint32_t>(elements_.size());
        elements_[cnt].neighbor3 = static_cast<std::uint32_t>(elements_.size());
        elements_[cnt].neighbor4 = static_cast<std::uint32_t>(elements_.size());
    }

    for (std::uint32_t cnt = 0; cnt != elements_.size(); ++cnt)
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
    std::istream& nodes_stream,
    std::istream& elements_stream,
    double const max,
    double const kappa,
    IMatterCPtr matter)
    : max_{ max }
    , kappa_{ kappa }
    , matter_{ std::move(matter) }
{
    readNodes(nodes_stream);
    calculateNodesRhoKappa(kappa);
    readElements(elements_stream);
    calculateElementSizes();
    findElementsNeighbours();
}


TetrahedralGrid::TetrahedralGrid(
    std::string const& nodes_file,
    std::string const& elements_file,
    std::string const& binary_file,
    double const max,
    double const kappa,
    IMatterCPtr matter)
    : max_{ max }
    , kappa_{ kappa }
    , matter_{ std::move(matter) }
{
    std::ifstream nodes(nodes_file);
    DATA_ASSERT(nodes.is_open(), nodes_file + " should exist.");
    readNodes(nodes);
    calculateNodesRhoKappa(kappa);

    std::ifstream elements(elements_file);
    DATA_ASSERT(elements.is_open(), elements_file + " should exist.");

    readElements(elements);
    calculateElementSizes();
    findElementsNeighbours();

    if (!binary_file.empty())
    {
        saveGridToBinaryFile(binary_file);
    }
}


TetrahedralGrid::TetrahedralGrid(
    const std::string &binary_file,
    double max,
    double kappa,
    IMatterCPtr matter)
    : max_{ max }
    , kappa_{ kappa }
    , matter_{ std::move(matter) }
{
    std::ifstream gridBin(binary_file, std::ios::binary);
    DATA_ASSERT(gridBin.is_open(), binary_file + " should exist.");

    std::uint32_t dotsSize{0};
    gridBin.read((char *)&dotsSize, 4);

    dots_.reserve(dotsSize);

    for (std::uint32_t i=0; i!=dotsSize; ++i)
    {
        double x, y, z;
        gridBin.read((char *)&x, sizeof(double));
        gridBin.read((char *)&y, sizeof(double));
        gridBin.read((char *)&z, sizeof(double));
        dots_.emplace_back(Vector3d{x,y,z});
    }

    std::uint32_t elementsSize{0};
    gridBin.read((char *)&elementsSize, 4);
    elements_.reserve(elementsSize);

    for (std::uint32_t i=0; i!=elementsSize; ++i)
    {
        std::uint32_t dot[4];
        std::uint32_t neighbors[4];
        double d[3];
        double size;

        gridBin.read((char *)&dot[0], 4);
        gridBin.read((char *)&dot[1], 4);
        gridBin.read((char *)&dot[2], 4);
        gridBin.read((char *)&dot[3], 4);

        gridBin.read((char *)&neighbors[0], 4);
        gridBin.read((char *)&neighbors[1], 4);
        gridBin.read((char *)&neighbors[2], 4);
        gridBin.read((char *)&neighbors[3], 4);

        gridBin.read((char *)&d[0], sizeof(double));
        gridBin.read((char *)&d[1], sizeof(double));
        gridBin.read((char *)&d[2], sizeof(double));

        gridBin.read((char *)&size, sizeof(double));

        elements_.emplace_back(dot[0], dot[1], dot[2], dot[3],
            neighbors[0], neighbors[1], neighbors[2], neighbors[3],
            d[0], d[1], d[2], size);

        Vector3d middle = 0.25*(dots_[elements_[i].dot1].pos()+dots_[elements_[i].dot2].pos()+dots_[elements_[i].dot3].pos()+dots_[elements_[i].dot4].pos());
        elements_[i].fEmpty = matter_->density({middle.x(), middle.y(), middle.z()}) == 0.0;
    }

    calculateNodesRhoKappa(kappa);
}



TetrahedralGrid::~TetrahedralGrid() = default;


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


double TetrahedralGrid::rhoInDot(const Vector3d& dot, const Tetrahedron& el) const
{
    double const L1 = distanceToPlane(dot, dots_[el.dot2].pos(), dots_[el.dot3].pos(), dots_[el.dot4].pos()) * el.d1;
    double const L2 = distanceToPlane(dot, dots_[el.dot1].pos(), dots_[el.dot3].pos(), dots_[el.dot4].pos()) * el.d2;
    double const L3 = distanceToPlane(dot, dots_[el.dot1].pos(), dots_[el.dot2].pos(), dots_[el.dot4].pos()) * el.d3;
    double const L4 = 1.0 - L1 - L2 - L3;

    return L1*dots_[el.dot1].rhoKappa() + L2*dots_[el.dot2].rhoKappa() + L3*dots_[el.dot3].rhoKappa() + L4*dots_[el.dot4].rhoKappa();
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
        const auto& rawElement = elements_[ph.cellId()];

        PlaneTetrahedron const element{
            dots_[rawElement.dot1],
            dots_[rawElement.dot2],
            dots_[rawElement.dot3],
            dots_[rawElement.dot4],
            rawElement};

        std::pair<double, std::uint64_t> const dcell = element.cellDistance(ph);
        double rho=0.0;
        if (!element.fEmpty())
        {
            rho = element.rhoInDot(ph.pos() + 0.5 * dcell.first * ph.dir().vector());
        }
        taucell=dcell.first * rho;
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
        const auto& rawElement = elements_[ph.cellId()];

        PlaneTetrahedron const element{
            dots_[rawElement.dot1],
            dots_[rawElement.dot2],
            dots_[rawElement.dot3],
            dots_[rawElement.dot4],
            rawElement};

        std::pair<double, std::uint64_t> const dcell = element.cellDistance(ph);

        double rho1=0.0, rho2=0.0;
        if (!element.fEmpty())
        {
            rho1 = element.rhoInDot(ph.pos());
            rho2 = element.rhoInDot(ph.pos() + dcell.first * ph.dir().vector());
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


int TetrahedralGrid::movePhotonAtRandomDepth(Photon &ph, Random *ran) const
{
    double const tau = -std::log(ran->Get());
    return movePhotonAtDepth(ph, tau, 0.0);
}


double TetrahedralGrid::computeMatterMass() const
{
    double m = 0;
    for(std::uint32_t id=0; id != elements_.size(); ++id)
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

            Vector3d middle = 0.25*(d1.pos()+d2.pos()+d3.pos()+d4.pos());
            m += V * matter_->density({middle.x(), middle.y(), middle.z()});
        }
    }
    return m * GPerCm3_MSunPerAU3;
}


std::uint64_t TetrahedralGrid::cellId(const Vector3d &position) const
{
    for (std::uint64_t i = 0; i != elements_.size(); ++i)
    {
        Vector3d const d1 = dots_[elements_[i].dot1].pos();
        Vector3d const d2 = dots_[elements_[i].dot2].pos();
        Vector3d const d3 = dots_[elements_[i].dot3].pos();
        Vector3d const d4 = dots_[elements_[i].dot4].pos();
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

    return static_cast<std::uint64_t>(elements_.size());
}

void TetrahedralGrid::peeloff(Photon ph, Observer &observer, const IDustCPtr &dust) const
{
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


void TetrahedralGrid::saveGridToBinaryFile(const std::string &file)
{
    std::ofstream gridBin(file, std::ios::binary);
    auto size = static_cast<uint32_t>(dots_.size());
    gridBin.write((char *)&size, 4);

    for (const auto& dot : dots_)
    {
        double const x = dot.x();
        gridBin.write((char *)&x, sizeof(double));

        double const y = dot.y();
        gridBin.write((char *)&y, sizeof(double));

        double const z = dot.z();
        gridBin.write((char *)&z, sizeof(double));
    }

    auto elementsSize = static_cast<uint32_t>(elements_.size());
    gridBin.write((char *)&elementsSize, 4);

    for (const auto& element : elements_)
    {
        gridBin.write((char *)&element.dot1, 4);
        gridBin.write((char *)&element.dot2, 4);
        gridBin.write((char *)&element.dot3, 4);
        gridBin.write((char *)&element.dot4, 4);

        gridBin.write((char *)&element.neighbor1, 4);
        gridBin.write((char *)&element.neighbor2, 4);
        gridBin.write((char *)&element.neighbor3, 4);
        gridBin.write((char *)&element.neighbor4, 4);

        gridBin.write((char *)&element.d1, sizeof(double));
        gridBin.write((char *)&element.d2, sizeof(double));
        gridBin.write((char *)&element.d3, sizeof(double));

        auto z = element.size();
        gridBin.write((char *)&z, sizeof(double));
    }
}
