#ifndef TETRAHEDRAL_GRID_HPP_
#define TETRAHEDRAL_GRID_HPP_

#include "IGrid.hpp"
#include "Vector3d.hpp"
#include <vector>
#include <string>

// Tetrahedral Grid
class TetrahedralGrid : public IGrid
{
public:
    TetrahedralGrid(
        std::string const& nodes_file,
        std::string const& elements_file,
        double max,
        double kappa,
        IMatterCPtr matter);

    ~TetrahedralGrid() override;

    double calculateRealTau(const Vector3d& v, double kappa) const;

    double findOpticalDepth(Photon ph) const override;
    int movePhotonAtDepth(Photon &ph, double tau, double tauold) const override;
    int movePhotonAtRandomDepth(Photon& ph, Random *ran) const override;
    void peeloff(Photon ph, Observer& observer, DustCRef dust) const override;
    double computeMatterMass() const override;
    std::uint64_t cellId(const Vector3d& position) const override;

private:
    class Node;
    class Tetrahedron;

    std::vector<Node> dots_;
    std::vector<Tetrahedron> elements_;
    double const max_;
    IMatterCPtr matter_;

    double maxDistance(Photon const &ph) const;
    std::pair<double, std::uint64_t> cellDistance(Photon &ph) const;
    double timeToPlane(Vector3d const &dot1, Vector3d const &dot2, Vector3d const &dot3, Photon const &ph) const;
    double distanceToPlane(Vector3d const &dot, Vector3d const &dot1, Vector3d const &dot2, Vector3d const &dot3) const;
    void readNodes(std::string const & file);
    void calculateNodesRhoKappa(double kappa);
    void readElements(std::string const & file);
    void calculateElementSizes();
    void findElementsNeighbours();
    double rhoInDot(const Vector3d& dot, const Tetrahedron& el) const;
};

#endif
