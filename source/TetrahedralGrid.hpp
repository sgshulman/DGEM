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
        std::istream& nodes_stream,
        std::istream& elements_stream,
        double max,
        double kappa,
        IMatterCPtr matter);

    TetrahedralGrid(
        std::string const& nodes_file,
        std::string const& elements_file,
        std::string const& binary_file,
        double max,
        double kappa,
        IMatterCPtr matter);

    TetrahedralGrid(
        std::string const& binary_file,
        double max,
        double kappa,
        IMatterCPtr matter);

    ~TetrahedralGrid() override;

    double findRealOpticalDepth(Vector3d const& position, Vector3d const& direction) const override;
    double findOpticalDepth(Photon ph) const override;
    double movePhotonAtDistance(Photon& ph, double distance) const override;
    int movePhotonAtDepth(Photon &ph, double tau, double tauold) const override;
    int movePhotonAtRandomDepth(Photon& ph, Random *ran) const override;
    void peeloff(Photon ph, Observer& observer, IDustCRef dust) const override;
    void peeloff(Photon ph, Observer& observer, IDustCRef dust, Vector3d const& pos1, Vector3d const& pos2) const override;
    double computeMatterMass() const override;
    double max() const override;
    std::uint64_t cellId(const Vector3d& position) const override;
    bool inside(const Photon& ph) const override;
private:
    class Node;
    class Tetrahedron;
    class PlaneTetrahedron;

    std::vector<Node> dots_;
    std::vector<Tetrahedron> elements_;
    double const max_;
    double const kappa_;
    IMatterCPtr matter_;

    double maxDistance(Photon const &ph) const;
    void readNodes(std::istream& nodes);
    void calculateNodesRhoKappa(double kappa);
    void readElements(std::istream& elements);
    void calculateElementSizes();
    void findElementsNeighbours();
    void saveGridToBinaryFile(std::string const & file);
    double rhoInDot(const Vector3d& dot, const Tetrahedron& el) const;
    inline bool inside_inner(std::uint64_t cellId) const;
};

#endif
