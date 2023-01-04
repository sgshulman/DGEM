#include "../third-party/catch2/catch.hpp"
#include "../TetrahedralGrid.hpp"
#include "../IMatter.hpp"
#include "../Units.hpp"
#include "../Photon.hpp"

#include <sstream>
#include <iostream>

namespace
{
    // d = d0 + a * |x| + b * |y| + c * |y|
    class PolynomialMatter final: public IMatter
    {
    public:
        PolynomialMatter(double d0, double a, double b, double c)
            : d0_{d0}
            , a_{a}
            , b_{b}
            , c_{c}
        {}

        double density(Vector3d const& position) const override
        {
            return d0_ + a_ * std::abs(position.x()) + b_ * std::abs(position.y()) + c_ * std::abs(position.z());
        }

    private:
        double d0_;
        double a_;
        double b_;
        double c_;
    };

    Photon photon(TetrahedralGrid const& grid, Vector3d const &position, Vector3d const &dir)
    {
        return {position, grid.cellId(position), Direction3d(dir), 1., 0};
    }

    void createGrid(std::stringstream* nodes, std::stringstream* elements, std::uint32_t const n, double const r)
    {
        *nodes << (n+1)*(n+1)*(n+1) << std::endl;

        std::size_t nodeId{ 1 };
        double const dr = 2 * r / n;

        for (std::uint32_t i = 0; i != n+1; ++i)
        {
            double const x = -r + dr * i;
            for (std::uint32_t j = 0; j != n+1; ++j)
            {
                double const y = -r + dr * j;
                for (std::uint32_t k = 0; k != n+1; ++k)
                {
                    double const z = -r + dr * k;
                    *nodes << nodeId << "\t" << x << "\t" << y << "\t" << z << "\n";
                    ++nodeId;
                }
            }
        }

        *elements << 6*n*n*n << "\n";

        std::size_t const dx = (n+1)*(n+1);
        std::size_t const dy = (n+1);
        std::size_t const dz = 1;

        for (std::size_t i = 0; i != n; ++i)
        {
            for (std::size_t j = 0; j != n; ++j)
            {
                for (std::size_t k = 0; k != n; ++k)
                {
                    std::size_t const corner = 1 + i*dx + j*dy + k*dz;

                    *elements << corner << "\t" << corner + dz << "\t" << corner + dx << "\t" << corner + dx + dy << "\n";
                    *elements << corner + dz << "\t" << corner + dx + dz << "\t" << corner + dx + dy << "\t" << corner + dx + dy + dz << "\n";
                    *elements << corner + dz << "\t" << corner + dx + dy << "\t" << corner + dy << "\t" << corner + dy + dz << "\n";
                    *elements << corner + dz << "\t" << corner + dx + dz << "\t" << corner + dx << "\t" << corner + dx + dy << "\n";
                    *elements << corner << "\t" << corner + dz << "\t" << corner + dx + dy << "\t" << corner + dy << "\n";
                    *elements << corner + dz << "\t" << corner + dx + dy + dz << "\t" << corner + dx + dy << "\t" << corner + dy + dz << "\n";
                }
            }
        }
    }
} // namespace

TEST_CASE("Tetrahedral Grid", "[grid]")
{
    SECTION("Uniform dust. findOpticalDepth")
    {
        std::stringstream nodes;
        std::stringstream elements;
        createGrid(&nodes, &elements, 10, 100.);

        IMatterCPtr matter = std::make_shared<PolynomialMatter>(1., 0., 0., 0.);
        TetrahedralGrid grid(nodes, elements, 100., 1. / AU_Cm, matter);

        REQUIRE(Approx(8e6 * GPerCm3_MSunPerAU3) == grid.computeMatterMass());

        Vector3d position;
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, { 1, 0, 0}))) == 100.);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {-1, 0, 0}))) == 100.);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {0,  1, 0}))) == 100.);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {0, -1, 0}))) == 100.);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {0, 0,  1}))) == 100.);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {0, 0, -1}))) == 100.);

        double const depth2{ std::sqrt(2.) * 100. };
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, { 1,  1, 0}))) == depth2);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, { 1, -1, 0}))) == depth2);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {-1,  1, 0}))) == depth2);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {-1, -1, 0}))) == depth2);

        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, { 1, 0,  1}))) == depth2);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, { 1, 0, -1}))) == depth2);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {-1, 0,  1}))) == depth2);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {-1, 0, -1}))) == depth2);

        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {0,  1,  1}))) == depth2);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {0,  1, -1}))) == depth2);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {0, -1,  1}))) == depth2);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {0, -1, -1}))) == depth2);

        double const depth3{ std::sqrt(3.) * 100. };
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, { 1,  1,  1}))) == depth3);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, { 1,  1, -1}))) == depth3);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, { 1, -1,  1}))) == depth3);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, { 1, -1, -1}))) == depth3);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {-1,  1,  1}))) == depth3);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {-1,  1, -1}))) == depth3);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {-1, -1,  1}))) == depth3);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {-1, -1, -1}))) == depth3);
    }

    SECTION("Uniform dust. movePhotonAtDistance")
    {
        std::stringstream nodes;
        std::stringstream elements;
        createGrid(&nodes, &elements, 10, 100.);

        IMatterCPtr matter = std::make_shared<PolynomialMatter>(1., 0., 0., 0.);
        TetrahedralGrid grid(nodes, elements, 100., 1. / AU_Cm, matter);

        Vector3d position;
        std::uint64_t cellId = grid.cellId(position);

        Photon ph1(position, cellId, Direction3d{ {1, 0, 0} }, 1, 1);
        REQUIRE(Approx(grid.movePhotonAtDistance(ph1, 0.1)) == 0.1);

        Photon ph2(position, cellId, Direction3d{ {1, 1, 1} }, 1, 1);
        REQUIRE(Approx(grid.movePhotonAtDistance(ph2, 1)) == 1);

        Photon ph3(position, cellId, Direction3d{ {0, -1, -1} }, 1, 1);
        REQUIRE(Approx(grid.movePhotonAtDistance(ph3, 10)) == 10);

        Photon ph4(position, cellId, Direction3d{ {-1, 1, 0.5} }, 1, 1);
        REQUIRE(Approx(grid.movePhotonAtDistance(ph4, 100)) == 100);
    }

    SECTION("Uniform dust. movePhotonAtDepth")
    {
        std::stringstream nodes;
        std::stringstream elements;
        createGrid(&nodes, &elements, 10, 100.);

        IMatterCPtr matter = std::make_shared<PolynomialMatter>(1., 0., 0., 0.);
        TetrahedralGrid grid(nodes, elements, 100., 1. / AU_Cm, matter);

        Vector3d position;
        std::uint64_t cellId = grid.cellId(position);

        Photon ph1(position, cellId, Direction3d{ {1, 0, 0} }, 1, 1);
        bool const inside1 = grid.movePhotonAtDepth(ph1, 0.1, 0.0);
        REQUIRE(Approx((ph1.pos()-position).norm()) == 0.1);
        REQUIRE(inside1);

        Photon ph2(position, cellId, Direction3d{ {1, 1, 1} }, 1, 1);
        bool const inside2 = grid.movePhotonAtDepth(ph2, 1, 0.0);
        REQUIRE(Approx((ph2.pos()-position).norm()) == 1);
        REQUIRE(inside2);

        Photon ph3(position, cellId, Direction3d{ {0, -1, -1} }, 1, 1);
        bool const inside3 = grid.movePhotonAtDepth(ph3, 10, 0.0);
        REQUIRE(Approx((ph3.pos()-position).norm()) == 10);
        REQUIRE(inside3);

        Photon ph4(position, cellId, Direction3d{ {-1, 1, 0.5} }, 1, 1);
        bool const inside4 = grid.movePhotonAtDepth(ph4, 100, 0.0);
        REQUIRE(Approx((ph4.pos()-position).norm()) == 100);
        REQUIRE(inside4);

        Photon ph5(position, cellId, Direction3d{ {1, 0, 0} }, 1, 1);
        bool const inside5 = grid.movePhotonAtDepth(ph5, 300, 0.0);
        REQUIRE(Approx((ph5.pos()-position).norm()) == 100);
        REQUIRE(!inside5);
    }

    SECTION("Polynomial dust density. findOpticalDepth")
    {
        std::stringstream nodes;
        std::stringstream elements;
        createGrid(&nodes, &elements, 10, 100.);

        IMatterCPtr matter = std::make_shared<PolynomialMatter>(1., -0.001, -0.003, -0.006);
        TetrahedralGrid grid(nodes, elements, 100., 1. / AU_Cm, matter);

        REQUIRE(Approx(4e6 * GPerCm3_MSunPerAU3) == grid.computeMatterMass());

        Vector3d position;
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, { 1, 0, 0}))) == 95.);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {-1, 0, 0}))) == 95.);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {0,  1, 0}))) == 85.);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {0, -1, 0}))) == 85.);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {0, 0,  1}))) == 70.);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {0, 0, -1}))) == 70.);

        double const depth2xy{ std::sqrt(2.) * 80. };
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, { 1,  1, 0}))) == depth2xy);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, { 1, -1, 0}))) == depth2xy);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {-1,  1, 0}))) == depth2xy);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {-1, -1, 0}))) == depth2xy);

        double const depth2xz{ std::sqrt(2.) * 65. };
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, { 1, 0,  1}))) == depth2xz);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, { 1, 0, -1}))) == depth2xz);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {-1, 0,  1}))) == depth2xz);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {-1, 0, -1}))) == depth2xz);

        double const depth2yz{ std::sqrt(2.) * 55. };
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {0,  1,  1}))) == depth2yz);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {0,  1, -1}))) == depth2yz);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {0, -1,  1}))) == depth2yz);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {0, -1, -1}))) == depth2yz);

        double const depth3{ std::sqrt(3.) * 50. };
        REQUIRE(Approx(grid.findRealOpticalDepth(position, Vector3d{ 1,  1,  1})) == depth3);

        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, { 1,  1,  1}))) == depth3);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, { 1,  1, -1}))) == depth3);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, { 1, -1,  1}))) == depth3);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, { 1, -1, -1}))) == depth3);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {-1,  1,  1}))) == depth3);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {-1,  1, -1}))) == depth3);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {-1, -1,  1}))) == depth3);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {-1, -1, -1}))) == depth3);
    }

    SECTION("Polynomial dust density. movePhotonAtDistance")
    {
        std::stringstream nodes;
        std::stringstream elements;
        createGrid(&nodes, &elements, 10, 100.);

        IMatterCPtr matter = std::make_shared<PolynomialMatter>(1., -0.001, -0.003, -0.006);
        TetrahedralGrid grid(nodes, elements, 100., 1. / AU_Cm, matter);

        Vector3d position;
        std::uint64_t cellId = grid.cellId(position);

        Photon ph1(position, cellId, Direction3d{ {1, 0, 0} }, 1, 1);
        REQUIRE(Approx(grid.movePhotonAtDistance(ph1, 0.1)) == 0.5 * 0.1 * (1 + 0.9999));

        Photon ph2(position, cellId, Direction3d{ {0, 1, 0} }, 1, 1);
        REQUIRE(Approx(grid.movePhotonAtDistance(ph2, 1)) == 0.5 * (1 + 0.997));

        double const sqrt2 = std::sqrt(2.0);
        Photon ph3(position, cellId, Direction3d{ {0, -1, -1} }, 1, 1);
        REQUIRE(Approx(grid.movePhotonAtDistance(ph3, 10 * sqrt2)) == 0.5 * 10 * sqrt2 * (1 + 0.91));

        Photon ph4(position, cellId, Direction3d{ {-1, 0, 1} }, 1, 1);
        REQUIRE(Approx(grid.movePhotonAtDistance(ph4, 100 * sqrt2)) == 0.5 * 100 * sqrt2 * (1 + 0.3));
    }

    SECTION("Polynomial dust density. movePhotonAtDistance")
    {
        std::stringstream nodes;
        std::stringstream elements;
        createGrid(&nodes, &elements, 10, 100.);

        IMatterCPtr matter = std::make_shared<PolynomialMatter>(1., -0.001, -0.003, -0.006);
        TetrahedralGrid grid(nodes, elements, 100., 1. / AU_Cm, matter);

        Vector3d position;
        std::uint64_t cellId = grid.cellId(position);

        Photon ph1(position, cellId, Direction3d{ {1, 0, 0} }, 1, 1);
        grid.movePhotonAtDepth(ph1, 0.5 * 0.1 * (1 + 0.9999), 0.0);
        REQUIRE(Approx((ph1.pos()-position).norm()) == 0.1);

        Photon ph2(position, cellId, Direction3d{ {0, 1, 0} }, 1, 1);
        grid.movePhotonAtDepth(ph2, 0.5 * (1 + 0.997), 0.0);
        REQUIRE(Approx((ph2.pos()-position).norm()) == 1.0);

        double const sqrt2 = std::sqrt(2.0);
        Photon ph3(position, cellId, Direction3d{ {0, -1, -1} }, 1, 1);
        grid.movePhotonAtDepth(ph3, 0.5 * 10 * sqrt2 * (1 + 0.91), 0.0);
        REQUIRE(Approx((ph3.pos()-position).norm()) == 10 * sqrt2);

        Photon ph4(position, cellId, Direction3d{ {-1, 0, 1} }, 1, 1);
        grid.movePhotonAtDepth(ph4, 0.5 * 100 * sqrt2 * (1 + 0.3), 0.0);
        REQUIRE(Approx((ph4.pos()-position).norm()) == 100 * sqrt2);
    }
}
