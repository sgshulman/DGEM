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
    class PolynomialMatter : public IMatter
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

    void createGrid(std::stringstream* nodes, std::stringstream* elements, std::size_t const n, double const r)
    {
        *nodes << (n+1)*(n+1)*(n+1) << std::endl;

        std::size_t nodeId{ 1 };
        double const dr = 2 * r / n;

        for (std::size_t i = 0; i != n+1; ++i)
        {
            double const x = -r + dr * i;
            for (std::size_t j = 0; j != n+1; ++j)
            {
                double const y = -r + dr * j;
                for (std::size_t k = 0; k != n+1; ++k)
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
}

TEST_CASE("Tetrahedral Grid", "[grid]")
{
    SECTION("Uniform dust")
    {
        std::stringstream nodes;
        std::stringstream elements;
        createGrid(&nodes, &elements, 10, 100.);

        IMatterCPtr matter = std::make_shared<PolynomialMatter>(1., 0., 0., 0.);
        TetrahedralGrid grid(nodes, elements, 100., 1. / AU_Cm, matter);

        REQUIRE(Approx(8e6 * GPerCm3_MSunPerAU3) == grid.computeMatterMass());

        Vector3d position;
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, { 1, 0, 0}))).epsilon(0.02) == 100.);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {-1, 0, 0}))).epsilon(0.02) == 100.);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {0,  1, 0}))).epsilon(0.02) == 100.);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {0, -1, 0}))).epsilon(0.02) == 100.);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {0, 0,  1}))).epsilon(0.02) == 100.);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {0, 0, -1}))).epsilon(0.02) == 100.);

        double const depth2{ std::sqrt(2.) * 100. };
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, { 1,  1, 0}))).epsilon(0.04) == depth2);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, { 1, -1, 0}))).epsilon(0.04) == depth2);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {-1,  1, 0}))).epsilon(0.04) == depth2);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {-1, -1, 0}))).epsilon(0.04) == depth2);

        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, { 1, 0,  1}))).epsilon(0.04) == depth2);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, { 1, 0, -1}))).epsilon(0.04) == depth2);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {-1, 0,  1}))).epsilon(0.04) == depth2);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {-1, 0, -1}))).epsilon(0.04) == depth2);

        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {0,  1,  1}))).epsilon(0.04) == depth2);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {0,  1, -1}))).epsilon(0.04) == depth2);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {0, -1,  1}))).epsilon(0.04) == depth2);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {0, -1, -1}))).epsilon(0.04) == depth2);

        double const depth3{ std::sqrt(3.) * 100. };
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, { 1,  1,  1}))).epsilon(0.06) == depth3);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, { 1,  1, -1}))).epsilon(0.06) == depth3);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, { 1, -1,  1}))).epsilon(0.06) == depth3);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, { 1, -1, -1}))).epsilon(0.06) == depth3);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {-1,  1,  1}))).epsilon(0.06) == depth3);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {-1,  1, -1}))).epsilon(0.06) == depth3);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {-1, -1,  1}))).epsilon(0.06) == depth3);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {-1, -1, -1}))).epsilon(0.06) == depth3);
    }

    SECTION("Polynomial dust density")
    {
        std::stringstream nodes;
        std::stringstream elements;
        createGrid(&nodes, &elements, 10, 100.);

        IMatterCPtr matter = std::make_shared<PolynomialMatter>(1., -0.001, -0.003, -0.006);
        TetrahedralGrid grid(nodes, elements, 100., 1. / AU_Cm, matter);

        REQUIRE(Approx(4e6 * GPerCm3_MSunPerAU3) == grid.computeMatterMass());

        Vector3d position;
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, { 1, 0, 0}))).epsilon(0.03) == 95.);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {-1, 0, 0}))).epsilon(0.03) == 95.);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {0,  1, 0}))).epsilon(0.03) == 85.);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {0, -1, 0}))).epsilon(0.03) == 85.);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {0, 0,  1}))).epsilon(0.03) == 70.);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {0, 0, -1}))).epsilon(0.03) == 70.);

        double const depth2xy{ std::sqrt(2.) * 80. };
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, { 1,  1, 0}))).epsilon(0.04) == depth2xy);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, { 1, -1, 0}))).epsilon(0.04) == depth2xy);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {-1,  1, 0}))).epsilon(0.04) == depth2xy);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {-1, -1, 0}))).epsilon(0.04) == depth2xy);

        double const depth2xz{ std::sqrt(2.) * 65. };
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, { 1, 0,  1}))).epsilon(0.04) == depth2xz);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, { 1, 0, -1}))).epsilon(0.04) == depth2xz);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {-1, 0,  1}))).epsilon(0.04) == depth2xz);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {-1, 0, -1}))).epsilon(0.04) == depth2xz);

        double const depth2yz{ std::sqrt(2.) * 55. };
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {0,  1,  1}))).epsilon(0.04) == depth2yz);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {0,  1, -1}))).epsilon(0.04) == depth2yz);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {0, -1,  1}))).epsilon(0.04) == depth2yz);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {0, -1, -1}))).epsilon(0.04) == depth2yz);

        double const depth3{ std::sqrt(3.) * 50. };
        REQUIRE(Approx(grid.findRealOpticalDepth(position, Vector3d{ 1,  1,  1})) == depth3);

        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, { 1,  1,  1}))).epsilon(0.04) == depth3);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, { 1,  1, -1}))).epsilon(0.04) == depth3);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, { 1, -1,  1}))).epsilon(0.04) == depth3);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, { 1, -1, -1}))).epsilon(0.04) == depth3);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {-1,  1,  1}))).epsilon(0.04) == depth3);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {-1,  1, -1}))).epsilon(0.04) == depth3);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {-1, -1,  1}))).epsilon(0.04) == depth3);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {-1, -1, -1}))).epsilon(0.04) == depth3);
    }
}
