#include "../third-party/catch2/catch.hpp"
#include "../CartesianGrid.hpp"
#include "../IMatter.hpp"
#include "../Units.hpp"
#include "../Photon.hpp"

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

    Photon photon(CartesianGrid const& grid, Vector3d const &position, Vector3d const &dir)
    {
        return {position, grid.cellId(position), Direction3d(dir), 1., 0};
    }
}

TEST_CASE("Cartezian Grid", "[grid]")
{
    SECTION("Uniform dust")
    {
        IMatterCPtr matter = std::make_shared<PolynomialMatter>(1., 0., 0., 0.);
        CartesianGrid grid(100., 100., 100., 1. / AU_Cm, 201, 201, 201, matter);

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

    SECTION("Polynomial dust density")
    {
        IMatterCPtr matter = std::make_shared<PolynomialMatter>(1., -0.001, -0.003, -0.006);
        CartesianGrid grid(100., 100., 100., 1. / AU_Cm, 201, 201, 201, matter);

        REQUIRE(Approx(4e6 * GPerCm3_MSunPerAU3).epsilon(0.0001) == grid.computeMatterMass());

        Vector3d position;
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, { 1, 0, 0}))) == 95.);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {-1, 0, 0}))) == 95.);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {0,  1, 0}))) == 85.);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {0, -1, 0}))) == 85.);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {0, 0,  1}))) == 70.);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {0, 0, -1}))) == 70.);

        double const depth2xy{ std::sqrt(2.) * 80. };
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, { 1,  1, 0}))).epsilon(0.0001) == depth2xy);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, { 1, -1, 0}))).epsilon(0.0001) == depth2xy);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {-1,  1, 0}))).epsilon(0.0001) == depth2xy);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {-1, -1, 0}))).epsilon(0.0001) == depth2xy);

        double const depth2xz{ std::sqrt(2.) * 65. };
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, { 1, 0,  1}))).epsilon(0.0001) == depth2xz);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, { 1, 0, -1}))).epsilon(0.0001) == depth2xz);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {-1, 0,  1}))).epsilon(0.0001) == depth2xz);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {-1, 0, -1}))).epsilon(0.0001) == depth2xz);

        double const depth2yz{ std::sqrt(2.) * 55. };
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {0,  1,  1}))).epsilon(0.0001) == depth2yz);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {0,  1, -1}))).epsilon(0.0001) == depth2yz);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {0, -1,  1}))).epsilon(0.0001) == depth2yz);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {0, -1, -1}))).epsilon(0.0001) == depth2yz);

        double const depth3{ std::sqrt(3.) * 50. };
        REQUIRE(Approx(grid.findRealOpticalDepth(position, Vector3d{ 1,  1,  1})) == depth3);

        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, { 1,  1,  1}))).epsilon(0.0001) == depth3);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, { 1,  1, -1}))).epsilon(0.0001) == depth3);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, { 1, -1,  1}))).epsilon(0.0001) == depth3);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, { 1, -1, -1}))).epsilon(0.0001) == depth3);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {-1,  1,  1}))).epsilon(0.0001) == depth3);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {-1,  1, -1}))).epsilon(0.0001) == depth3);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {-1, -1,  1}))).epsilon(0.0001) == depth3);
        REQUIRE(Approx(grid.findOpticalDepth(photon(grid, position, {-1, -1, -1}))).epsilon(0.0001) == depth3);
    }
}
