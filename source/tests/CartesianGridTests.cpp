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
    SECTION("Inside")
    {
        IMatterCPtr matter = std::make_shared<PolynomialMatter>(1., 0., 0., 0.);
        CartesianGrid grid(10., 20., 100., 1. / AU_Cm, 30, 20, 200, matter);

        {
            REQUIRE(grid.inside(photon(grid, { 9.9,  19.9,  99.9}, {  1, 0, 0})));
            REQUIRE(grid.inside(photon(grid, { 9.9,  19.9, -99.9}, {  0, 1, 0})));
            REQUIRE(grid.inside(photon(grid, { 9.9, -19.9,  99.9}, {  0, 0, 1})));
            REQUIRE(grid.inside(photon(grid, { 9.9, -19.9, -99.9}, { -1, 0, 0})));

            REQUIRE(grid.inside(photon(grid, { 9.9,  19.9,  99.9}, { 0, -1,  0})));
            REQUIRE(grid.inside(photon(grid, { 9.9,  19.9, -99.9}, { 0,  0, -1})));
            REQUIRE(grid.inside(photon(grid, {-9.9, -19.9,  99.9}, { 1, -1,  0})));
            REQUIRE(grid.inside(photon(grid, {-9.9, -19.9, -99.9}, { 0,  1, -1})));
        }

        {
            REQUIRE(!grid.inside(photon(grid, { 10.1, -20.1, 100.1}, { 1,  0,  0})));
            REQUIRE(!grid.inside(photon(grid, { 10.1,  11.1, 100.1}, { 0,  1,  0})));
            REQUIRE(!grid.inside(photon(grid, { 10.1,  20.1, 100.1}, { 0,  0,  1})));
            REQUIRE(!grid.inside(photon(grid, {  9.9, -20.1, 100.1}, {-1,  0,  0})));
            REQUIRE(!grid.inside(photon(grid, {  0.0,  -9.9, 100.1}, { 0, -1,  0})));
            REQUIRE(!grid.inside(photon(grid, { -9.9, -20.1, 100.1}, { 0,  0, -1})));
            REQUIRE(!grid.inside(photon(grid, {-10.1,  20.1, 100.1}, { 1, -1,  0})));
            REQUIRE(!grid.inside(photon(grid, {-10.1,   9.9, 100.1}, { 0,  1, -1})));
            REQUIRE(!grid.inside(photon(grid, {-10.1,  20.1, 100.1}, {-1,  0,  1})));

            REQUIRE(!grid.inside(photon(grid, { 10.1, -20.1,  80}, {  1, 0, 0})));
            REQUIRE(!grid.inside(photon(grid, { 10.1,  11.1, -80}, {  0, 1, 0})));
            REQUIRE(!grid.inside(photon(grid, { 10.1,  20.1,  90}, {  0, 0, 1})));
            REQUIRE(!grid.inside(photon(grid, {  9.9, -20.1, -90}, { -1, 0, 0})));

            REQUIRE(!grid.inside(photon(grid, { -9.9,  20.1, -15}, { 0,  0, -1})));
            REQUIRE(!grid.inside(photon(grid, {-10.1, -20.1,  15}, { 1, -1,  0})));
            REQUIRE(!grid.inside(photon(grid, {-10.1, -19.9, -60}, { 0,  1, -1})));
            REQUIRE(!grid.inside(photon(grid, {-10.1,  20.1,  60}, { -1,  0, 1})));

            REQUIRE(!grid.inside(photon(grid, { 10.1, -20.1, -100.1}, { 1,  0,  0})));
            REQUIRE(!grid.inside(photon(grid, { 10.1,  11.1, -100.1}, { 0,  1,  0})));
            REQUIRE(!grid.inside(photon(grid, { 10.1,  20.1, -100.1}, { 0,  0,  1})));
            REQUIRE(!grid.inside(photon(grid, {  9.9, -20.1, -100.1}, {-1,  0,  0})));
            REQUIRE(!grid.inside(photon(grid, {  0.0,  -9.9, -100.1}, { 0, -1,  0})));
            REQUIRE(!grid.inside(photon(grid, { -9.9,  20.1, -100.1}, { 0,  0, -1})));
            REQUIRE(!grid.inside(photon(grid, {-10.1, -20.1, -100.1}, { 1, -1,  0})));
            REQUIRE(!grid.inside(photon(grid, {-10.1,   9.9, -100.1}, { 0,  1, -1})));
            REQUIRE(!grid.inside(photon(grid, {-10.1,  20.1, -100.1}, {-1,  0,  1})));
        }
    }

    SECTION("Uniform dust. findOpticalDepth")
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

    SECTION("Uniform dust. movePhotonAtDistance")
    {
        IMatterCPtr matter = std::make_shared<PolynomialMatter>(1., 0., 0., 0.);
        CartesianGrid grid(100., 100., 100., 1. / AU_Cm, 201, 201, 201, matter);

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
        IMatterCPtr matter = std::make_shared<PolynomialMatter>(1., 0., 0., 0.);
        CartesianGrid grid(100., 100., 100., 1. / AU_Cm, 201, 201, 201, matter);

        Vector3d position;
        std::uint64_t cellId = grid.cellId(position);

        Photon ph1(position, cellId, Direction3d{ {1, 0, 0} }, 1, 1);
        grid.movePhotonAtDepth(ph1, 0.1, 0.0);
        REQUIRE(Approx((ph1.pos()-position).norm()) == 0.1);

        Photon ph2(position, cellId, Direction3d{ {1, 1, 1} }, 1, 1);
        grid.movePhotonAtDepth(ph2, 1, 0.0);
        REQUIRE(Approx((ph2.pos()-position).norm()) == 1);

        Photon ph3(position, cellId, Direction3d{ {0, -1, -1} }, 1, 1);
        grid.movePhotonAtDepth(ph3, 10, 0.0);
        REQUIRE(Approx((ph3.pos()-position).norm()) == 10);

        Photon ph4(position, cellId, Direction3d{ {-1, 1, 0.5} }, 1, 1);
        grid.movePhotonAtDepth(ph4, 100, 0.0);
        REQUIRE(Approx((ph4.pos()-position).norm()) == 100);
    }

    SECTION("Polynomial dust density. findOpticalDepth")
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

    SECTION("Polynomial dust density. movePhotonAtDistance")
    {
        IMatterCPtr matter = std::make_shared<PolynomialMatter>(1., -0.001, -0.003, -0.006);
        CartesianGrid grid(100., 100., 100., 1. / AU_Cm, 201, 201, 201, matter);

        Vector3d position;
        std::uint64_t cellId = grid.cellId(position);

        Photon ph1(position, cellId, Direction3d{ {1, 0, 0} }, 1, 1);
        REQUIRE(Approx(grid.movePhotonAtDistance(ph1, 0.1)).epsilon(0.0001) == 0.5 * 0.1 * (1 + 0.9999));

        Photon ph2(position, cellId, Direction3d{ {0, 1, 0} }, 1, 1);
        REQUIRE(Approx(grid.movePhotonAtDistance(ph2, 1)) == 0.5 * (1 + 0.997));

        double const sqrt2 = std::sqrt(2.0);
        Photon ph3(position, cellId, Direction3d{ {0, -1, -1} }, 1, 1);
        REQUIRE(Approx(grid.movePhotonAtDistance(ph3, 10 * sqrt2)) == 0.5 * 10 * sqrt2 * (1 + 0.91));

        Photon ph4(position, cellId, Direction3d{ {-1, 0, 1} }, 1, 1);
        REQUIRE(Approx(grid.movePhotonAtDistance(ph4, 100 * sqrt2)).epsilon(0.0001) == 0.5 * 100 * sqrt2 * (1 + 0.3));
    }

    SECTION("Polynomial dust density. movePhotonAtDepth")
    {
        IMatterCPtr matter = std::make_shared<PolynomialMatter>(1., -0.001, -0.003, -0.006);
        CartesianGrid grid(100., 100., 100., 1. / AU_Cm, 201, 201, 201, matter);

        Vector3d position;
        std::uint64_t cellId = grid.cellId(position);

        Photon ph1(position, cellId, Direction3d{ {1, 0, 0} }, 1, 1);
        grid.movePhotonAtDepth(ph1, 0.5 * 0.1 * (1 + 0.9999), 0.0);
        REQUIRE(Approx((ph1.pos()-position).norm()).epsilon(0.0001) == 0.1);

        Photon ph2(position, cellId, Direction3d{ {0, 1, 0} }, 1, 1);
        grid.movePhotonAtDepth(ph2, 0.5 * (1 + 0.997), 0.0);
        REQUIRE(Approx((ph2.pos()-position).norm()) == 1.0);

        double const sqrt2 = std::sqrt(2.0);
        Photon ph3(position, cellId, Direction3d{ {0, -1, -1} }, 1, 1);
        grid.movePhotonAtDepth(ph3, 0.5 * 10 * sqrt2 * (1 + 0.91), 0.0);
        REQUIRE(Approx((ph3.pos()-position).norm()) == 10 * sqrt2);

        Photon ph4(position, cellId, Direction3d{ {-1, 0, 1} }, 1, 1);
        grid.movePhotonAtDepth(ph4, 0.5 * 100 * sqrt2 * (1 + 0.3), 0.0);
        REQUIRE(Approx((ph4.pos()-position).norm()).epsilon(0.0001) == 100 * sqrt2);
    }
}
