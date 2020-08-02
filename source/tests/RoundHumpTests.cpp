#include "../third-party/catch2/catch.hpp"
#include "../RoundHump.hpp"
#include "../Predefines.hpp"

TEST_CASE("Round Hump", "[matter]")
{
    SECTION("Zero Hump")
    {
        IDiskHumpCPtr const hump = std::make_shared<RoundHump const>(0., 1., 1.);

        REQUIRE(Approx(1.) == hump->hump(1., Vector3d{1., 0., 0.}));
        REQUIRE(Approx(1.) == hump->hump(1., Vector3d{1., 1., 0.}));
    }

    SECTION("Hump")
    {
        IDiskHumpCPtr const hump = std::make_shared<RoundHump const>(1., 1., 1.);

        REQUIRE(Approx(2.) == hump->hump(1., Vector3d{1., 0., 0.}));
        REQUIRE(Approx(1. + std::exp(-2.)) == hump->hump(1., Vector3d{3., 0., 0.}));
        REQUIRE(Approx(1. + std::exp(-2.)) == hump->hump(1., Vector3d{1., 2., 1.}));

        IDiskHumpCPtr const hump2 = std::make_shared<RoundHump const>(2., 1., 1.);

        REQUIRE(Approx(3.) == hump2->hump(1., Vector3d{1., 0., 0.}));
        REQUIRE(Approx(1. + 2. * std::exp(-2.)) == hump2->hump(1., Vector3d{3., 0., 0.}));
        REQUIRE(Approx(2. * (1. + 2. * std::exp(-2.))) == hump2->hump(2., Vector3d{1., -2., 1.}));

        IDiskHumpCPtr const hump3 = std::make_shared<RoundHump const>(2., 2., 2.);

        REQUIRE(Approx(3.) == hump3->hump(1., Vector3d{2., 0., 0.}));
        REQUIRE(Approx(1. + 2. * std::exp(-1.)) == hump3->hump(1., Vector3d{4., 0., 0.}));
        REQUIRE(Approx(2. * (1. + 2. * std::exp(-1.))) == hump3->hump(2., Vector3d{2., -2., 1.}));
    }
}
