#include "../third-party/catch2/catch.hpp"
#include "../WhiteDust.hpp"
#include "../Predefines.hpp"

TEST_CASE("White Dust", "[dust]")
{
    double const albedo{ 0.545 };
    double const g{ 0.48 };
    double const pl{ 0.38 };

    IDustCPtr dust = std::make_shared<WhiteDust>(albedo, g, pl, 0.0, 1.);

    SECTION("Forward scattering")
    {
        double p1, p2, p3, p4;
        dust->scatteringMatrixElements(p1, p2, p3, p4, 1.);
        REQUIRE(Approx(p1) == (1. + g) / (1. - g) / (1. - g));
        REQUIRE(Approx(p2) == 0.0);
        REQUIRE(Approx(p3) == p1);
    }

    SECTION("Backward scattering")
    {
        double p1, p2, p3, p4;
        dust->scatteringMatrixElements(p1, p2, p3, p4, -1.);
        REQUIRE(Approx(p1) == (1. - g) / (1. + g) / (1. + g));
        REQUIRE(Approx(p2) == 0.0);
        REQUIRE(Approx(p3) == -p1);
    }

    SECTION("Orthogonal scattering")
    {
        double p1, p2, p3, p4;
        dust->scatteringMatrixElements(p1, p2, p3, p4, 0);
        REQUIRE(Approx(p1) == (1. - g * g) / std::pow(1. + g * g, 1.5));
        REQUIRE(Approx(p2) == -p1 * pl);
        REQUIRE(Approx(p3) == 0.0);
    }

    SECTION("60 degrees scattering")
    {
        double p1, p2, p3, p4;
        dust->scatteringMatrixElements(p1, p2, p3, p4, 0.5);
        REQUIRE(Approx(p1) == (1. - g * g) / std::pow(1. -g + g * g, 1.5));
        REQUIRE(Approx(p2) == -p1 * pl * 3 / 5);
        REQUIRE(Approx(p3) == p1 * 0.8);
    }
}
