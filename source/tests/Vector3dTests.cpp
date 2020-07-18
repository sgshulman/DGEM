#include "../third-party/catch2/catch.hpp"
#include "../Vector3d.hpp"
#include "../model.hpp"

TEST_CASE("Vector3d", "[vector3d]")
{
    SECTION("Union vector")
    {
        Vector3d northPole{0., 0.};
        REQUIRE(Approx(1.) == northPole.norm());
        REQUIRE(Approx(1.) == northPole.z());

        Vector3d southPole{0., PI};
        REQUIRE(Approx(1.) == southPole.norm());
        REQUIRE(Approx(-1.) == southPole.z());

        Vector3d equator{PI / 2., PI / 2};
        REQUIRE(Approx(1.) == equator.norm());
        REQUIRE(Approx(1.) == equator.y());

        Vector3d v1{5 * PI / 4., PI / 4};
        REQUIRE(Approx(1.) == v1.norm());
        REQUIRE(Approx(-0.5) == v1.x());
        REQUIRE(Approx(-0.5) == v1.y());
        REQUIRE(Approx(0.5 * std::sqrt(2.)) == v1.z());
    }
}
