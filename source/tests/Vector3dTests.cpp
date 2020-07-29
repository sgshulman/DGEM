#include "../third-party/catch2/catch.hpp"
#include "../Vector3d.hpp"
#include "../MathUtils.hpp"
#include "TestUtils.hpp"

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

    SECTION("Setters")
    {
        Vector3d v{1., 0., 0.};
        REQUIRE(Approx(1.) == v.norm());

        v.x() = 3;
        REQUIRE(Approx(3.) == v.norm());

        v.y() = 4;
        REQUIRE(Approx(5.) == v.norm());

        v.z() = 5;
        REQUIRE(Approx(5. * std::sqrt(2.)) == v.norm());
    }

    SECTION("Normalized")
    {
        Vector3d v{3., 0., 4.};
        REQUIRE(Approx(5.) == v.norm());

        Vector3d normalized{v.normalized()};
        REQUIRE(Approx(1.) == normalized.norm());
        REQUIRE(Approx(0.6) == normalized.x());
        REQUIRE(Approx(0.8) == normalized.z());
    }

    SECTION("Operators")
    {
        REQUIRE(equal(Vector3d{1., 3., 5.} + Vector3d{2., -3., -6}, {3., 0., -1}));
        REQUIRE(equal(Vector3d{1., 3., 5.} - Vector3d{2., -3., -6}, {-1., 6., 11}));

        REQUIRE(equal(0. * Vector3d{1., 3., 5.}, {0., 0., 0.}));
        REQUIRE(equal(1. * Vector3d{1., 3., 5.}, {1., 3., 5.}));
        REQUIRE(equal(-2. * Vector3d{1., 3., 5.}, {-2., -6., -10.}));
        REQUIRE(equal(2. * Vector3d{0., 0., 0.}, {0., 0., 0.}));

        REQUIRE(equal(Vector3d{1., 3., 5.} * 0., {0., 0., 0.}));
        REQUIRE(equal(Vector3d{1., 3., 5.} * 1., {1., 3., 5.}));
        REQUIRE(equal(Vector3d{1., 3., 5.} * -2., {-2., -6., -10.}));
        REQUIRE(equal(Vector3d{0., 0., 0.} * 2., {0., 0., 0.}));

        REQUIRE(equal(Vector3d{1., 3., 5.} / 1., {1., 3., 5.}));
        REQUIRE(equal(Vector3d{1., 3., 5.} / -0.5, {-2., -6., -10.}));
        REQUIRE(equal(Vector3d{0., 0., 0.} / 0.5, {0., 0., 0.}));

        REQUIRE(Approx(Vector3d{1., 3., 0.} * Vector3d{0., 0., 6.}).margin(1e-12) == 0.);
        REQUIRE(Approx(Vector3d{3., 4., 0.} * Vector3d{3., 4., 0.}) == 25.);
        REQUIRE(Approx(Vector3d{1., 3., -5.} * Vector3d{2., 4., 6.}) == -16.);
    }

    SECTION("Vector Product")
    {
        REQUIRE(equal(vectorProduct({1., 3., 5.}, {2., 6., 10.}), {0., 0., 0.}));
        REQUIRE(equal(vectorProduct({1., 0., 0.}, {0., 1., 0.}), {0., 0., 1.}));

        REQUIRE(Approx(vectorProduct({1., 0., 0.}, {0., 3., 4.}).norm()) == 5.);
        REQUIRE(Approx(vectorProduct({0., 1., 0.}, {-3., 0., 4.}).norm()) == 5.);
        REQUIRE(Approx(vectorProduct({0., 0., 1.}, {-3., -4., 0.}).norm()) == 5.);
    }

    SECTION("Triple Product")
    {
        REQUIRE(Approx(tripleProduct({1., 3., 0.}, {7., -8., 0.}, {5., 0., 0.})) == 0.);
        REQUIRE(Approx(tripleProduct({1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.})) == 1.);
        REQUIRE(Approx(tripleProduct({0., 2., 0.}, {-1., 0., 0.}, {0., 0., 1.})) == 2.);
        REQUIRE(Approx(tripleProduct({1., 0., 0.}, {0., 1., 0.}, {0., 0., -1.})) == -1.);
    }
}
