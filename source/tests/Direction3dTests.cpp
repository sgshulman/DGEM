#include "../third-party/catch2/catch.hpp"
#include "../Direction3d.hpp"

TEST_CASE("Direction3d", "[vector3d]")
{
    SECTION("Vector constructor")
    {
        Direction3d northPole{Vector3d{0., 0., 2.}};
        REQUIRE(Approx(northPole.vector().norm()) == 1.);
        REQUIRE(Approx(northPole.z()) == 1.);
        REQUIRE(Approx(northPole.sinTheta()).margin(1e-12) == 0.);
        REQUIRE(Approx(northPole.cosTheta()) == 1.);

        Direction3d southPole{Vector3d{0., 0., -2.}};
        REQUIRE(Approx(southPole.vector().norm()) == 1.);
        REQUIRE(Approx(southPole.z()) == -1.);
        REQUIRE(Approx(southPole.sinTheta()).margin(1e-12) == 0.);
        REQUIRE(Approx(southPole.cosTheta()) == -1.);

        Direction3d north{Vector3d{1., 0., 0.}};
        REQUIRE(Approx(north.vector().norm()) == 1.);
        REQUIRE(Approx(north.x()) == 1.);
        REQUIRE(Approx(north.phi()) == 0.);
        REQUIRE(Approx(north.sinTheta()) == 1.);
        REQUIRE(Approx(north.cosTheta()).margin(1e-12) == 0.);

        Direction3d east{Vector3d{0., 1., 0.}};
        REQUIRE(Approx(east.vector().norm()) == 1.);
        REQUIRE(Approx(east.y()) == 1.);
        REQUIRE(Approx(east.phi()) == PI / 2.);
        REQUIRE(Approx(east.sinTheta()) == 1.);
        REQUIRE(Approx(east.cosTheta()).margin(1e-12) == 0.);

        Direction3d south{Vector3d{-1., 0., 0.}};
        REQUIRE(Approx(south.vector().norm()) == 1.);
        REQUIRE(Approx(south.x()) == -1.);
        REQUIRE(Approx(south.phi()) == PI);
        REQUIRE(Approx(south.sinTheta()) == 1.);
        REQUIRE(Approx(south.cosTheta()).margin(1e-12) == 0.);

        Direction3d northEast{Vector3d{1., 1., 0.}};
        REQUIRE(Approx(northEast.vector().norm()) == 1.);
        REQUIRE(Approx(northEast.x()) == 0.5 * std::sqrt(2.));
        REQUIRE(Approx(northEast.phi()) == PI / 4.);
        REQUIRE(Approx(northEast.sinTheta()) == 1.);
        REQUIRE(Approx(northEast.cosTheta()).margin(1e-12) == 0.);
    }

    SECTION("Coordinates constructor")
    {
        Direction3d northPole{0., 0.};
        REQUIRE(Approx(northPole.vector().norm()) == 1.);
        REQUIRE(Approx(northPole.z()) == 1.);
        REQUIRE(Approx(northPole.sinTheta()).margin(1e-12) == 0.);
        REQUIRE(Approx(northPole.cosTheta()) == 1.);

        Direction3d southPole{0., 5. * PI};
        REQUIRE(Approx(southPole.vector().norm()) == 1.);
        REQUIRE(Approx(southPole.z()) == -1.);
        REQUIRE(Approx(southPole.sinTheta()).margin(1e-12) == 0.);
        REQUIRE(Approx(southPole.cosTheta()) == -1.);

        Direction3d north{0., PI / 2.};
        REQUIRE(Approx(north.vector().norm()) == 1.);
        REQUIRE(Approx(north.x()) == 1.);
        REQUIRE(Approx(north.phi()) == 0.);
        REQUIRE(Approx(north.sinTheta()) == 1.);
        REQUIRE(Approx(north.cosTheta()).margin(1e-12) == 0.);

        Direction3d east{PI / 2., PI / 2.};
        REQUIRE(Approx(east.vector().norm()) == 1.);
        REQUIRE(Approx(east.y()) == 1.);
        REQUIRE(Approx(east.phi()) == PI / 2.);
        REQUIRE(Approx(east.sinTheta()) == 1.);
        REQUIRE(Approx(east.cosTheta()).margin(1e-12) == 0.);

        Direction3d south{3. * PI, PI / 2.};
        REQUIRE(Approx(south.vector().norm()) == 1.);
        REQUIRE(Approx(south.x()) == -1.);
        REQUIRE(Approx(south.phi()) == PI);
        REQUIRE(Approx(south.sinTheta()) == 1.);
        REQUIRE(Approx(south.cosTheta()).margin(1e-12) == 0.);

        Direction3d northEast{PI / 4., PI / 2.};
        REQUIRE(Approx(northEast.vector().norm()) == 1.);
        REQUIRE(Approx(northEast.x()) == 0.5 * std::sqrt(2.));
        REQUIRE(Approx(northEast.phi()) == PI / 4.);
        REQUIRE(Approx(northEast.sinTheta()) == 1.);
        REQUIRE(Approx(northEast.cosTheta()).margin(1e-12) == 0.);
    }

    SECTION("Rotation")
    {
        SECTION("Upward direction")
        {
            Direction3d const result1 = Direction3d(0., 0.).rotate({PI / 2, 0.});
            REQUIRE(Approx(result1.phi()) == PI / 2);
            REQUIRE(Approx(result1.cosTheta()) == 1.);

            Direction3d const result2 = Direction3d(0., 0.).rotate({PI / 2, PI / 2});
            REQUIRE(Approx(result2.phi()) == PI / 2);
            REQUIRE(Approx(result2.sinTheta()) == 1.);

            Direction3d const result3 = Direction3d(0., 0.).rotate({2 * PI / 3, 3 * PI / 4});
            REQUIRE(Approx(result3.phi()) == 2 * PI / 3);
            REQUIRE(Approx(result3.sinTheta()) == std::sqrt(0.5));
            REQUIRE(Approx(result3.cosTheta()) == -std::sqrt(0.5));

            Direction3d const result4 = Direction3d(0., 0.).rotate({7 * PI / 6, PI});
            REQUIRE(Approx(result4.phi()) == 7 * PI / 6);
            REQUIRE(Approx(result4.cosTheta()) == -1.);
        }

        SECTION("Tangent direction")
        {
            Direction3d const result1 = Direction3d(0, PI / 2).rotate({PI / 2, 0.});
            REQUIRE(Approx(result1.phi()) == PI / 2);
            REQUIRE(Approx(result1.sinTheta()) == 1.);

            Direction3d const result2 = Direction3d(0, PI / 2).rotate({0, PI / 2});
            REQUIRE(Approx(result2.phi()) == 0);
            REQUIRE(Approx(result2.cosTheta()) == 1.);

            Direction3d const result3 = Direction3d(0, PI / 2).rotate({PI / 2, PI / 2});
            REQUIRE(Approx(result3.phi()) == 0);
            REQUIRE(Approx(result3.cosTheta()) == 1.);

            Direction3d const result4 = Direction3d(PI / 2, PI / 2).rotate({0, PI / 2});
            REQUIRE(Approx(result4.phi()) == PI / 2);
            REQUIRE(Approx(result4.sinTheta()) == 1.);

            Direction3d const result5 = Direction3d(0, PI / 2).rotate({2 * PI / 3, 3 * PI / 4});
            REQUIRE(Approx(result5.phi()) == 2 * PI / 3);
            REQUIRE(Approx(result5.sinTheta()) == std::sqrt(0.5));
            REQUIRE(Approx(result5.cosTheta()) == std::sqrt(0.5));

            Direction3d const result6 = Direction3d(-PI / 3, PI / 2).rotate({7 * PI / 6, PI});
            REQUIRE(Approx(result6.phi()) == 5 * PI / 6);
            REQUIRE(Approx(result6.sinTheta()) == 1.0);
        }
    }
}
