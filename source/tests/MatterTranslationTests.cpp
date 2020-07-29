#include "../third-party/catch2/catch.hpp"
#include "../MatterTranslation.hpp"
#include "../MathUtils.hpp"
#include "TestUtils.hpp"

TEST_CASE("Matter Translation", "[matter]")
{
    SECTION("Stationary")
    {
        MatterTranslation translation{0., 0., 0., Vector3d{0, 0, 0}};

        Vector3d position{0, 0, 0};
        REQUIRE(equal(position, translation(position)));

        Vector3d positionX{10., 0, 0};
        REQUIRE(equal(positionX, translation(positionX)));

        Vector3d positionY{0, 10., 0};
        REQUIRE(equal(positionY, translation(positionY)));

        Vector3d positionZ{0, 0, 10.};
        REQUIRE(equal(positionZ, translation(positionZ)));
    }

    SECTION("Translation")
    {
        MatterTranslation translation{0., 0., 0., Vector3d{1, 2, -4}};

        REQUIRE(equal(Vector3d{-1., -2., 4}, translation(Vector3d{0, 0, 0})));
        REQUIRE(equal(Vector3d{0, 0, 0}, translation(Vector3d{1., 2., -4})));
        REQUIRE(equal(Vector3d{9., -2., 4}, translation(Vector3d{10., 0, 0})));
    }

    SECTION("Nutation")
    {
        MatterTranslation translation{0., PI / 2., 0., Vector3d{0, 0, 0}};

        REQUIRE(equal(Vector3d{0., 0., 0.}, translation(Vector3d{0., 0., 0.})));
        REQUIRE(equal(Vector3d{1., 0., 0.}, translation(Vector3d{1., 0., 0.})));
        REQUIRE(equal(Vector3d{0, 1., 0.}, translation(Vector3d{0., 0., 1.})));
    }

    SECTION("Precession")
    {
        MatterTranslation translation{PI / 4., 0., 0., Vector3d{0, 0, 0}};

        REQUIRE(equal(Vector3d{0., 0., 0.}, translation(Vector3d{0., 0., 0.})));
        REQUIRE(equal(Vector3d{0., 0., 1.}, translation(Vector3d{0., 0., 1.})));
        REQUIRE(equal(Vector3d{0.5 * std::sqrt(2.), 0.5 * std::sqrt(2.), 0.}, translation(Vector3d{0., 1., 0.})));
    }

    SECTION("Intrinsic Rotation")
    {
        MatterTranslation translation{0., 0., -PI / 3., Vector3d{0, 0, 0}};

        REQUIRE(equal(Vector3d{0., 0., 0.}, translation(Vector3d{0., 0., 0.})));
        REQUIRE(equal(Vector3d{0., 0., 1.}, translation(Vector3d{0., 0., 1.})));
        REQUIRE(equal(Vector3d{-0.5 * std::sqrt(3.), 0.5, 0.}, translation(Vector3d{0., 1., 0.})));
    }

    SECTION("Full Translation")
    {
        MatterTranslation translation{-PI / 2.,PI / 3., PI / 6., Vector3d{10, 5, -5}};

        REQUIRE(equal(Vector3d{0., 0., 0.}, translation(Vector3d{10, 5, -5})));

        // (0, 0, 1) -> (0, 0, 1) -> (0, -sqrt(3)/2, 1/2) -> (sqrt(3)/4, -3/4, 1/2) -> (10. + sqrt(3)/4, 4.25, -4.5)
        REQUIRE(equal(Vector3d{0., 0., 1.}, translation(Vector3d{10.+std::sqrt(3.)/4, 4.25, -4.5})));

        // (0, 1, 0) -> (1, 0, 0) -> (1, 0, 0) -> (sqrt(3)/2, 0.5, 0)  -> (10. + sqrt(3)/2, 5.5, -5)
        REQUIRE(equal(Vector3d{0., 1., 0.}, translation(Vector3d{10.+std::sqrt(3.)/2, 5.5, -5.})));

        // (1, 0, 0) -> (0, -1, 0) -> (0, -1/2, -sqrt(3)/2) -> (1/4, -sqrt(3)/4, -sqrt(3)/2) -> (10. + sqrt(3)/4, 4.25, -5 -sqrt(3)/2)
        REQUIRE(equal(Vector3d{1., 0., 0.}, translation(Vector3d{10.+0.25, 5-std::sqrt(3.)/4, -5 - std::sqrt(3.) / 2})));
    }
}
