#include "../third-party/catch2/catch.hpp"
#include "../MathUtils.hpp"

TEST_CASE("Math Utils", "[math_utils]")
{
    SECTION("Radians")
    {
        REQUIRE(Approx(radians(180.)) == PI);
        REQUIRE(Approx(radians(360.)) == 2 * PI);
        REQUIRE(Approx(radians(-45.)) == -PI / 4.);
    }

    SECTION("Degrees")
    {
        REQUIRE(Approx(degrees(0.)).margin(1e-12) == 0.);
        REQUIRE(Approx(degrees(PI)) == 180.);
        REQUIRE(Approx(degrees(-PI / 4.)) == -45.);
        REQUIRE(Approx(degrees(radians(67.))) == 67.);
    }

    SECTION("Norm Angle")
    {
        REQUIRE(Approx(normAngle(0.)).margin(1e-12) == 0.);
        REQUIRE(Approx(normAngle(2. * PI)).margin(1e-12) == 0.);
        REQUIRE(Approx(normAngle(-2. * PI)).margin(1e-12) == 0.);
        REQUIRE(Approx(normAngle(4. * PI)).margin(1e-12) == 0.);

        REQUIRE(Approx(normAngle(PI)) == PI);
        REQUIRE(Approx(normAngle(3. * PI)) == PI);
        REQUIRE(Approx(normAngle(-3. * PI)) == PI);
        REQUIRE(Approx(normAngle(11. * PI)) == PI);
    }
}