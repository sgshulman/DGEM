#include "../third-party/catch2/catch.hpp"
#include "../Directions.hpp"
#include "../Direction3d.hpp"

namespace
{
    double sph_legendre(std::int32_t const l, std::int32_t const m, double const theta, double const phi)
    {
        if (m == 0)
        {
            return std::sqrt((2*l + 1)/(4*PI)) * std::assoc_legendre(l, m, std::cos(theta));
        } else if (m > 0) {
            return std::sqrt(2) * std::sph_legendre(l, m, theta) * std::cos(m * phi);
        }
        return std::sqrt(2) * std::sph_legendre(l, -m, theta) * std::sin(-m * phi);
    }

    double integrate(Directions const& directions, std::int32_t const l, std::int32_t const m)
    {
        double result{};
        double const w = 4 * PI / directions.number();

        for (std::uint64_t i=0; i!=directions.number(); ++i)
        {
            auto const& direction = Direction3d(directions.direction(i));
            double const y = sph_legendre(l, m, std::acos(direction.cosTheta()), direction.phi());
            result += y * y * directions.w(i) * w;
        }
        return result;
    }

    double dir_diff(Vector3d const& direction1, Vector3d const& direction2)
    {
        return (direction1 - direction2).norm();
    }
}

TEST_CASE("Icosahedron Integration", "[sphere integration]")
{
    SECTION("Icosahedron 1")
    {
        Directions d(1, false);
        REQUIRE(Approx(integrate(d, 0, 0)) == 1.);

        REQUIRE(Approx(integrate(d, 1, 1)) == 1.);
        REQUIRE(Approx(integrate(d, 1, 0)) == 1.);
        REQUIRE(Approx(integrate(d, 1, -1)) == 1.);

        REQUIRE(Approx(integrate(d, 2, 2)) == 1.);
        REQUIRE(Approx(integrate(d, 2, 1)) == 1.);
        REQUIRE(Approx(integrate(d, 2, 0)) == 1.);
        REQUIRE(Approx(integrate(d, 2, -1)) == 1.);
        REQUIRE(Approx(integrate(d, 2, -2)) == 1.);

        REQUIRE(Approx(integrate(d, 3, 3)).epsilon(0.04) == 1.);
        REQUIRE(Approx(integrate(d, 3, 2)).epsilon(0.3) == 1.);
        REQUIRE(Approx(integrate(d, 3, 1)).epsilon(0.6) == 1.);
        REQUIRE(Approx(integrate(d, 3, -1)).epsilon(0.6) == 1.);
        REQUIRE(Approx(integrate(d, 3, -2)).epsilon(0.3) == 1.);
        REQUIRE(Approx(integrate(d, 3, -3)).epsilon(0.04) == 1.);
    }

    SECTION("Icosahedron 4")
    {
        Directions d(4, false);
        REQUIRE(Approx(integrate(d, 0, 0)) == 1.);

        REQUIRE(Approx(integrate(d, 1, 1)) == 1.);
        REQUIRE(Approx(integrate(d, 1, 0)) == 1.);
        REQUIRE(Approx(integrate(d, 1, -1)) == 1.);

        REQUIRE(Approx(integrate(d, 2, 2)) == 1.);
        REQUIRE(Approx(integrate(d, 2, 1)) == 1.);
        REQUIRE(Approx(integrate(d, 2, 0)) == 1.);
        REQUIRE(Approx(integrate(d, 2, -1)) == 1.);
        REQUIRE(Approx(integrate(d, 2, -2)) == 1.);

        REQUIRE(Approx(integrate(d, 3, 3)) == 1.);
        REQUIRE(Approx(integrate(d, 3, 2)).epsilon(3e-5) == 1.);
        REQUIRE(Approx(integrate(d, 3, 1)).epsilon(6e-5) == 1.);
        REQUIRE(Approx(integrate(d, 3, 0)).epsilon(1e-4) == 1.);
        REQUIRE(Approx(integrate(d, 3, -1)).epsilon(6e-5) == 1.);
        REQUIRE(Approx(integrate(d, 3, -2)).epsilon(3e-5) == 1.);
        REQUIRE(Approx(integrate(d, 3, -3)) == 1.);

        REQUIRE(Approx(integrate(d, 4, 4)) == 1.);
        REQUIRE(Approx(integrate(d, 4, 3)).epsilon(5e-5) == 1.);
        REQUIRE(Approx(integrate(d, 4, 2)).epsilon(6e-5) == 1.);
        REQUIRE(Approx(integrate(d, 4, 1)) == 1.);
        REQUIRE(Approx(integrate(d, 4, 0)).epsilon(5e-5) == 1.);
        REQUIRE(Approx(integrate(d, 4, 1)) == 1.);
        REQUIRE(Approx(integrate(d, 4, -2)).epsilon(6e-5) == 1.);
        REQUIRE(Approx(integrate(d, 4, -3)).epsilon(5e-5) == 1.);
        REQUIRE(Approx(integrate(d, 4, -4)) == 1.);

        REQUIRE(Approx(integrate(d, 5, 5)).epsilon(4e-3) == 1.);
        REQUIRE(Approx(integrate(d, 5, 4)).epsilon(2e-4) == 1.);
        REQUIRE(Approx(integrate(d, 5, 3)).epsilon(5e-4) == 1.);
        REQUIRE(Approx(integrate(d, 5, 2)).epsilon(1e-3) == 1.);
        REQUIRE(Approx(integrate(d, 5, 1)).epsilon(2e-3) == 1.);
        REQUIRE(Approx(integrate(d, 5, 0)).epsilon(3e-3) == 1.);
        REQUIRE(Approx(integrate(d, 5, 1)).epsilon(2e-3) == 1.);
        REQUIRE(Approx(integrate(d, 5, -2)).epsilon(1e-3) == 1.);
        REQUIRE(Approx(integrate(d, 5, -3)).epsilon(5e-4) == 1.);
        REQUIRE(Approx(integrate(d, 5, -4)).epsilon(2e-4) == 1.);
        REQUIRE(Approx(integrate(d, 5, -5)).epsilon(4e-3) == 1.);
    }

    SECTION("Icosahedron 7")
    {
        Directions d(7, false);
        REQUIRE(Approx(integrate(d, 0, 0)) == 1.);

        REQUIRE(Approx(integrate(d, 1, 1)) == 1.);
        REQUIRE(Approx(integrate(d, 1, 0)) == 1.);
        REQUIRE(Approx(integrate(d, 1, -1)) == 1.);

        REQUIRE(Approx(integrate(d, 2, 2)) == 1.);
        REQUIRE(Approx(integrate(d, 2, 1)) == 1.);
        REQUIRE(Approx(integrate(d, 2, 0)) == 1.);
        REQUIRE(Approx(integrate(d, 2, -1)) == 1.);
        REQUIRE(Approx(integrate(d, 2, -2)) == 1.);

        REQUIRE(Approx(integrate(d, 3, 3)) == 1.);
        REQUIRE(Approx(integrate(d, 3, 2)) == 1.);
        REQUIRE(Approx(integrate(d, 3, 1)) == 1.);
        REQUIRE(Approx(integrate(d, 3, 0)) == 1.);
        REQUIRE(Approx(integrate(d, 3, -1)) == 1.);
        REQUIRE(Approx(integrate(d, 3, -2)) == 1.);
        REQUIRE(Approx(integrate(d, 3, -3)) == 1.);

        REQUIRE(Approx(integrate(d, 4, 4)) == 1.);
        REQUIRE(Approx(integrate(d, 4, 3)) == 1.);
        REQUIRE(Approx(integrate(d, 4, 2)) == 1.);
        REQUIRE(Approx(integrate(d, 4, 1)) == 1.);
        REQUIRE(Approx(integrate(d, 4, 0)) == 1.);
        REQUIRE(Approx(integrate(d, 4, -1)) == 1.);
        REQUIRE(Approx(integrate(d, 4, -2)) == 1.);
        REQUIRE(Approx(integrate(d, 4, -3)) == 1.);
        REQUIRE(Approx(integrate(d, 4, -4)) == 1.);

        REQUIRE(Approx(integrate(d, 5, 5)).epsilon(5e-5) == 1.);
        REQUIRE(Approx(integrate(d, 5, 4)) == 1.);
        REQUIRE(Approx(integrate(d, 5, 3)) == 1.);
        REQUIRE(Approx(integrate(d, 5, 2)).epsilon(2e-5) == 1.);
        REQUIRE(Approx(integrate(d, 5, 1)).epsilon(3e-5) == 1.);
        REQUIRE(Approx(integrate(d, 5, 0)).epsilon(4e-5) == 1.);
        REQUIRE(Approx(integrate(d, 5, -1)).epsilon(3e-5) == 1.);
        REQUIRE(Approx(integrate(d, 5, -2)).epsilon(2e-5) == 1.);
        REQUIRE(Approx(integrate(d, 5, -3)) == 1.);
        REQUIRE(Approx(integrate(d, 5, -4)) == 1.);
        REQUIRE(Approx(integrate(d, 5, -5)).epsilon(5e-5) == 1.);

        REQUIRE(Approx(integrate(d, 6, 6)) == 1.);
        REQUIRE(Approx(integrate(d, 6, 5)).epsilon(7e-5) == 1.);
        REQUIRE(Approx(integrate(d, 6, 4)) == 1.);
        REQUIRE(Approx(integrate(d, 6, 3)) == 1.);
        REQUIRE(Approx(integrate(d, 6, 2)) == 1.);
        REQUIRE(Approx(integrate(d, 6, 1)).epsilon(4e-5) == 1.);
        REQUIRE(Approx(integrate(d, 6, 0)).epsilon(5e-5) == 1.);
        REQUIRE(Approx(integrate(d, 6, -1)).epsilon(4e-5) == 1.);
        REQUIRE(Approx(integrate(d, 6, -2)) == 1.);
        REQUIRE(Approx(integrate(d, 6, -3)) == 1.);
        REQUIRE(Approx(integrate(d, 6, -4)) == 1.);
        REQUIRE(Approx(integrate(d, 6, -5)).epsilon(7e-5) == 1.);
        REQUIRE(Approx(integrate(d, 6, -6)) == 1.);
    }
}


TEST_CASE("HEALPix Integration", "[sphere integration]")
{
    SECTION("HEALPix 2")
    {
        Directions d(2, true);
        REQUIRE(Approx(integrate(d, 0, 0)) == 1.);

        REQUIRE(Approx(integrate(d, 1, 1)).epsilon(0.02) == 1.);
        REQUIRE(Approx(integrate(d, 1, 0)).epsilon(0.03) == 1.);
        REQUIRE(Approx(integrate(d, 1, -1)).epsilon(0.02) == 1.);

        REQUIRE(Approx(integrate(d, 2, 2)).epsilon(0.001)  == 1.);
        REQUIRE(Approx(integrate(d, 2, 1)).epsilon(0.04) == 1.);
        REQUIRE(Approx(integrate(d, 2, 0)).epsilon(0.09) == 1.);
        REQUIRE(Approx(integrate(d, 2, -1)).epsilon(0.04) == 1.);
        REQUIRE(Approx(integrate(d, 2, -2)).epsilon(0.02) == 1.);

        REQUIRE(Approx(integrate(d, 3, 3)).epsilon(0.004) == 1.);
        REQUIRE(Approx(integrate(d, 3, 2)).epsilon(0.02) == 1.);
        REQUIRE(Approx(integrate(d, 3, 1)).epsilon(0.02) == 1.);
        REQUIRE(Approx(integrate(d, 3, 0)).epsilon(0.12) == 1.);
        REQUIRE(Approx(integrate(d, 3, -1)).epsilon(0.02) == 1.);
        REQUIRE(Approx(integrate(d, 3, -2)).epsilon(0.08) == 1.);
        REQUIRE(Approx(integrate(d, 3, -3)).epsilon(0.004) == 1.);
    }

    SECTION("HEALPix 3")
    {
        Directions d(3, true);
        REQUIRE(Approx(integrate(d, 0, 0)) == 1.);

        REQUIRE(Approx(integrate(d, 1, 1)).epsilon(0.006) == 1.);
        REQUIRE(Approx(integrate(d, 1, 0)).epsilon(0.02) == 1.);
        REQUIRE(Approx(integrate(d, 1, -1)).epsilon(0.006) == 1.);

        REQUIRE(Approx(integrate(d, 2, 2)).epsilon(0.003)  == 1.);
        REQUIRE(Approx(integrate(d, 2, 1)).epsilon(0.02) == 1.);
        REQUIRE(Approx(integrate(d, 2, 0)).epsilon(0.04) == 1.);
        REQUIRE(Approx(integrate(d, 2, -1)).epsilon(0.02) == 1.);
        REQUIRE(Approx(integrate(d, 2, -2)).epsilon(0.005) == 1.);

        REQUIRE(Approx(integrate(d, 3, 3)).epsilon(0.002) == 1.);
        REQUIRE(Approx(integrate(d, 3, 2)).epsilon(0.008) == 1.);
        REQUIRE(Approx(integrate(d, 3, 1)).epsilon(0.008) == 1.);
        REQUIRE(Approx(integrate(d, 3, 0)).epsilon(0.05) == 1.);
        REQUIRE(Approx(integrate(d, 3, -1)).epsilon(0.008) == 1.);
        REQUIRE(Approx(integrate(d, 3, -2)).epsilon(0.02) == 1.);
        REQUIRE(Approx(integrate(d, 3, -3)).epsilon(0.002) == 1.);
    }

    SECTION("HEALPix 20")
    {
        Directions d(20, true);
        REQUIRE(Approx(integrate(d, 0, 0)) == 1.);

        REQUIRE(Approx(integrate(d, 1, 1)).epsilon(2e-4) == 1.);
        REQUIRE(Approx(integrate(d, 1, 0)).epsilon(3e-4) == 1.);
        REQUIRE(Approx(integrate(d, 1, -1)).epsilon(2e-4) == 1.);

        REQUIRE(Approx(integrate(d, 2, 2)).epsilon(1e-4) == 1.);
        REQUIRE(Approx(integrate(d, 2, 1)).epsilon(3e-4) == 1.);
        REQUIRE(Approx(integrate(d, 2, 0)).epsilon(7e-4) == 1.);
        REQUIRE(Approx(integrate(d, 2, -1)).epsilon(3e-4) == 1.);
        REQUIRE(Approx(integrate(d, 2, -2)).epsilon(1e-4) == 1.);

        REQUIRE(Approx(integrate(d, 3, 3)).epsilon(6e-5) == 1.);
        REQUIRE(Approx(integrate(d, 3, 2)).epsilon(3e-4) == 1.);
        REQUIRE(Approx(integrate(d, 3, 1)).epsilon(2e-4) == 1.);
        REQUIRE(Approx(integrate(d, 3, 0)).epsilon(1e-3) == 1.);
        REQUIRE(Approx(integrate(d, 3, -1)).epsilon(2e-4) == 1.);
        REQUIRE(Approx(integrate(d, 3, -2)).epsilon(3e-4) == 1.);
        REQUIRE(Approx(integrate(d, 3, -3)).epsilon(6e-5) == 1.);

        REQUIRE(Approx(integrate(d, 4, 4)).epsilon(4e-5) == 1.);
        REQUIRE(Approx(integrate(d, 4, 3)).epsilon(3e-4) == 1.);
        REQUIRE(Approx(integrate(d, 4, 2)).epsilon(3e-4) == 1.);
        REQUIRE(Approx(integrate(d, 4, 1)) == 1.);
        REQUIRE(Approx(integrate(d, 4, 0)).epsilon(2e-3) == 1.);
        REQUIRE(Approx(integrate(d, 4, -1)) == 1.);
        REQUIRE(Approx(integrate(d, 4, -2)).epsilon(3e-4) == 1.);
        REQUIRE(Approx(integrate(d, 4, -3)).epsilon(3e-4) == 1.);
        REQUIRE(Approx(integrate(d, 4, -4)).epsilon(4e-5) == 1.);

        REQUIRE(Approx(integrate(d, 5, 5)).epsilon(2e-5) == 1.);
        REQUIRE(Approx(integrate(d, 5, 4)).epsilon(2e-4) == 1.);
        REQUIRE(Approx(integrate(d, 5, 3)).epsilon(4e-4) == 1.);
        REQUIRE(Approx(integrate(d, 5, 2)).epsilon(8e-5) == 1.);
        REQUIRE(Approx(integrate(d, 5, 1)).epsilon(2e-4)  == 1.);
        REQUIRE(Approx(integrate(d, 5, 0)).epsilon(2e-3)  == 1.);
        REQUIRE(Approx(integrate(d, 5, -1)).epsilon(2e-4)  == 1.);
        REQUIRE(Approx(integrate(d, 5, -2)).epsilon(8e-5) == 1.);
        REQUIRE(Approx(integrate(d, 5, -3)).epsilon(4e-4) == 1.);
        REQUIRE(Approx(integrate(d, 5, -4)).epsilon(2e-4) == 1.);
        REQUIRE(Approx(integrate(d, 5, -5)).epsilon(2e-5) == 1.);
    }

    SECTION("HEALPix 80")
    {
        Directions d(80, true);
        REQUIRE(Approx(integrate(d, 0, 0)) == 1.);

        REQUIRE(Approx(integrate(d, 1, 1)) == 1.);
        REQUIRE(Approx(integrate(d, 1, 0)).epsilon(2e-5) == 1.);
        REQUIRE(Approx(integrate(d, 1, -1)) == 1.);

        REQUIRE(Approx(integrate(d, 2, 2)) == 1.);
        REQUIRE(Approx(integrate(d, 2, 1)).epsilon(2e-5) == 1.);
        REQUIRE(Approx(integrate(d, 2, 0)).epsilon(5e-5) == 1.);
        REQUIRE(Approx(integrate(d, 2, -1)).epsilon(2e-5) == 1.);
        REQUIRE(Approx(integrate(d, 2, -2)) == 1.);

        REQUIRE(Approx(integrate(d, 3, 3)) == 1.);
        REQUIRE(Approx(integrate(d, 3, 2)).epsilon(2e-5) == 1.);
        REQUIRE(Approx(integrate(d, 3, 1)) == 1.);
        REQUIRE(Approx(integrate(d, 3, 0)).epsilon(6e-5) == 1.);
        REQUIRE(Approx(integrate(d, 3, -1)) == 1.);
        REQUIRE(Approx(integrate(d, 3, -2)).epsilon(2e-5) == 1.);
        REQUIRE(Approx(integrate(d, 3, -3)) == 1.);

        REQUIRE(Approx(integrate(d, 4, 4)) == 1.);
        REQUIRE(Approx(integrate(d, 4, 3)).epsilon(2e-5) == 1.);
        REQUIRE(Approx(integrate(d, 4, 2)).epsilon(2e-5) == 1.);
        REQUIRE(Approx(integrate(d, 4, 1)) == 1.);
        REQUIRE(Approx(integrate(d, 4, 0)).epsilon(7e-5) == 1.);
        REQUIRE(Approx(integrate(d, 4, -1)) == 1.);
        REQUIRE(Approx(integrate(d, 4, -2)).epsilon(2e-5) == 1.);
        REQUIRE(Approx(integrate(d, 4, -3)).epsilon(2e-5) == 1.);
        REQUIRE(Approx(integrate(d, 4, -4)) == 1.);

        REQUIRE(Approx(integrate(d, 5, 5)) == 1.);
        REQUIRE(Approx(integrate(d, 5, 4)) == 1.);
        REQUIRE(Approx(integrate(d, 5, 3)).epsilon(3e-5) == 1.);
        REQUIRE(Approx(integrate(d, 5, 2)) == 1.);
        REQUIRE(Approx(integrate(d, 5, 1)) == 1.);
        REQUIRE(Approx(integrate(d, 5, 0)).epsilon(1e-4)  == 1.);
        REQUIRE(Approx(integrate(d, 5, -1)) == 1.);
        REQUIRE(Approx(integrate(d, 5, -2)) == 1.);
        REQUIRE(Approx(integrate(d, 5, -3)).epsilon(3e-5) == 1.);
        REQUIRE(Approx(integrate(d, 5, -4)) == 1.);
        REQUIRE(Approx(integrate(d, 5, -5)) == 1.);

        REQUIRE(Approx(integrate(d, 6, 6)).epsilon(3e-5) == 1.);
        REQUIRE(Approx(integrate(d, 6, 5)) == 1.);
        REQUIRE(Approx(integrate(d, 6, 4)).epsilon(3e-5) == 1.);
        REQUIRE(Approx(integrate(d, 6, 3)).epsilon(2e-5) == 1.);
        REQUIRE(Approx(integrate(d, 6, 2)) == 1.);
        REQUIRE(Approx(integrate(d, 6, 1)).epsilon(2e-5) == 1.);
        REQUIRE(Approx(integrate(d, 6, 0)).epsilon(2e-4) == 1.);
        REQUIRE(Approx(integrate(d, 6, -1)).epsilon(2e-5) == 1.);
        REQUIRE(Approx(integrate(d, 6, -2)) == 1.);
        REQUIRE(Approx(integrate(d, 6, -3)).epsilon(2e-5) == 1.);
        REQUIRE(Approx(integrate(d, 6, -4)).epsilon(3e-5) == 1.);
        REQUIRE(Approx(integrate(d, 6, -5)) == 1.);
        REQUIRE(Approx(integrate(d, 6, -6)) == 1.);
    }
}


TEST_CASE("Isolatitude 3-6 Integration", "[sphere integration]")
{
    SECTION("Sphere 2")
    {
        Directions d(3, 6, 2, true);
        REQUIRE(Approx(integrate(d, 0, 0)) == 1.);

        REQUIRE(Approx(integrate(d, 0, 0)) == 1.);

        REQUIRE(Approx(integrate(d, 1, 1)).epsilon(0.02) == 1.);
        REQUIRE(Approx(integrate(d, 1, 0)).epsilon(0.03) == 1.);
        REQUIRE(Approx(integrate(d, 1, -1)).epsilon(0.02) == 1.);

        REQUIRE(Approx(integrate(d, 2, 2)).epsilon(0.008)  == 1.);
        REQUIRE(Approx(integrate(d, 2, 1)).epsilon(0.04) == 1.);
        REQUIRE(Approx(integrate(d, 2, 0)).epsilon(0.09) == 1.);
        REQUIRE(Approx(integrate(d, 2, -1)).epsilon(0.04) == 1.);
        REQUIRE(Approx(integrate(d, 2, -2)).epsilon(0.008) == 1.);

        REQUIRE(Approx(integrate(d, 3, 3)).epsilon(0.002) == 1.);
        REQUIRE(Approx(integrate(d, 3, 2)).epsilon(0.04) == 1.);
        REQUIRE(Approx(integrate(d, 3, 1)).epsilon(0.02) == 1.);
        REQUIRE(Approx(integrate(d, 3, 0)).epsilon(0.12) == 1.);
        REQUIRE(Approx(integrate(d, 3, -1)).epsilon(0.02) == 1.);
        REQUIRE(Approx(integrate(d, 3, -2)).epsilon(0.04) == 1.);
        REQUIRE(Approx(integrate(d, 3, -3)).epsilon(0.005) == 1.);
    }

    SECTION("Sphere 3")
    {
        Directions d(3, 6, 3, true);
        REQUIRE(Approx(integrate(d, 0, 0)) == 1.);

        REQUIRE(Approx(integrate(d, 1, 1)).epsilon(0.006) == 1.);
        REQUIRE(Approx(integrate(d, 1, 0)).epsilon(0.02) == 1.);
        REQUIRE(Approx(integrate(d, 1, -1)).epsilon(0.006) == 1.);

        REQUIRE(Approx(integrate(d, 2, 2)).epsilon(0.004)  == 1.);
        REQUIRE(Approx(integrate(d, 2, 1)).epsilon(0.02) == 1.);
        REQUIRE(Approx(integrate(d, 2, 0)).epsilon(0.04) == 1.);
        REQUIRE(Approx(integrate(d, 2, -1)).epsilon(0.02) == 1.);
        REQUIRE(Approx(integrate(d, 2, -2)).epsilon(0.004) == 1.);

        REQUIRE(Approx(integrate(d, 3, 3)).epsilon(0.002) == 1.);
        REQUIRE(Approx(integrate(d, 3, 2)).epsilon(0.02) == 1.);
        REQUIRE(Approx(integrate(d, 3, 1)).epsilon(0.008) == 1.);
        REQUIRE(Approx(integrate(d, 3, 0)).epsilon(0.05) == 1.);
        REQUIRE(Approx(integrate(d, 3, -1)).epsilon(0.008) == 1.);
        REQUIRE(Approx(integrate(d, 3, -2)).epsilon(0.02) == 1.);
        REQUIRE(Approx(integrate(d, 3, -3)).epsilon(0.003) == 1.);
    }
}


TEST_CASE("Isolatitude 2-4 Integration", "[sphere integration]")
{
    SECTION("Sphere 3")
    {
        Directions d(2, 4, 3, true);
        REQUIRE(Approx(integrate(d, 0, 0)) == 1.);

        REQUIRE(Approx(integrate(d, 0, 0)) == 1.);

        REQUIRE(Approx(integrate(d, 1, 1)).epsilon(0.02) == 1.);
        REQUIRE(Approx(integrate(d, 1, 0)).epsilon(0.03) == 1.);
        REQUIRE(Approx(integrate(d, 1, -1)).epsilon(0.02) == 1.);

        REQUIRE(Approx(integrate(d, 2, 2)).epsilon(0.01)  == 1.);
        REQUIRE(Approx(integrate(d, 2, 1)).epsilon(0.02) == 1.);
        REQUIRE(Approx(integrate(d, 2, 0)).epsilon(0.06) == 1.);
        REQUIRE(Approx(integrate(d, 2, -1)).epsilon(0.02) == 1.);
        REQUIRE(Approx(integrate(d, 2, -2)).epsilon(0.01) == 1.);

        REQUIRE(Approx(integrate(d, 3, 3)).epsilon(0.008) == 1.);
        REQUIRE(Approx(integrate(d, 3, 2)).epsilon(0.03) == 1.);
        REQUIRE(Approx(integrate(d, 3, 1)).epsilon(0.004) == 1.);
        REQUIRE(Approx(integrate(d, 3, 0)).epsilon(0.06) == 1.);
        REQUIRE(Approx(integrate(d, 3, -1)).epsilon(0.004) == 1.);
        REQUIRE(Approx(integrate(d, 3, -2)).epsilon(0.03) == 1.);
        REQUIRE(Approx(integrate(d, 3, -3)).epsilon(0.008) == 1.);
    }

    SECTION("Sphere 20")
    {
        Directions d(2, 4, 20, true);
        REQUIRE(Approx(integrate(d, 0, 0)) == 1.);

        REQUIRE(Approx(integrate(d, 1, 1)).epsilon(3e-4) == 1.);
        REQUIRE(Approx(integrate(d, 1, 0)).epsilon(5e-4) == 1.);
        REQUIRE(Approx(integrate(d, 1, -1)).epsilon(3e-4) == 1.);

        REQUIRE(Approx(integrate(d, 2, 2)).epsilon(3e-4)  == 1.);
        REQUIRE(Approx(integrate(d, 2, 1)).epsilon(3e-4) == 1.);
        REQUIRE(Approx(integrate(d, 2, 0)).epsilon(2e-3) == 1.);
        REQUIRE(Approx(integrate(d, 2, -1)).epsilon(3e-4) == 1.);
        REQUIRE(Approx(integrate(d, 2, -2)).epsilon(3e-4) == 1.);

        REQUIRE(Approx(integrate(d, 3, 3)).epsilon(2e-4) == 1.);
        REQUIRE(Approx(integrate(d, 3, 2)).epsilon(4e-4) == 1.);
        REQUIRE(Approx(integrate(d, 3, 1)).epsilon(2e-5) == 1.);
        REQUIRE(Approx(integrate(d, 3, 0)).epsilon(2e-3) == 1.);
        REQUIRE(Approx(integrate(d, 3, -1)).epsilon(2e-5) == 1.);
        REQUIRE(Approx(integrate(d, 3, -2)).epsilon(4e-4) == 1.);
        REQUIRE(Approx(integrate(d, 3, -3)).epsilon(2e-4) == 1.);

        REQUIRE(Approx(integrate(d, 4, 4)).epsilon(2e-4) == 1.);
        REQUIRE(Approx(integrate(d, 4, 3)).epsilon(5e-4) == 1.);
        REQUIRE(Approx(integrate(d, 4, 2)).epsilon(1e-4) == 1.);
        REQUIRE(Approx(integrate(d, 4, 1)).epsilon(2e-4) == 1.);
        REQUIRE(Approx(integrate(d, 4, 0)).epsilon(2e-3) == 1.);
        REQUIRE(Approx(integrate(d, 4, -1)).epsilon(2e-4) == 1.);
        REQUIRE(Approx(integrate(d, 4, -2)).epsilon(1e-4) == 1.);
        REQUIRE(Approx(integrate(d, 4, -3)).epsilon(5e-4) == 1.);
        REQUIRE(Approx(integrate(d, 4, -4)).epsilon(2e-4) == 1.);

        REQUIRE(Approx(integrate(d, 5, 5)).epsilon(2e-4) == 1.);
        REQUIRE(Approx(integrate(d, 5, 4)).epsilon(5e-4) == 1.);
        REQUIRE(Approx(integrate(d, 5, 3)).epsilon(3e-4) == 1.);
        REQUIRE(Approx(integrate(d, 5, 2)).epsilon(7e-5) == 1.);
        REQUIRE(Approx(integrate(d, 5, 1)).epsilon(3e-4) == 1.);
        REQUIRE(Approx(integrate(d, 5, 0)).epsilon(3e-3) == 1.);
        REQUIRE(Approx(integrate(d, 5, -1)).epsilon(3e-4) == 1.);
        REQUIRE(Approx(integrate(d, 5, -2)).epsilon(7e-5) == 1.);
        REQUIRE(Approx(integrate(d, 5, -3)).epsilon(3e-4) == 1.);
        REQUIRE(Approx(integrate(d, 5, -4)).epsilon(5e-4) == 1.);
        REQUIRE(Approx(integrate(d, 5, -5)).epsilon(2e-4) == 1.);
    }
}


TEST_CASE("Isolatitude 2-4 Directions", "[directions]")
{
    SECTION("Icosahedron 1")
    {
        Directions d(2, 4, 1, true);

        REQUIRE(dir_diff(d.direction(0), {  PI/4, PI/3}) < std::numeric_limits<float>::epsilon());
        REQUIRE(dir_diff(d.direction(1), {3*PI/4, PI/3}) < std::numeric_limits<float>::epsilon());
        REQUIRE(dir_diff(d.direction(2), {5*PI/4, PI/3}) < std::numeric_limits<float>::epsilon());
        REQUIRE(dir_diff(d.direction(3), {7*PI/4, PI/3}) < std::numeric_limits<float>::epsilon());

        REQUIRE(dir_diff(d.direction(4), {   0.0, 2*PI/3}) < std::numeric_limits<float>::epsilon());
        REQUIRE(dir_diff(d.direction(5), {  PI/2, 2*PI/3}) < std::numeric_limits<float>::epsilon());
        REQUIRE(dir_diff(d.direction(6), {    PI, 2*PI/3}) < std::numeric_limits<float>::epsilon());
        REQUIRE(dir_diff(d.direction(7), {3*PI/2, 2*PI/3}) < std::numeric_limits<float>::epsilon());
    }

    SECTION("Icosahedron 2")
    {
        Directions d(2, 4, 2, true);

        REQUIRE(dir_diff(d.direction(0), {  PI/4, acos(7./8)}) < std::numeric_limits<float>::epsilon());
        REQUIRE(dir_diff(d.direction(1), {3*PI/4, acos(7./8)}) < std::numeric_limits<float>::epsilon());
        REQUIRE(dir_diff(d.direction(2), {5*PI/4, acos(7./8)}) < std::numeric_limits<float>::epsilon());
        REQUIRE(dir_diff(d.direction(3), {7*PI/4, acos(7./8)}) < std::numeric_limits<float>::epsilon());

        REQUIRE(dir_diff(d.direction( 4), {   PI/8, PI/3}) < std::numeric_limits<float>::epsilon());
        REQUIRE(dir_diff(d.direction( 5), { 3*PI/8, PI/3}) < std::numeric_limits<float>::epsilon());
        REQUIRE(dir_diff(d.direction( 6), { 5*PI/8, PI/3}) < std::numeric_limits<float>::epsilon());
        REQUIRE(dir_diff(d.direction( 7), { 7*PI/8, PI/3}) < std::numeric_limits<float>::epsilon());
        REQUIRE(dir_diff(d.direction( 8), { 9*PI/8, PI/3}) < std::numeric_limits<float>::epsilon());
        REQUIRE(dir_diff(d.direction( 9), {11*PI/8, PI/3}) < std::numeric_limits<float>::epsilon());
        REQUIRE(dir_diff(d.direction(10), {13*PI/8, PI/3}) < std::numeric_limits<float>::epsilon());
        REQUIRE(dir_diff(d.direction(11), {15*PI/8, PI/3}) < std::numeric_limits<float>::epsilon());

        REQUIRE(dir_diff(d.direction(12), {   0.0, PI/2}) < std::numeric_limits<float>::epsilon());
        REQUIRE(dir_diff(d.direction(13), {  PI/4, PI/2}) < std::numeric_limits<float>::epsilon());
        REQUIRE(dir_diff(d.direction(14), {  PI/2, PI/2}) < std::numeric_limits<float>::epsilon());
        REQUIRE(dir_diff(d.direction(15), {3*PI/4, PI/2}) < std::numeric_limits<float>::epsilon());
        REQUIRE(dir_diff(d.direction(16), {    PI, PI/2}) < std::numeric_limits<float>::epsilon());
        REQUIRE(dir_diff(d.direction(17), {5*PI/4, PI/2}) < std::numeric_limits<float>::epsilon());
        REQUIRE(dir_diff(d.direction(18), {3*PI/2, PI/2}) < std::numeric_limits<float>::epsilon());
        REQUIRE(dir_diff(d.direction(19), {7*PI/4, PI/2}) < std::numeric_limits<float>::epsilon());

        REQUIRE(dir_diff(d.direction(20), {   PI/8, 2*PI/3}) < std::numeric_limits<float>::epsilon());
        REQUIRE(dir_diff(d.direction(21), { 3*PI/8, 2*PI/3}) < std::numeric_limits<float>::epsilon());
        REQUIRE(dir_diff(d.direction(22), { 5*PI/8, 2*PI/3}) < std::numeric_limits<float>::epsilon());
        REQUIRE(dir_diff(d.direction(23), { 7*PI/8, 2*PI/3}) < std::numeric_limits<float>::epsilon());
        REQUIRE(dir_diff(d.direction(24), { 9*PI/8, 2*PI/3}) < std::numeric_limits<float>::epsilon());
        REQUIRE(dir_diff(d.direction(25), {11*PI/8, 2*PI/3}) < std::numeric_limits<float>::epsilon());
        REQUIRE(dir_diff(d.direction(26), {13*PI/8, 2*PI/3}) < std::numeric_limits<float>::epsilon());
        REQUIRE(dir_diff(d.direction(27), {15*PI/8, 2*PI/3}) < std::numeric_limits<float>::epsilon());

        INFO(d.direction(28).x());
        INFO(d.direction(28).y());
        INFO(d.direction(28).z());

        INFO(Vector3d(   0.0, acos(-7./8)).x());
        INFO(Vector3d(   0.0, acos(-7./8)).y());
        INFO(Vector3d(   0.0, acos(-7./8)).z());

        REQUIRE(dir_diff(d.direction(28), {   0.0, acos(-7./8)}) < std::numeric_limits<float>::epsilon());
        REQUIRE(dir_diff(d.direction(29), {  PI/2, acos(-7./8)}) < std::numeric_limits<float>::epsilon());
        REQUIRE(dir_diff(d.direction(30), {    PI, acos(-7./8)}) < std::numeric_limits<float>::epsilon());
        REQUIRE(dir_diff(d.direction(31), {3*PI/2, acos(-7./8)}) < std::numeric_limits<float>::epsilon());
    }
}


TEST_CASE("HEALPix 2 ring and nested schemes", "[HEALPix schemes]")
{
    Directions ring(3, 4, 2, true);
    Directions nested(3, 4, 2, false);

    REQUIRE(dir_diff(ring.direction( 0), nested.direction( 3)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 1), nested.direction( 7)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 2), nested.direction(11)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 3), nested.direction(15)) < std::numeric_limits<float>::epsilon());

    REQUIRE(dir_diff(ring.direction( 4), nested.direction( 2)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 5), nested.direction( 1)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 6), nested.direction( 6)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 7), nested.direction( 5)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 8), nested.direction(10)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 9), nested.direction( 9)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(10), nested.direction(14)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(11), nested.direction(13)) < std::numeric_limits<float>::epsilon());

    REQUIRE(dir_diff(ring.direction(12), nested.direction(19)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(13), nested.direction( 0)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(14), nested.direction(23)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(15), nested.direction( 4)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(16), nested.direction(27)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(17), nested.direction( 8)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(18), nested.direction(31)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(19), nested.direction(12)) < std::numeric_limits<float>::epsilon());

    REQUIRE(dir_diff(ring.direction(20), nested.direction(17)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(21), nested.direction(22)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(22), nested.direction(21)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(23), nested.direction(26)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(24), nested.direction(25)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(25), nested.direction(30)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(26), nested.direction(29)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(27), nested.direction(18)) < std::numeric_limits<float>::epsilon());

    REQUIRE(dir_diff(ring.direction(28), nested.direction(16)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(29), nested.direction(35)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(30), nested.direction(20)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(31), nested.direction(39)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(32), nested.direction(24)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(33), nested.direction(43)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(34), nested.direction(28)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(35), nested.direction(47)) < std::numeric_limits<float>::epsilon());

    REQUIRE(dir_diff(ring.direction(36), nested.direction(34)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(37), nested.direction(33)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(38), nested.direction(38)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(39), nested.direction(37)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(40), nested.direction(42)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(41), nested.direction(41)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(42), nested.direction(46)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(43), nested.direction(45)) < std::numeric_limits<float>::epsilon());

    REQUIRE(dir_diff(ring.direction(44), nested.direction(32)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(45), nested.direction(36)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(46), nested.direction(40)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(47), nested.direction(44)) < std::numeric_limits<float>::epsilon());
}


TEST_CASE("HEALPix 4 ring and nested schemes", "[HEALPix schemes]")
{
    Directions ring(3, 4, 4, true);
    Directions nested(3, 4, 4, false);

    REQUIRE(dir_diff(ring.direction(  0), nested.direction( 15)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(  1), nested.direction( 31)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(  2), nested.direction( 47)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(  3), nested.direction( 63)) < std::numeric_limits<float>::epsilon());

    REQUIRE(dir_diff(ring.direction(  4), nested.direction( 14)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(  5), nested.direction( 13)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(  6), nested.direction( 30)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(  7), nested.direction( 29)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(  8), nested.direction( 46)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(  9), nested.direction( 45)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 10), nested.direction( 62)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 11), nested.direction( 61)) < std::numeric_limits<float>::epsilon());

    REQUIRE(dir_diff(ring.direction( 12), nested.direction( 11)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 13), nested.direction( 12)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 14), nested.direction(  7)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 15), nested.direction( 27)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 16), nested.direction( 28)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 17), nested.direction( 23)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 18), nested.direction( 43)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 19), nested.direction( 44)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 20), nested.direction( 39)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 21), nested.direction( 59)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 22), nested.direction( 60)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 23), nested.direction( 55)) < std::numeric_limits<float>::epsilon());

    REQUIRE(dir_diff(ring.direction( 24), nested.direction( 10)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 25), nested.direction(  9)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 26), nested.direction(  6)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 27), nested.direction(  5)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 28), nested.direction( 26)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 29), nested.direction( 25)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 30), nested.direction( 22)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 31), nested.direction( 21)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 32), nested.direction( 42)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 33), nested.direction( 41)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 34), nested.direction( 38)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 35), nested.direction( 37)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 36), nested.direction( 58)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 37), nested.direction( 57)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 38), nested.direction( 54)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 39), nested.direction( 53)) < std::numeric_limits<float>::epsilon());

    REQUIRE(dir_diff(ring.direction( 40), nested.direction( 79)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 41), nested.direction(  8)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 42), nested.direction(  3)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 43), nested.direction(  4)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 44), nested.direction( 95)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 45), nested.direction( 24)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 46), nested.direction( 19)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 47), nested.direction( 20)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 48), nested.direction(111)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 49), nested.direction( 40)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 50), nested.direction( 35)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 51), nested.direction( 36)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 52), nested.direction(127)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 53), nested.direction( 56)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 54), nested.direction( 51)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 55), nested.direction( 52)) < std::numeric_limits<float>::epsilon());

    REQUIRE(dir_diff(ring.direction( 56), nested.direction( 77)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 57), nested.direction(  2)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 58), nested.direction(  1)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 59), nested.direction( 94)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 60), nested.direction( 93)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 61), nested.direction( 18)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 62), nested.direction( 17)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 63), nested.direction(110)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 64), nested.direction(109)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 65), nested.direction( 34)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 66), nested.direction( 33)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 67), nested.direction(126)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 68), nested.direction(125)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 69), nested.direction( 50)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 70), nested.direction( 49)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 71), nested.direction( 78)) < std::numeric_limits<float>::epsilon());

    REQUIRE(dir_diff(ring.direction( 72), nested.direction( 76)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 73), nested.direction( 71)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 74), nested.direction(  0)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 75), nested.direction( 91)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 76), nested.direction( 92)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 77), nested.direction( 87)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 78), nested.direction( 16)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 79), nested.direction(107)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 80), nested.direction(108)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 81), nested.direction(103)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 82), nested.direction( 32)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 83), nested.direction(123)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 84), nested.direction(124)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 85), nested.direction(119)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 86), nested.direction( 48)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 87), nested.direction( 75)) < std::numeric_limits<float>::epsilon());

    REQUIRE(dir_diff(ring.direction( 88), nested.direction( 70)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 89), nested.direction( 69)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 90), nested.direction( 90)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 91), nested.direction( 89)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 92), nested.direction( 86)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 93), nested.direction( 85)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 94), nested.direction(106)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 95), nested.direction(105)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 96), nested.direction(102)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 97), nested.direction(101)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 98), nested.direction(122)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction( 99), nested.direction(121)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(100), nested.direction(118)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(101), nested.direction(117)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(102), nested.direction( 74)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(103), nested.direction( 73)) < std::numeric_limits<float>::epsilon());

    REQUIRE(dir_diff(ring.direction(104), nested.direction( 67)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(105), nested.direction( 68)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(106), nested.direction(143)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(107), nested.direction( 88)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(108), nested.direction( 83)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(109), nested.direction( 84)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(110), nested.direction(159)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(111), nested.direction(104)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(112), nested.direction( 99)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(113), nested.direction(100)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(114), nested.direction(175)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(115), nested.direction(120)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(116), nested.direction(115)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(117), nested.direction(116)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(118), nested.direction(191)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(119), nested.direction( 72)) < std::numeric_limits<float>::epsilon());

    REQUIRE(dir_diff(ring.direction(120), nested.direction( 65)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(121), nested.direction(142)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(122), nested.direction(141)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(123), nested.direction( 82)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(124), nested.direction( 81)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(125), nested.direction(158)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(126), nested.direction(157)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(127), nested.direction( 98)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(128), nested.direction( 97)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(129), nested.direction(174)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(130), nested.direction(173)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(131), nested.direction(114)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(132), nested.direction(113)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(133), nested.direction(190)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(134), nested.direction(189)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(135), nested.direction( 66)) < std::numeric_limits<float>::epsilon());

    REQUIRE(dir_diff(ring.direction(136), nested.direction( 64)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(137), nested.direction(139)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(138), nested.direction(140)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(139), nested.direction(135)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(140), nested.direction( 80)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(141), nested.direction(155)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(142), nested.direction(156)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(143), nested.direction(151)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(144), nested.direction( 96)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(145), nested.direction(171)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(146), nested.direction(172)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(147), nested.direction(167)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(148), nested.direction(112)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(149), nested.direction(187)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(150), nested.direction(188)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(151), nested.direction(183)) < std::numeric_limits<float>::epsilon());

    REQUIRE(dir_diff(ring.direction(152), nested.direction(138)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(153), nested.direction(137)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(154), nested.direction(134)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(155), nested.direction(133)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(156), nested.direction(154)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(157), nested.direction(153)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(158), nested.direction(150)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(159), nested.direction(149)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(160), nested.direction(170)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(161), nested.direction(169)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(162), nested.direction(166)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(163), nested.direction(165)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(164), nested.direction(186)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(165), nested.direction(185)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(166), nested.direction(182)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(167), nested.direction(181)) < std::numeric_limits<float>::epsilon());

    REQUIRE(dir_diff(ring.direction(168), nested.direction(136)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(169), nested.direction(131)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(170), nested.direction(132)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(171), nested.direction(152)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(172), nested.direction(147)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(173), nested.direction(148)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(174), nested.direction(168)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(175), nested.direction(163)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(176), nested.direction(164)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(177), nested.direction(184)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(178), nested.direction(179)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(179), nested.direction(180)) < std::numeric_limits<float>::epsilon());

    REQUIRE(dir_diff(ring.direction(180), nested.direction(130)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(181), nested.direction(129)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(182), nested.direction(146)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(183), nested.direction(145)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(184), nested.direction(162)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(185), nested.direction(161)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(186), nested.direction(178)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(187), nested.direction(177)) < std::numeric_limits<float>::epsilon());

    REQUIRE(dir_diff(ring.direction(188), nested.direction(128)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(189), nested.direction(144)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(190), nested.direction(160)) < std::numeric_limits<float>::epsilon());
    REQUIRE(dir_diff(ring.direction(191), nested.direction(176)) < std::numeric_limits<float>::epsilon());
}

