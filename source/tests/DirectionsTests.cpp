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

TEST_CASE("Icosahedron Directions", "[sphere integration]")
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

TEST_CASE("HEALPix Directions", "[sphere integration]")
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


TEST_CASE("Isolatitude 3-6 Directions", "[sphere integration]")
{
    SECTION("Sphere 2")
    {
        Directions d(6, 2, true);
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
        Directions d(6, 3, true);
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


TEST_CASE("HEALPix 2 grid and nested schemes", "[HEALPix schemes]")
{
    Directions ring(4, 2, true);
    Directions nested(4, 2, false);

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
