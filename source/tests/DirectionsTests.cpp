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
}

TEST_CASE("Directions", "[sphere integration]")
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
        REQUIRE(Approx(integrate(d, 5, 1)).epsilon(3e-5)  == 1.);
        REQUIRE(Approx(integrate(d, 5, 0)).epsilon(4e-5)  == 1.);
        REQUIRE(Approx(integrate(d, 5, -1)).epsilon(3e-5)  == 1.);
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