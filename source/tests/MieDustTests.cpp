#include "../third-party/catch2/catch.hpp"
#include "../MieDust.hpp"

#include <sstream>

#include "../MathUtils.hpp"
#include "../Predefines.hpp"

TEST_CASE("Mie Dust", "[dust]")
{
    SECTION("Rayleigh scattering")
    {
        std::stringstream data;

        for (std::size_t i = 0; i <= 180; i += 2)
        {
            auto const theta = static_cast<double>(i);
            double const cosTheta = std::cos(radians(theta));
            data << theta << "\t" << 1.5 * cosTheta * cosTheta << "\t" << 1.5 << "\t";
            data << 1.5 * cosTheta << "\t" << 0.0 << "\n";
        }

        IDustCPtr dust = std::make_shared<MieDust>(0.5, std::move(data));

        SECTION("MatrixElements: cosine upper bound")
        {
            double p1, p2, p3, p4;
            dust->scatteringMatrixElements(p1, p2, p3, p4, 1.001);
            REQUIRE(Approx(p1).epsilon(0.001) == 1.5);
            REQUIRE(Approx(p2).margin(1e-14) == 0.0);
            REQUIRE(Approx(p3).epsilon(0.001) == 1.5);
            REQUIRE(Approx(p4).margin(1e-14) == 0.0);

            REQUIRE(Approx(dust->fraction(1.001)).epsilon(0.001) == 1.5);
        }

        SECTION("MatrixElements")
        {
            for (std::size_t i = 0; i <= 180; ++i)
            {
                auto const theta = static_cast<double>(i);
                double const cosTheta = std::cos(radians(theta));

                double p1, p2, p3, p4;
                dust->scatteringMatrixElements(p1, p2, p3, p4, cosTheta);

                REQUIRE(Approx(p1).epsilon(0.001) == 0.75 * (cosTheta*cosTheta + 1));
                REQUIRE(Approx(p2).margin(1e-14).epsilon(0.005) == 0.75 * (cosTheta*cosTheta - 1));
                REQUIRE(Approx(p3).margin(1e-14).epsilon(0.005) == 1.5 * cosTheta);
                REQUIRE(Approx(p4).margin(1e-14) == 0.0);

                REQUIRE(Approx(dust->fraction(cosTheta)).epsilon(0.001) == 0.75 * (cosTheta * cosTheta + 1));
            }
        }

        SECTION("MatrixElements: cosine lower bound")
        {
            double p1, p2, p3, p4;
            dust->scatteringMatrixElements(p1, p2, p3, p4, -1.001);
            REQUIRE(Approx(p1).epsilon(0.001) == 1.5);
            REQUIRE(Approx(p2).margin(1e-14) == 0.0);
            REQUIRE(Approx(p3).epsilon(0.001) == -1.5);
            REQUIRE(Approx(p4).margin(1e-14) == 0.0);

            REQUIRE(Approx(dust->fraction(-1.001)).epsilon(0.001) == 1.5);
        }

        SECTION("cosRandomTheta")
        {
            REQUIRE(Approx(dust->cosRandomTheta(-0.0001)).epsilon(0.001) == 1.0);

            for (int i = 0; i <= 200; ++i)
            {
                auto const fraction = i / 200.;
                // Cardano's formula to get cosTheta from fraction
                double const q = 4. * fraction - 2.;
                double const Q = 1 + q * q;
                double const cosTheta = std::cbrt(-q + std::sqrt(Q)) + std::cbrt(-q - std::sqrt(Q));

                REQUIRE(Approx(dust->cosRandomTheta(fraction)).margin(1e-14).epsilon(0.001) == cosTheta);
            }

            REQUIRE(Approx(dust->cosRandomTheta(1.0001)).epsilon(0.001) == -1.0);
        }
    }
}
