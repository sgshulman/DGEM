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

        {
            double p1, p2, p3, p4;
            dust->scatteringMatrixElements(p1, p2, p3, p4, 1.);
            REQUIRE(Approx(p1).epsilon(0.001) == 1.5);
            REQUIRE(Approx(p2).margin(1e-12) == 0.0);
            REQUIRE(Approx(p3).epsilon(0.001) == 1.5);
            REQUIRE(Approx(p4).margin(1e-12) == 0.0);
        }
    }

}
