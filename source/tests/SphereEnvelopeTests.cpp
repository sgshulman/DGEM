#include "../third-party/catch2/catch.hpp"
#include "../SphereEnvelope.hpp"

TEST_CASE("Sphere Envelope", "[matter]")
{
    SECTION("Union Sphere")
    {
        SphereEnvelope envelope{1., 3., 1., 1., 0.};

        REQUIRE(Approx(0.).margin(1e-12) == envelope.density(0.9, 0., 0.));

        REQUIRE(Approx(1.) == envelope.density(1., 0., 0.));
        REQUIRE(Approx(1.) == envelope.density(0., 1., 0.));
        REQUIRE(Approx(1.) == envelope.density(0., 0., 1.));

        REQUIRE(Approx(1.) == envelope.density(-2., 0., 0.));
        REQUIRE(Approx(1.) == envelope.density(0., -2., 0.));
        REQUIRE(Approx(1.) == envelope.density(0., 0., -2.));

        REQUIRE(Approx(1.) == envelope.density(-3., 0., 0.));
        REQUIRE(Approx(1.) == envelope.density(0., -3., 0.));
        REQUIRE(Approx(1.) == envelope.density(0., 0., -3.));

        REQUIRE(Approx(0.).margin(1e-12) == envelope.density(3.1, 0., 0.));
        REQUIRE(Approx(0.).margin(1e-12) == envelope.density(0.0, -3.1, 0.));
        REQUIRE(Approx(0.).margin(1e-12) == envelope.density(-1., 0., 3.));
    }

    SECTION("Alpha = 2")
    {
        SphereEnvelope envelope{1., 3., 1., 1., 2.};

        REQUIRE(Approx(0.).margin(1e-12) == envelope.density(0.9, 0., 0.));

        REQUIRE(Approx(1.) == envelope.density(1., 0., 0.));
        REQUIRE(Approx(1.) == envelope.density(0., 1., 0.));
        REQUIRE(Approx(1.) == envelope.density(0., 0., 1.));

        REQUIRE(Approx(0.25) == envelope.density(-2., 0., 0.));
        REQUIRE(Approx(0.25) == envelope.density(0., -2., 0.));
        REQUIRE(Approx(0.25) == envelope.density(0., 0., -2.));

        REQUIRE(Approx(1. / 9) == envelope.density(-3., 0., 0.));
        REQUIRE(Approx(1. / 9) == envelope.density(0., -3., 0.));
        REQUIRE(Approx(1. / 9) == envelope.density(0., 0., -3.));

        REQUIRE(Approx(0.).margin(1e-12) == envelope.density(3.1, 0., 0.));
        REQUIRE(Approx(0.).margin(1e-12) == envelope.density(0.0, -3.1, 0.));
        REQUIRE(Approx(0.).margin(1e-12) == envelope.density(-1., 0., 3.));
    }
}
