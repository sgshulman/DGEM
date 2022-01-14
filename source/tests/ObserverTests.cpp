#include "../third-party/catch2/catch.hpp"
#include "../Observer.hpp"

namespace
{
    double luminosity(Observer const& obs, std::uint32_t const nX = 4, std::uint32_t const nY = 4)
    {
        double sum = 0.0;

        for (std::uint32_t i=0; i != nX; ++i)
        {
            for (std::uint32_t j=0; j != nY; ++j)
            {
                sum += obs.totalLuminosity(i, j);
            }
        }

        return sum;
    }
}


TEST_CASE("Observer. bin", "[observer]")
{
    SECTION("Usual photons")
    {
        // align image and world axes
        Observer observer(radians(-90.0), 0., 2., 0., 4, 4);

        observer.bin(Photon{{ -1.5, -1.5, 0.}, 0, Direction3d{0., 0., 1.}, 1.0, 0});
        REQUIRE(Approx(luminosity(observer)) == 1.0);
        REQUIRE(Approx(observer.totalLuminosity(0, 0)) == 1.0);

        observer.bin(Photon{{ 1.5, 1.5, 0.}, 0, Direction3d{0., 0., 1.}, 1.0, 1});
        REQUIRE(Approx(luminosity(observer)) == 2.0);
        REQUIRE(Approx(observer.totalLuminosity(3, 3)) == 1.0);

        observer.bin(Photon{{ -1.5, 1.5, 0.}, 0, Direction3d{0., 0., 1.}, 1.0, 2});
        REQUIRE(Approx(luminosity(observer)) == 3.0);
        REQUIRE(Approx(observer.totalLuminosity(0, 3)) == 1.0);

        observer.bin(Photon{{ 0.5, -0.5, 0.}, 0, Direction3d{0., 0., 1.}, 1.0, 3});
        REQUIRE(Approx(luminosity(observer)) == 4.0);
        REQUIRE(Approx(observer.totalLuminosity(2, 1)) == 1.0);

        // Photons are outside of the image
        observer.bin(Photon{{ 2.1,  2.1, 0.}, 0, Direction3d{0., 0., 1.}, 1.0, 0});
        observer.bin(Photon{{ 2.1, -2.1, 0.}, 0, Direction3d{0., 0., 1.}, 1.0, 1});
        observer.bin(Photon{{-2.1,  2.1, 0.}, 0, Direction3d{0., 0., 1.}, 1.0, 2});
        observer.bin(Photon{{-2.1, -2.1, 0.}, 0, Direction3d{0., 0., 1.}, 1.0, 3});
        observer.bin(Photon{{ 0.0,  2.1, 0.}, 0, Direction3d{0., 0., 1.}, 1.0, 2});
        observer.bin(Photon{{ 0.0, -2.1, 0.}, 0, Direction3d{0., 0., 1.}, 1.0, 3});
        observer.bin(Photon{{ 2.1,  0.0, 0.}, 0, Direction3d{0., 0., 1.}, 1.0, 2});
        observer.bin(Photon{{-2.1,  0.0, 0.}, 0, Direction3d{0., 0., 1.}, 1.0, 3});
        REQUIRE(Approx(luminosity(observer)) == 4.0);
    }

    SECTION("Border photons")
    {
        // align image and world axes
        Observer observer(radians(-90.0), 0., 2., 0., 4, 4);

        // direct photons are not split between pixels
        observer.bin(Photon{{ 0., 0., 0.}, 0, Direction3d{0., 0., 1.}, 1.0, 0});
        REQUIRE(Approx(luminosity(observer)) == 1.0);
        REQUIRE(Approx(observer.totalLuminosity(2, 2)) == 1.0);

        // scattered photons on pixel borders are distributed between pixels
        observer.bin(Photon{{ 0., 0., 0.}, 0, Direction3d{0., 0., 1.}, 1.0, 1});
        REQUIRE(Approx(luminosity(observer)) == 2.0);
        REQUIRE(Approx(observer.totalLuminosity(1, 1)) == 0.25);
        REQUIRE(Approx(observer.totalLuminosity(1, 2)) == 0.25);
        REQUIRE(Approx(observer.totalLuminosity(2, 1)) == 0.25);
        REQUIRE(Approx(observer.totalLuminosity(2, 2)) == 1.25);

        observer.bin(Photon{{ 1., -1.2, 0.}, 0, Direction3d{0., 0., 1.}, 1.0, 2});
        REQUIRE(Approx(luminosity(observer)) == 3.0);
        REQUIRE(Approx(observer.totalLuminosity(2, 0)) == 0.5);
        REQUIRE(Approx(observer.totalLuminosity(3, 0)) == 0.5);

        observer.bin(Photon{{ -1.5, -1, 0.}, 0, Direction3d{0., 0., 1.}, 1.0, 3});
        REQUIRE(Approx(luminosity(observer)) == 4.0);
        REQUIRE(Approx(observer.totalLuminosity(0, 0)) == 0.5);
        REQUIRE(Approx(observer.totalLuminosity(0, 1)) == 0.5);
    }

    SECTION("Section photons")
    {
        // align image and world axes
        Observer observer(radians(-90.0), 0., 2., 0., 4, 4);

        observer.bin(Photon{{ -1.5, -1., 0.}, 0, Direction3d{0., 0., 1.}, 1.0, 0}, {-1.5, -1.1, 0.0}, {-1.5, -0.1, 0.0});
        REQUIRE(Approx(luminosity(observer)) == 1.0);
        REQUIRE(Approx(observer.totalLuminosity(0, 0)) == 0.1);
        REQUIRE(Approx(observer.totalLuminosity(0, 1)) == 0.9);

        observer.bin(Photon{{ 1, -0.5, 0.}, 0, Direction3d{0., 0., 1.}, 1.0, 1}, {0.9, -0.5, 0.0}, {1.3, -0.5, 0.0});
        REQUIRE(Approx(luminosity(observer)) == 2.0);
        REQUIRE(Approx(observer.totalLuminosity(2, 1)) == 0.25);
        REQUIRE(Approx(observer.totalLuminosity(3, 1)) == 0.75);

        observer.bin(Photon{{ 1.0, 1.2, 0.}, 0, Direction3d{0., 0., 1.}, 1.0, 2}, {0.6, 0.8, 0.0}, {1.4, 1.6, 0.0});
        REQUIRE(Approx(luminosity(observer)) == 3.0);
        REQUIRE(Approx(observer.totalLuminosity(2, 2)) == 0.25);
        REQUIRE(Approx(observer.totalLuminosity(2, 3)) == 0.25);
        REQUIRE(Approx(observer.totalLuminosity(3, 3)) == 0.5);

        // One half of the section is outside of the image
        observer.bin(Photon{{ 0.1,  2.0, 0.}, 0, Direction3d{0., 0., 1.}, 1.0, 2}, {0.1, 2.5, 0.0}, {0.1, 1.5, 0.0});
        REQUIRE(Approx(luminosity(observer)) == 3.5);

        observer.bin(Photon{{ 0.1, -2.0, 0.}, 0, Direction3d{0., 0., 1.}, 1.0, 2}, {0.1, -1.5, 0.0}, {0.1, -2.5, 0.0});
        REQUIRE(Approx(luminosity(observer)) == 4.0);

        observer.bin(Photon{{ 2.0,  0.1, 0.}, 0, Direction3d{0., 0., 1.}, 1.0, 2}, {2.5, 0.1, 0.0}, {1.5, 0.1, 0.0});
        REQUIRE(Approx(luminosity(observer)) == 4.5);

        observer.bin(Photon{{-2.0,  0.1, 0.}, 0, Direction3d{0., 0., 1.}, 1.0, 2}, {-1.5, 0.1, 0.0}, {-2.5, 0.1, 0.0});
        REQUIRE(Approx(luminosity(observer)) == 5.0);
    }

    SECTION("Section border photons")
    {
        // align image and world axes
        Observer observer(radians(-90.0), 0., 2., 0., 4, 4);

        observer.bin(Photon{{ -1.0, 0.0, 0.}, 0, Direction3d{0., 0., 1.}, 1.0, 0}, {-1.0, -1.5, 0.0}, {-1.0, 1.0, 0.0});
        REQUIRE(Approx(luminosity(observer)) == 1.0);
        REQUIRE(Approx(observer.totalLuminosity(0, 0)) == 0.1);
        REQUIRE(Approx(observer.totalLuminosity(1, 0)) == 0.1);
        REQUIRE(Approx(observer.totalLuminosity(0, 1)) == 0.2);
        REQUIRE(Approx(observer.totalLuminosity(1, 1)) == 0.2);
        REQUIRE(Approx(observer.totalLuminosity(0, 2)) == 0.2);
        REQUIRE(Approx(observer.totalLuminosity(1, 2)) == 0.2);

        observer.bin(Photon{{ 0.3, 1.0, 0.}, 0, Direction3d{0., 0., 1.}, 1.0, 1}, {0.2, 1.0, 0.0}, {0.4, 1.0, 0.0});
        REQUIRE(Approx(luminosity(observer)) == 2.0);
        REQUIRE(Approx(observer.totalLuminosity(2, 2)) == 0.5);
        REQUIRE(Approx(observer.totalLuminosity(2, 3)) == 0.5);
    }
}


TEST_CASE("Observer. bin mask", "[observer]")
{
    SECTION("Usual photons")
    {
        // align image and world axes
        Observer observer(radians(-90.0), 0.0, 2., 0.5, 4, 4);

        observer.bin(Photon{{ 0.4, 0.4, 0.}, 0, Direction3d{0., 0., 1.}, 1.0, 1});
        REQUIRE(Approx(luminosity(observer)) == 1.0);

        observer.bin(Photon{{ -0.4, 0., 0.}, 0, Direction3d{0., 0., 1.}, 1.0, 1});
        REQUIRE(Approx(luminosity(observer)) == 1.0);

        observer.bin(Photon{{ -0.3, -0.3, 0.}, 0, Direction3d{0., 0., 1.}, 1.0, 1});
        REQUIRE(Approx(luminosity(observer)) == 1.0);

        observer.bin(Photon{{ 0., -0.51, 0.}, 0, Direction3d{0., 0., 1.}, 1.0, 1});
        REQUIRE(Approx(luminosity(observer)) == 2.0);
    }

    SECTION("Section photons")
    {
        // align image and world axes
        Observer observer(radians(-90.0), 0.0, 2., 0.5, 4, 4);

        SECTION("Half - along the center - one side")
        {
            observer.bin(Photon{{ 0.0, 0.5, 0.}, 0, Direction3d{0., 0., 1.}, 1.0, 1}, {0., 0.1, 0.0}, {0.0, 0.9, 0.0});
            REQUIRE(Approx(luminosity(observer)) == 0.5);

            observer.bin(Photon{{ 0.0, -0.5, 0.}, 0, Direction3d{0., 0., 1.}, 1.0, 1}, {0., -0.1, 0.0}, {0.0, -0.9, 0.0});
            REQUIRE(Approx(luminosity(observer)) == 1.0);

            observer.bin(Photon{{ 0.5, 0.0, 0.}, 0, Direction3d{0., 0., 1.}, 1.0, 1}, {0.1, 0.0, 0.0}, {0.9, 0.0, 0.0});
            REQUIRE(Approx(luminosity(observer)) == 1.5);

            observer.bin(Photon{{ -0.5, -0.0, 0.}, 0, Direction3d{0., 0., 1.}, 1.0, 1}, {-0.1, 0.0, 0.0}, {-0.9, 0.0, 0.0});
            REQUIRE(Approx(luminosity(observer)) == 2.0);
        }

        SECTION("Half - along the center - one side - reversed points order")
        {
            observer.bin(Photon{{ 0.0, 0.5, 0.}, 0, Direction3d{0., 0., 1.}, 1.0, 1}, {0.0, 0.9, 0.0}, {0., 0.1, 0.0});
            REQUIRE(Approx(luminosity(observer)) == 0.5);

            observer.bin(Photon{{ 0.0, -0.5, 0.}, 0, Direction3d{0., 0., 1.}, 1.0, 1}, {0.0, -0.9, 0.0}, {0., -0.1, 0.0});
            REQUIRE(Approx(luminosity(observer)) == 1.0);

            observer.bin(Photon{{ 0.5, 0.0, 0.}, 0, Direction3d{0., 0., 1.}, 1.0, 1}, {0.9, 0.0, 0.0}, {0.1, 0.0, 0.0});
            REQUIRE(Approx(luminosity(observer)) == 1.5);

            observer.bin(Photon{{ -0.5, -0.0, 0.}, 0, Direction3d{0., 0., 1.}, 1.0, 1}, {-0.9, 0.0, 0.0}, {-0.1, 0.0, 0.0});
            REQUIRE(Approx(luminosity(observer)) == 2.0);
        }

        SECTION("Half - along the center - two side")
        {
            observer.bin(Photon{{ 0.0, 0.5, 0.}, 0, Direction3d{0., 0., 1.}, 1.0, 1}, {0., -0.1, 0.0}, {0.0, 1.1, 0.0});
            REQUIRE(Approx(luminosity(observer)) == 0.5);

            observer.bin(Photon{{ 0.0, -0.5, 0.}, 0, Direction3d{0., 0., 1.}, 1.0, 1}, {0., 0.1, 0.0}, {0.0, -1.1, 0.0});
            REQUIRE(Approx(luminosity(observer)) == 1.0);

            observer.bin(Photon{{ 0.5, 0.0, 0.}, 0, Direction3d{0., 0., 1.}, 1.0, 1}, {-0.1, 0.0, 0.0}, {1.1, 0.0, 0.0});
            REQUIRE(Approx(luminosity(observer)) == 1.5);

            observer.bin(Photon{{ -0.5, -0.0, 0.}, 0, Direction3d{0., 0., 1.}, 1.0, 1}, {0.1, 0.0, 0.0}, {-1.1, 0.0, 0.0});
            REQUIRE(Approx(luminosity(observer)) == 2.0);

            observer.bin(Photon{{ 0.0, 0.5, 0.}, 0, Direction3d{0., 0., 1.}, 1.0, 1}, {0.0, 1.1, 0.0}, {0., -0.1, 0.0});
            REQUIRE(Approx(luminosity(observer)) == 2.5);
        }

        SECTION("Half - shifted")
        {
            observer.bin(Photon{{ 0.0, 0.5, 0.}, 0, Direction3d{0., 0., 1.}, 1.0, 1}, {0.1, 0.4, 0.0}, {-0.1, 0.6, 0.0});
            REQUIRE(Approx(luminosity(observer)) == 0.5);

            observer.bin(Photon{{ 0.0, 0.5, 0.}, 0, Direction3d{0., 0., 1.}, 1.0, 1}, {-0.1, 0.6, 0.0}, {0.1, 0.4, 0.0});
            REQUIRE(Approx(luminosity(observer)) == 1.0);
        }

        SECTION("Outside - along the center")
        {
            observer.bin(Photon{{ 0.0, 1.5, 0.}, 0, Direction3d{0., 0., 1.}, 1.0, 1}, {0.0, 1.4, 0.0}, {0.0, 1.6, 0.0});
            REQUIRE(Approx(luminosity(observer)) == 1.0);

            observer.bin(Photon{{ 1.5, 1.5, 0.}, 0, Direction3d{0., 0., 1.}, 1.0, 1}, {1.4, 1.4, 0.0}, {1.6, 1.6, 0.0});
            REQUIRE(Approx(luminosity(observer)) == 2.0);
        }

        SECTION("Outside - shifted")
        {
            observer.bin(Photon{{ 0.2, 1.5, 0.}, 0, Direction3d{0., 0., 1.}, 1.0, 1}, {0.2, 1.4, 0.0}, {0.2, 1.6, 0.0});
            REQUIRE(Approx(luminosity(observer)) == 1.0);

            observer.bin(Photon{{ 1.0, 0.9, 0.}, 0, Direction3d{0., 0., 1.}, 1.0, 1}, {1.1, 1.0, 0.0}, {0.9, 0.8, 0.0});
            REQUIRE(Approx(luminosity(observer)) == 2.0);
        }

        SECTION("Outside - far away")
        {
            observer.bin(Photon{{ 1.0, 1.0, 0.}, 0, Direction3d{0., 0., 1.}, 1.0, 1}, {0.9, 1.1, 0.0}, {1.1, 0.9, 0.0});
            REQUIRE(Approx(luminosity(observer)) == 1.0);

            observer.bin(Photon{{ 1.0, -1.0, 0.}, 0, Direction3d{0., 0., 1.}, 1.0, 1}, {1.0, 0.0, 0.0}, {1.0, -2.0, 0.0});
            REQUIRE(Approx(luminosity(observer)) == 2.0);
        }

        SECTION("Inner")
        {
            observer.bin(Photon{{ 0.0, 0.0, 0.}, 0, Direction3d{0., 0., 1.}, 1.0, 1}, {-0.4, 0.0, 0.0}, {0.4, 0.0, 0.0});
            REQUIRE(luminosity(observer) < std::numeric_limits<double>::epsilon());

            observer.bin(Photon{{ 0.0, 0.0, 0.}, 0, Direction3d{0., 0., 1.}, 1.0, 1}, {0.0, -0.4, 0.0}, {0.0, 0.4, 0.0});
            REQUIRE(luminosity(observer) < std::numeric_limits<double>::epsilon());

            observer.bin(Photon{{ 0.1, 0.1, 0.}, 0, Direction3d{0., 0., 1.}, 1.0, 1}, {0.1, 0.05, 0.0}, {0.1, 0.15, 0.0});
            REQUIRE(luminosity(observer) < std::numeric_limits<double>::epsilon());

            observer.bin(Photon{{ 0.2, 0.2, 0.}, 0, Direction3d{0., 0., 1.}, 1.0, 1}, {0.3, 0.1, 0.0}, {0.1, 0.3, 0.0});
            REQUIRE(luminosity(observer) < std::numeric_limits<double>::epsilon());
        }

        SECTION("Single point")
        {
            observer.bin(Photon{{ 0.5, 0.5, 0.}, 0, Direction3d{0., 0., 1.}, 1.0, 1}, {0.5, 0.5, 0.0}, {0.5, 0.5, 0.0});
            REQUIRE(Approx(luminosity(observer)) == 1.0);

            observer.bin(Photon{{ 0.2, 0.2, 0.}, 0, Direction3d{0., 0., 1.}, 1.0, 1}, {0.2, 0.2, 0.0}, {0.2, 0.2, 0.0});
            REQUIRE(Approx(luminosity(observer)) == 1.0);
        }

        SECTION("Long - along the center")
        {
            observer.bin(Photon{{ 0.0, 0.0, 0.}, 0, Direction3d{0., 0., 1.}, 1.0, 1}, {1.0, 0.0, 0.0}, {-1.0, 0.0, 0.0});
            REQUIRE(Approx(luminosity(observer)) == 0.5);

            observer.bin(Photon{{ 0.0, 0.0, 0.}, 0, Direction3d{0., 0., 1.}, 1.0, 1}, {0.0, -0.6, 0.0}, {0.0, 1.4, 0.0});
            REQUIRE(Approx(luminosity(observer)) == 1.0);
        }
    }
}


TEST_CASE("Observer. inFov", "[observer]")
{
    SECTION("Without mask")
    {
        Observer observer(0., 0., 100., 0., 200, 200);

        REQUIRE(observer.inFov({  0.,   0.,    0.}));
        REQUIRE(observer.inFov({  0., -0.1,   10.}));
        REQUIRE(observer.inFov({  0.,   1.,  100.}));
        REQUIRE(observer.inFov({ 10.,   0., -200.}));
        REQUIRE(observer.inFov({  0.,  50.,  250.}));
        REQUIRE(observer.inFov({-75.,   0., -300.}));

        REQUIRE(observer.inFov({ 99., 0., 0.}));
        REQUIRE(observer.inFov({-99., 0., 0.}));
        REQUIRE(observer.inFov({0.,  99., 0.}));
        REQUIRE(observer.inFov({0., -99., 0.}));

        REQUIRE(!observer.inFov({ 101., 0., 0.}));
        REQUIRE(!observer.inFov({-101., 0., 0.}));
        REQUIRE(!observer.inFov({0.,  101., 0.}));
        REQUIRE(!observer.inFov({0., -101., 0.}));
    }

    SECTION("With mask")
    {
        Observer observer(0., 0., 100., 10., 200, 200);

        REQUIRE(!observer.inFov({ 0., 0., 0.}));

        REQUIRE(!observer.inFov({ 9., 0., 0.}));
        REQUIRE(!observer.inFov({-9., 0., 0.}));
        REQUIRE(!observer.inFov({ 0., 9., 0.}));
        REQUIRE(!observer.inFov({ 0.,-9., 0.}));

        REQUIRE(!observer.inFov({ 7.,  7., 0.}));
        REQUIRE(!observer.inFov({ 7., -7., 0.}));
        REQUIRE(!observer.inFov({-7.,  7., 0.}));
        REQUIRE(!observer.inFov({-7., -7., 0.}));

        REQUIRE(observer.inFov({ 11., 0., 0.}));
        REQUIRE(observer.inFov({-11., 0., 0.}));
        REQUIRE(observer.inFov({0.,  11., 0.}));
        REQUIRE(observer.inFov({0., -11., 0.}));

        REQUIRE(observer.inFov({ 8.,  8., 0.}));
        REQUIRE(observer.inFov({ 8., -8., 0.}));
        REQUIRE(observer.inFov({-8.,  8., 0.}));
        REQUIRE(observer.inFov({-8., -8., 0.}));

        REQUIRE(observer.inFov({ 99., 0., 0.}));
        REQUIRE(observer.inFov({-99., 0., 0.}));
        REQUIRE(observer.inFov({0.,  99., 0.}));
        REQUIRE(observer.inFov({0., -99., 0.}));

        REQUIRE(!observer.inFov({ 101., 0., 0.}));
        REQUIRE(!observer.inFov({-101., 0., 0.}));
        REQUIRE(!observer.inFov({0.,  101., 0.}));
        REQUIRE(!observer.inFov({0., -101., 0.}));
    }
}