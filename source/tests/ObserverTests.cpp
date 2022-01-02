#include "../third-party/catch2/catch.hpp"
#include "../Observer.hpp"

TEST_CASE("Observer. inFov", "[observer]")
{
    SECTION("Without mask")
    {
        Observer observer(0., 0., 100., 0., 200, 200);

        REQUIRE(observer.inFov(Photon{{  0.,   0.,     0.}, 0, Direction3d{0., 0., 1.}, 1, 0}));
        REQUIRE(observer.inFov(Photon{{  0., -0.1,    10.}, 0, Direction3d{0., 0., 1.}, 1, 0}));
        REQUIRE(observer.inFov(Photon{{  0.,    1.,  100.}, 0, Direction3d{0., 0., 1.}, 1, 0}));
        REQUIRE(observer.inFov(Photon{{ 10.,    0., -200.}, 0, Direction3d{0., 0., 1.}, 1, 0}));
        REQUIRE(observer.inFov(Photon{{  0.,   50.,  250.}, 0, Direction3d{0., 0., 1.}, 1, 0}));
        REQUIRE(observer.inFov(Photon{{-75.,    0., -300.}, 0, Direction3d{0., 0., 1.}, 1, 0}));

        REQUIRE(observer.inFov(Photon{{ 99., 0., 0.}, 0, Direction3d{0., 0., 1.}, 1, 0}));
        REQUIRE(observer.inFov(Photon{{-99., 0., 0.}, 0, Direction3d{0., 0., 1.}, 1, 0}));
        REQUIRE(observer.inFov(Photon{{0.,  99., 0.}, 0, Direction3d{0., 0., 1.}, 1, 0}));
        REQUIRE(observer.inFov(Photon{{0., -99., 0.}, 0, Direction3d{0., 0., 1.}, 1, 0}));

        REQUIRE(!observer.inFov(Photon{{ 101., 0., 0.}, 0, Direction3d{0., 0., 1.}, 1, 0}));
        REQUIRE(!observer.inFov(Photon{{-101., 0., 0.}, 0, Direction3d{0., 0., 1.}, 1, 0}));
        REQUIRE(!observer.inFov(Photon{{0.,  101., 0.}, 0, Direction3d{0., 0., 1.}, 1, 0}));
        REQUIRE(!observer.inFov(Photon{{0., -101., 0.}, 0, Direction3d{0., 0., 1.}, 1, 0}));
    }

    SECTION("With mask")
    {
        Observer observer(0., 0., 100., 10., 200, 200);

        REQUIRE(!observer.inFov(Photon{{  0.,   0.,     0.}, 0, Direction3d{0., 0., 1.}, 1, 0}));

        REQUIRE(!observer.inFov(Photon{{ 9., 0., 0.}, 0, Direction3d{0., 0., 1.}, 1, 0}));
        REQUIRE(!observer.inFov(Photon{{-9., 0., 0.}, 0, Direction3d{0., 0., 1.}, 1, 0}));
        REQUIRE(!observer.inFov(Photon{{0.,  9., 0.}, 0, Direction3d{0., 0., 1.}, 1, 0}));
        REQUIRE(!observer.inFov(Photon{{0., -9., 0.}, 0, Direction3d{0., 0., 1.}, 1, 0}));

        REQUIRE(!observer.inFov(Photon{{ 7.,  7., 0.}, 0, Direction3d{0., 0., 1.}, 1, 0}));
        REQUIRE(!observer.inFov(Photon{{ 7., -7., 0.}, 0, Direction3d{0., 0., 1.}, 1, 0}));
        REQUIRE(!observer.inFov(Photon{{-7.,  7., 0.}, 0, Direction3d{0., 0., 1.}, 1, 0}));
        REQUIRE(!observer.inFov(Photon{{-7., -7., 0.}, 0, Direction3d{0., 0., 1.}, 1, 0}));

        REQUIRE(observer.inFov(Photon{{ 11., 0., 0.}, 0, Direction3d{0., 0., 1.}, 1, 0}));
        REQUIRE(observer.inFov(Photon{{-11., 0., 0.}, 0, Direction3d{0., 0., 1.}, 1, 0}));
        REQUIRE(observer.inFov(Photon{{0.,  11., 0.}, 0, Direction3d{0., 0., 1.}, 1, 0}));
        REQUIRE(observer.inFov(Photon{{0., -11., 0.}, 0, Direction3d{0., 0., 1.}, 1, 0}));

        REQUIRE(observer.inFov(Photon{{ 8.,  8., 0.}, 0, Direction3d{0., 0., 1.}, 1, 0}));
        REQUIRE(observer.inFov(Photon{{ 8., -8., 0.}, 0, Direction3d{0., 0., 1.}, 1, 0}));
        REQUIRE(observer.inFov(Photon{{-8.,  8., 0.}, 0, Direction3d{0., 0., 1.}, 1, 0}));
        REQUIRE(observer.inFov(Photon{{-8., -8., 0.}, 0, Direction3d{0., 0., 1.}, 1, 0}));

        REQUIRE(observer.inFov(Photon{{ 99., 0., 0.}, 0, Direction3d{0., 0., 1.}, 1, 0}));
        REQUIRE(observer.inFov(Photon{{-99., 0., 0.}, 0, Direction3d{0., 0., 1.}, 1, 0}));
        REQUIRE(observer.inFov(Photon{{0.,  99., 0.}, 0, Direction3d{0., 0., 1.}, 1, 0}));
        REQUIRE(observer.inFov(Photon{{0., -99., 0.}, 0, Direction3d{0., 0., 1.}, 1, 0}));

        REQUIRE(!observer.inFov(Photon{{ 101., 0., 0.}, 0, Direction3d{0., 0., 1.}, 1, 0}));
        REQUIRE(!observer.inFov(Photon{{-101., 0., 0.}, 0, Direction3d{0., 0., 1.}, 1, 0}));
        REQUIRE(!observer.inFov(Photon{{0.,  101., 0.}, 0, Direction3d{0., 0., 1.}, 1, 0}));
        REQUIRE(!observer.inFov(Photon{{0., -101., 0.}, 0, Direction3d{0., 0., 1.}, 1, 0}));
    }
}