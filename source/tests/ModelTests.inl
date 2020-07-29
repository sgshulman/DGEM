#include "../third-party/catch2/catch.hpp"

namespace
{
    void compareDirection(const Direction3d& direction, const Vector3d& reference)
    {
        REQUIRE(Approx(1.) == reference.norm());
        REQUIRE(Approx(1.) == direction.vector().norm());

        REQUIRE(Approx(reference.x()).margin(1e-12) == direction.x());
        REQUIRE(Approx(reference.y()).margin(1e-12) == direction.y());
        REQUIRE(Approx(reference.z()).margin(1e-12) == direction.z());
    }
}

TEST_CASE("Parse Observers", "[model]")
{
    std::vector<Observer> observers;

    SECTION("Manual")
    {
        nlohmann::json observersJson = R"({
            "rimage": 800.0,
            "manual": [
              {
                "phi": 0.0,
                "theta": 45.0
              },
              {
                "phi": 90.0,
                "theta": 120.0
              }
            ]}
            )"_json;

        parseObservers(&observers, observersJson);

        REQUIRE(2 == observers.size());
        REQUIRE(Approx(0.).margin(1e-12) == observers.at(0).phi());
        REQUIRE(Approx(PI / 4.) == observers.at(0).theta());

        compareDirection(
            observers.at(0).direction(),
            Vector3d{std::sqrt(2.)/2, 0, std::sqrt(2.)/2});

        REQUIRE(Approx(PI / 2.) == observers.at(1).phi());
        REQUIRE(Approx(2 * PI / 3.) == observers.at(1).theta());
        compareDirection(observers.at(1).direction(), Vector3d{0.,  std::sqrt(3.)/2, -0.5});
    }

    SECTION("Parallel")
    {
        nlohmann::json observersJson = R"({
            "rimage": 800.0,
            "parallel": {
              "numberOfObservers" : 10,
              "theta" : 90.0
            }}
            )"_json;

        parseObservers(&observers, observersJson);

        REQUIRE(10 == observers.size());

        for (std::size_t i = 0; i != observers.size(); ++i)
        {
            REQUIRE(Approx(0.2 * PI * i).margin(1e-12) == observers.at(i).phi());
            REQUIRE(Approx(PI / 2.) == observers.at(i).theta());
        }
    }

    SECTION("Meridian")
    {
        nlohmann::json observersJson = R"({
            "rimage": 800.0,
            "meridian": {
              "numberOfObservers" : 5,
              "phi" : 0.0
            }}
            )"_json;

        parseObservers(&observers, observersJson);

        REQUIRE(5 == observers.size());

        REQUIRE(Approx(PI / 10.) == observers.at(0).theta());
        REQUIRE(Approx(3 * PI / 10.) == observers.at(1).theta());
        REQUIRE(Approx(PI / 2.) == observers.at(2).theta());
        REQUIRE(Approx(7 * PI / 10.) == observers.at(3).theta());
        REQUIRE(Approx(9 * PI / 10.) == observers.at(4).theta());
    }

    SECTION("Combined")
    {
        nlohmann::json observersJson = R"({
            "rimage": 800.0,
            "manual": [
            {
                "phi": 45.0,
                "theta": 45.0
            },
            {
                "phi": 60.0,
                "theta": 120.0
            }
            ],
            "parallel": {
                "numberOfObservers" : 10,
                "theta" : 90.0
            },
            "meridian": {
                "numberOfObservers" : 2,
                "phi" : 0.0
            }}
            )"_json;

        parseObservers(&observers, observersJson);

        REQUIRE(14 == observers.size());
    }
}

TEST_CASE("Parse Geometry", "[model]")
{
    SECTION("Flared Disk")
    {
        nlohmann::json flaredDiskJson = R"({
            "flaredDisk": {
              "rInner": 200,
              "rOuter": 800,
              "rho0": 2.4e-15,
              "h0": 7,
              "r0": 100,
              "alpha": 2.25,
              "beta": 1.25
            }}
            )"_json;

        const auto disk = parseGeometry(flaredDiskJson);
        REQUIRE(disk);
    }

    SECTION("Sphere Envelope")
    {
        nlohmann::json sphereEnvelopJson = R"({
            "sphereEnvelope": {
              "rInner": 500,
              "rOuter": 550,
              "rho0": 2.4e-17,
              "r0": 100,
              "alpha": 0
            }}
            )"_json;

        const auto sphere = parseGeometry(sphereEnvelopJson);
        REQUIRE(sphere);
    }

    SECTION("Maximum")
    {
        nlohmann::json matterArrayJson = R"({
            "max" : [
              {"sphereEnvelope": {
                "rInner": 500,
                "rOuter": 550,
                "rho0": 2.4e-17,
                "r0": 100,
                "alpha": 0
              }},
              {"flaredDisk": {
                "rInner": 200,
                "rOuter": 800,
                "rho0": 2.4e-15,
                "h0": 7,
                "r0": 100,
                "alpha": 2.25,
                "beta": 1.25
              }}
            ]}
            )"_json;

        const auto matterArray = parseGeometry(matterArrayJson);
        REQUIRE(matterArray);
    }

    SECTION("Flared Disk with the wind and translation")
    {
        nlohmann::json flaredDiskJson = R"({
            "flaredDisk": {
              "rInner": 0.1,
              "rOuter": 26,
              "rho0": 11.945e-8,
              "h0": 0.035,
              "r0": 1,
              "alpha": 3.04,
              "beta": 1.25,
              "safierWind" : {
                "model" : "C",
                "mOut" : 1e-8,
                "mStar" : 2.4,
                "h0" : 0.1,
                "rMin" : 0.1,
                "rMax" : 1
              },
              "translation" : {
                "precession" : 75,
                "nutation" : 15,
                "x" : 1,
                "y" : 0.5
              }
            }}
            )"_json;

        const auto disk = parseGeometry(flaredDiskJson);
        REQUIRE(disk);
    }
}
