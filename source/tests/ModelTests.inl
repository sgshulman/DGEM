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
}
