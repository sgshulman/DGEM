#include "../third-party/catch2/catch.hpp"

TEST_CASE("Parse Observers", "[model]")
{
    std::vector<Observer> observers;

    SECTION("Manual Positions")
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
                "theta": -60.0
              }
            ]}
            )"_json;

        parseObservers(&observers, observersJson);

        REQUIRE(2 == observers.size());
        REQUIRE(Approx(0.) == observers.at(0).phi());
        REQUIRE(Approx(PI / 4.) == observers.at(0).theta());

        REQUIRE(Approx(PI / 2.) == observers.at(1).phi());
        REQUIRE(Approx(-PI / 3.) == observers.at(1).theta());
    }
}