#include "../third-party/catch2/catch.hpp"
#include "../LEcuyer.hpp"

#include <sstream>

TEST_CASE("LEcuyer", "[LEcuyer saving]")
{
    LEcuyer longRandom(-1556);
    LEcuyer firstRandom(-1556);

    for (int i=0; i!=10; ++i)
    {
        (void)longRandom.Get();
        (void)firstRandom.Get();
    }

    std::stringstream stream;
    firstRandom.save(stream);

    LEcuyer secondRandom(-1556);
    secondRandom.load(stream);

    for (int i=0; i!=10; ++i)
    {
        REQUIRE(Approx(longRandom.Get()) == secondRandom.Get());
    }
}
