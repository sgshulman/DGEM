#include "../third-party/catch2/catch.hpp"
#include "../StdRandomGenerator.hpp"

#include <sstream>

TEST_CASE("StdRandomGenerator", "[StdRandomGenerator saving]")
{
    StdRandomGenerator<std::mt19937_64> longRandom(1U);
    StdRandomGenerator<std::mt19937_64> firstRandom(1U);

    for (int i=0; i!=10; ++i)
    {
        (void)longRandom.Get();
        (void)firstRandom.Get();
    }

    std::stringstream stream;
    firstRandom.save(stream);

    StdRandomGenerator<std::mt19937_64> secondRandom(1U);
    secondRandom.load(stream);

    for (int i=0; i!=10; ++i)
    {
        REQUIRE(Approx(longRandom.Get()) == secondRandom.Get());
    }
}
