#include "../third-party/catch2/catch.hpp"
#include "../Random.hpp"

#include <sstream>

TEST_CASE("Random", "[Saving]")
{
    Random longRandom(-1556);
    Random firstRandom(-1556);

    for (int i=0; i!=10; ++i)
    {
        (void)longRandom.Get();
        (void)firstRandom.Get();
    }

    std::stringstream stream;
    firstRandom.save(stream);

    Random secondRandom(-1556);
    secondRandom.load(stream);

    for (int i=0; i!=10; ++i)
    {
        REQUIRE(Approx(longRandom.Get()) == secondRandom.Get());
    }
}
