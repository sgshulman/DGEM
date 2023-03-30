#include "../third-party/catch2/catch.hpp"
#include "../StdRandomGenerator.hpp"

#include <sstream>

namespace
{
    template<typename GeneratorType>
    void TestGeneratorSaving()
    {
        StdRandomGenerator<GeneratorType> longRandom(1U);
        StdRandomGenerator<GeneratorType> firstRandom(1U);

        for (int i=0; i!=10; ++i)
        {
            (void)longRandom.Get();
            (void)firstRandom.Get();
        }

        std::stringstream stream;
        firstRandom.save(stream);

        StdRandomGenerator<GeneratorType> secondRandom(1U);
        secondRandom.load(stream);

        for (int i=0; i!=10; ++i)
        {
            REQUIRE(Approx(longRandom.Get()) == secondRandom.Get());
        }
    }
} // namespace


TEST_CASE("StdRandomGenerator", "[StdRandomGenerator saving]")
{
    SECTION("MinimumStandard")
    {
        TestGeneratorSaving<std::minstd_rand>();
    }

    SECTION("MersenneTwister")
    {
        TestGeneratorSaving<std::mt19937_64>();
    }

    SECTION("Ranlux48")
    {
        TestGeneratorSaving<std::ranlux48>();
    }
}
