#include "../third-party/catch2/catch.hpp"
#include "../Halton.hpp"

#include <sstream>

TEST_CASE("Halton", "[Halton]")
{
    SECTION("2D")
    {
        unsigned int const NUM_POINTS{ 30 };
        unsigned int const NUM_DIMENSIONS{ 2 };

        double const POINT_TABLE[NUM_POINTS * NUM_DIMENSIONS] = {
             1./ 2,     1./ 3,
             1./ 4,     2./ 3,
             3./ 4,     1./ 9,
             1./ 8,     4./ 9,
             5./ 8,     7./ 9,
             3./ 8,     2./ 9,
             7./ 8,     5./ 9,
             1./16,     8./ 9,
             9./16,     1./27,
             5./16,    10./27,
            13./16,    19./27,
             3./16,     4./27,
            11./16,    13./27,
             7./16,    22./27,
            15./16,     7./27,
             1./32,    16./27,
            17./32,    25./27,
             9./32,     2./27,
            25./32,    11./27,
             5./32,    20./27,
            21./32,     5./27,
            13./32,    14./27,
            29./32,    23./27,
             3./32,     8./27,
            19./32,    17./27,
            11./32,    26./27,
            27./32,     1./81,
             7./32,    28./81,
            23./32,    55./81,
            15./32,    10./81,
        };

        Halton generator(NUM_DIMENSIONS);

        for (unsigned int pointId = 0; pointId != NUM_POINTS; ++pointId)
        {
            for (unsigned int dim = 0; dim != NUM_DIMENSIONS; ++dim)
            {
                REQUIRE(POINT_TABLE[pointId*NUM_DIMENSIONS + dim] == generator.Get());
            }
        }
    }

    SECTION("5D")
    {
        unsigned int const NUM_POINTS{ 20 };
        unsigned int const NUM_DIMENSIONS{ 5 };

        double const POINT_TABLE[NUM_POINTS * NUM_DIMENSIONS] = {
             1./ 2,     1./ 3,     1./ 5,     1./ 7,      1./ 11,
             1./ 4,     2./ 3,     2./ 5,     2./ 7,      2./ 11,
             3./ 4,     1./ 9,     3./ 5,     3./ 7,      3./ 11,
             1./ 8,     4./ 9,     4./ 5,     4./ 7,      4./ 11,
             5./ 8,     7./ 9,     1./25,     5./ 7,      5./ 11,
             3./ 8,     2./ 9,     6./25,     6./ 7,      6./ 11,
             7./ 8,     5./ 9,    11./25,     1./49,      7./ 11,
             1./16,     8./ 9,    16./25,     8./49,      8./ 11,
             9./16,     1./27,    21./25,    15./49,      9./ 11,
             5./16,    10./27,     2./25,    22./49,     10./ 11,
            13./16,    19./27,     7./25,    29./49,      1./121,
             3./16,     4./27,    12./25,    36./49,     12./121,
            11./16,    13./27,    17./25,    43./49,     23./121,
             7./16,    22./27,    22./25,     2./49,     34./121,
            15./16,     7./27,     3./25,     9./49,     45./121,
             1./32,    16./27,     8./25,    16./49,     56./121,
            17./32,    25./27,    13./25,    23./49,     67./121,
             9./32,     2./27,    18./25,    30./49,     78./121,
            25./32,    11./27,    23./25,    37./49,     89./121,
             5./32,    20./27,     4./25,    44./49,    100./121,
        };

        Halton generator(NUM_DIMENSIONS);

        for (unsigned int pointId = 0; pointId != NUM_POINTS; ++pointId)
        {
            for (unsigned int dim = 0; dim != NUM_DIMENSIONS; ++dim)
            {
                REQUIRE(POINT_TABLE[pointId*NUM_DIMENSIONS + dim] == generator.Get());
            }
        }
    }

    SECTION("Saving")
    {
        Halton longGenerator(2);
        Halton firstGenerator(2);

        for (int i=0; i!=10; ++i)
        {
            (void)longGenerator.Get();
            (void)firstGenerator.Get();
        }

        std::stringstream stream;
        firstGenerator.save(stream);

        Halton secondGenerator(2);
        secondGenerator.load(stream);

        for (int i=0; i!=10; ++i)
        {
            REQUIRE(Approx(longGenerator.Get()) == secondGenerator.Get());
        }
    }

    SECTION("Skip")
    {
        unsigned int const NUM_POINTS{ 10 };
        unsigned int const NUM_DIMENSIONS{ 6 };

        double const POINT_TABLE[NUM_POINTS * NUM_DIMENSIONS] = {
            1./ 2,     1./ 3,    1./ 5,    1./ 7,    1./11,     1./13,
            1./ 4,     2./ 3,    2./ 5,    2./ 7,    2./11,     2./13,
            3./ 4,     1./ 9,    3./ 5,    3./ 7,    3./11,     3./13,
            1./ 8,     4./ 9,    4./ 5,    4./ 7,    4./11,     4./13,
            5./ 8,     7./ 9,    1./25,    5./ 7,    5./11,     5./13,
            3./ 8,     2./ 9,    6./25,    6./ 7,    6./11,     6./13,
            7./ 8,     5./ 9,   11./25,    1./49,    7./11,     7./13,
            1./16,     8./ 9,   16./25,    8./49,    8./11,     8./13,
            9./16,     1./27,   21./25,   15./49,    9./11,     9./13,
            5./16,    10./27,    2./25,   22./49,   10./11,    10./13,
        };

        Halton generator(NUM_DIMENSIONS);

        for (unsigned int pointId = 0; pointId != 5; ++pointId)
        {
            for (unsigned int dim = 0; dim != NUM_DIMENSIONS; ++dim)
            {
                REQUIRE(POINT_TABLE[pointId*NUM_DIMENSIONS + dim] == generator.Get());
            }
            generator.Skip();
        }

        REQUIRE(POINT_TABLE[5*NUM_DIMENSIONS + 0] == generator.Get());
        generator.Skip();

        REQUIRE(POINT_TABLE[6*NUM_DIMENSIONS + 0] == generator.Get());
        REQUIRE(POINT_TABLE[6*NUM_DIMENSIONS + 1] == generator.Get());
        generator.Skip();

        REQUIRE(POINT_TABLE[7*NUM_DIMENSIONS + 0] == generator.Get());
        REQUIRE(POINT_TABLE[7*NUM_DIMENSIONS + 1] == generator.Get());
        REQUIRE(POINT_TABLE[7*NUM_DIMENSIONS + 2] == generator.Get());
        generator.Skip();

        REQUIRE(POINT_TABLE[8*NUM_DIMENSIONS + 0] == generator.Get());
        REQUIRE(POINT_TABLE[8*NUM_DIMENSIONS + 1] == generator.Get());
        REQUIRE(POINT_TABLE[8*NUM_DIMENSIONS + 2] == generator.Get());
        REQUIRE(POINT_TABLE[8*NUM_DIMENSIONS + 3] == generator.Get());
        generator.Skip();

        REQUIRE(POINT_TABLE[9*NUM_DIMENSIONS + 0] == generator.Get());
        REQUIRE(POINT_TABLE[9*NUM_DIMENSIONS + 1] == generator.Get());
        REQUIRE(POINT_TABLE[9*NUM_DIMENSIONS + 2] == generator.Get());
        REQUIRE(POINT_TABLE[9*NUM_DIMENSIONS + 3] == generator.Get());
        REQUIRE(POINT_TABLE[9*NUM_DIMENSIONS + 4] == generator.Get());
        generator.Skip();
    }
}
