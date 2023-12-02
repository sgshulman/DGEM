#include "../third-party/catch2/catch.hpp"
#include "../Hammersley.hpp"

#include <sstream>

TEST_CASE("Hammersley", "[Hammersley]")
{
    SECTION("2D")
    {
        unsigned int const NUM_POINTS{ 30 };
        unsigned int const NUM_DIMENSIONS{ 2 };

        double const POINT_TABLE[NUM_POINTS * NUM_DIMENSIONS] = {
             1./60,   1./ 2,
             3./60,   1./ 4,
             5./60,   3./ 4,
             7./60,   1./ 8,
             9./60,   5./ 8,
            11./60,   3./ 8,
            13./60,   7./ 8,
            15./60,   1./16,
            17./60,   9./16,
            19./60,   5./16,
            21./60,  13./16,
            23./60,   3./16,
            25./60,  11./16,
            27./60,   7./16,
            29./60,  15./16,
            31./60,   1./32,
            33./60,  17./32,
            35./60,   9./32,
            37./60,  25./32,
            39./60,   5./32,
            41./60,  21./32,
            43./60,  13./32,
            45./60,  29./32,
            47./60,   3./32,
            49./60,  19./32,
            51./60,  11./32,
            53./60,  27./32,
            55./60,   7./32,
            57./60,  23./32,
            59./60,  15./32,
        };

        Hammersley generator(NUM_DIMENSIONS, NUM_POINTS);

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
             1./40,     1./ 2,     1./ 3,     1./ 5,     1./ 7,
             3./40,     1./ 4,     2./ 3,     2./ 5,     2./ 7,
             5./40,     3./ 4,     1./ 9,     3./ 5,     3./ 7,
             7./40,     1./ 8,     4./ 9,     4./ 5,     4./ 7,
             9./40,     5./ 8,     7./ 9,     1./25,     5./ 7,
            11./40,     3./ 8,     2./ 9,     6./25,     6./ 7,
            13./40,     7./ 8,     5./ 9,    11./25,     1./49,
            15./40,     1./16,     8./ 9,    16./25,     8./49,
            17./40,     9./16,     1./27,    21./25,    15./49,
            19./40,     5./16,    10./27,     2./25,    22./49,
            21./40,    13./16,    19./27,     7./25,    29./49,
            23./40,     3./16,     4./27,    12./25,    36./49,
            25./40,    11./16,    13./27,    17./25,    43./49,
            27./40,     7./16,    22./27,    22./25,     2./49,
            29./40,    15./16,     7./27,     3./25,     9./49,
            31./40,     1./32,    16./27,     8./25,    16./49,
            33./40,    17./32,    25./27,    13./25,    23./49,
            35./40,     9./32,     2./27,    18./25,    30./49,
            37./40,    25./32,    11./27,    23./25,    37./49,
            39./40,     5./32,    20./27,     4./25,    44./49,
        };

        Hammersley generator(NUM_DIMENSIONS, NUM_POINTS);

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
        Hammersley longGenerator(2, 100);
        Hammersley firstGenerator(2, 100);

        for (int i=0; i!=10; ++i)
        {
            (void)longGenerator.Get();
            (void)firstGenerator.Get();
        }

        std::stringstream stream;
        firstGenerator.save(stream);

        Hammersley secondGenerator(2, 100);
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
             1./20,     1./ 2,     1./ 3,    1./ 5,    1./ 7,    1./11,
             3./20,     1./ 4,     2./ 3,    2./ 5,    2./ 7,    2./11,
             5./20,     3./ 4,     1./ 9,    3./ 5,    3./ 7,    3./11,
             7./20,     1./ 8,     4./ 9,    4./ 5,    4./ 7,    4./11,
             9./20,     5./ 8,     7./ 9,    1./25,    5./ 7,    5./11,
            11./20,     3./ 8,     2./ 9,    6./25,    6./ 7,    6./11,
            13./20,     7./ 8,     5./ 9,   11./25,    1./49,    7./11,
            15./20,     1./16,     8./ 9,   16./25,    8./49,    8./11,
            17./20,     9./16,     1./27,   21./25,   15./49,    9./11,
            19./20,     5./16,    10./27,    2./25,   22./49,   10./11,
        };

        Hammersley generator(NUM_DIMENSIONS, NUM_POINTS);

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

    SECTION("Configuration")
    {
        unsigned int const NUM_DIMENSIONS{ 5 };
        Hammersley generator(NUM_DIMENSIONS, 100);
        REQUIRE(generator.GetConfiguration() == "Hammersley. dimension = 5. Number of points = 100");
    };
}
