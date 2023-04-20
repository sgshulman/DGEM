#include "../third-party/catch2/catch.hpp"
#include "../Faure.hpp"

#include <sstream>

TEST_CASE("Faure", "[Faure]")
{
    SECTION("2D")
    {
        unsigned int const NUM_POINTS{ 30 };
        unsigned int const NUM_DIMENSIONS{ 2 };

        double const POINT_TABLE[NUM_POINTS * NUM_DIMENSIONS] = {
             1./ 2,    	 1./ 2,
             1./ 4,	     3./ 4,
             3./ 4,	     1./ 4,
             1./ 8,	     5./ 8,
             5./ 8,	     1./ 8,
             3./ 8,	     3./ 8,
             7./ 8,	     7./ 8,
             1./16,	    15./16,
             9./16,	     7./16,
             5./16,	     3./16,
            13./16,	    11./16,
             3./16,	     5./16,
            11./16,	    13./16,
             7./16,	     9./16,
            15./16,	     1./16,
             1./32,	    17./32,
            17./32,	     1./32,
             9./32,	     9./32,
            25./32,	    25./32,
             5./32,	     5./32,
            21./32,	    21./32,
            13./32,	    29./32,
            29./32,	    13./32,
             3./32,	    15./32,
            19./32,	    31./32,
            11./32,     23./32,
            27./32,      7./32,
             7./32,     27./32,
            23./32,     11./32,
            15./32,      3./32
        };

        Faure generator(NUM_DIMENSIONS);

        for (unsigned int pointId = 0; pointId != NUM_POINTS; ++pointId)
        {
            for (unsigned int dim = 0; dim != NUM_DIMENSIONS; ++dim)
            {
                REQUIRE(Approx(POINT_TABLE[pointId*NUM_DIMENSIONS + dim]) == generator.Get());
            }
        }
    }

    SECTION("5D")
    {
        unsigned int const NUM_POINTS{ 24 };
        unsigned int const NUM_DIMENSIONS{ 5 };

        double const POINT_TABLE[NUM_POINTS * NUM_DIMENSIONS] = {
             1./ 5,     1./ 5,     1./ 5,     1./ 5,     1./ 5,
             2./ 5,	    2./ 5,	   2./ 5,     2./ 5,     2./ 5,
             3./ 5,	    3./ 5,	   3./ 5,     3./ 5,     3./ 5,
             4./ 5,	    4./ 5,	   4./ 5,     4./ 5,     4./ 5,

             1./25,	    6./25,	   11./25,   16./25,    21./25,
             6./25,	   11./25,	   16./25,   21./25,     1./25,
            11./25,	   16./25,	   21./25,    1./25,     6./25,
            16./25,	   21./25,	    1./25,    6./25,    11./25,
            21./25,	    1./25,	    6./25,   11./25,    16./25,

             2./25,	   12./25,     22./25,    7./25,    17./25,
             7./25,	   17./25,	    2./25,   12./25,    22./25,
            12./25,	   22./25,	    7./25,   17./25,     2./25,
            17./25,	    2./25,     12./25,   22./25,     7./25,
            22./25,	    7./25,     17./25,    2./25,    12./25,

             3./25,    18./25,      8./25,   23./25,    13./25,
             8./25,    23./25,     13./25,    3./25,    18./25,
            13./25,     3./25,     18./25,    8./25,    23./25,
            18./25,     8./25,     23./25,   13./25,     3./25,
            23./25,    13./25,      3./25,   18./25,     8./25,

             4./25,    24./25,     19./25,   14./25,     9./25,
             9./25,     4./25,     24./25,   19./25,    14./25,
            14./25,     9./25,      4./25,   24./25,    19./25,
            19./25,    14./25,      9./25,    4./25,    24./25,
            24./25,    19./25,     14./25,    9./25,     4./25
        };

        Faure generator(NUM_DIMENSIONS);

        for (unsigned int pointId = 0; pointId != NUM_POINTS; ++pointId)
        {
            for (unsigned int dim = 0; dim != NUM_DIMENSIONS; ++dim)
            {
                REQUIRE(Approx(POINT_TABLE[pointId*NUM_DIMENSIONS + dim]) == generator.Get());
            }
        }
    }

    SECTION("Saving")
    {
        Faure longGenerator(3);
        Faure firstGenerator(3);

        for (int i=0; i!=15; ++i)
        {
            (void)longGenerator.Get();
            (void)firstGenerator.Get();
        }

        std::stringstream stream;
        firstGenerator.save(stream);

        Faure secondGenerator(3);
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
             1./ 7,     1./ 7,     1./ 7,     1./ 7,     1./ 7,     1./ 7,
             2./ 7,     2./ 7,     2./ 7,     2./ 7,     2./ 7,     2./ 7,
             3./ 7,     3./ 7,     3./ 7,     3./ 7,     3./ 7,     3./ 7,
             4./ 7,     4./ 7,     4./ 7,     4./ 7,     4./ 7,     4./ 7,
             5./ 7,     5./ 7,     5./ 7,     5./ 7,     5./ 7,     5./ 7,
             6./ 7,     6./ 7,     6./ 7,     6./ 7,     6./ 7,     6./ 7,
             1./49,     8./49,    15./49,    22./49,    29./49,    36./49,
             8./49,    15./49,    22./49,    29./49, 	36./49,    43./49,
            15./49,    22./49,    29./49,    36./49,    43./49,     1./49,
            22./49,    29./49,    36./49,    43./49,     1./49,     8./49
        };

        Faure generator(NUM_DIMENSIONS);

        for (unsigned int pointId = 0; pointId != 5; ++pointId)
        {
            for (unsigned int dim = 0; dim != NUM_DIMENSIONS; ++dim)
            {
                REQUIRE(Approx(POINT_TABLE[pointId*NUM_DIMENSIONS + dim]) == generator.Get());
            }
            generator.Skip();
        }

        REQUIRE(Approx(POINT_TABLE[5*NUM_DIMENSIONS + 0]) == generator.Get());
        generator.Skip();

        REQUIRE(Approx(POINT_TABLE[6*NUM_DIMENSIONS + 0]) == generator.Get());
        REQUIRE(Approx(POINT_TABLE[6*NUM_DIMENSIONS + 1]) == generator.Get());
        generator.Skip();

        REQUIRE(Approx(POINT_TABLE[7*NUM_DIMENSIONS + 0]) == generator.Get());
        REQUIRE(Approx(POINT_TABLE[7*NUM_DIMENSIONS + 1]) == generator.Get());
        REQUIRE(Approx(POINT_TABLE[7*NUM_DIMENSIONS + 2]) == generator.Get());
        generator.Skip();

        REQUIRE(Approx(POINT_TABLE[8*NUM_DIMENSIONS + 0]) == generator.Get());
        REQUIRE(Approx(POINT_TABLE[8*NUM_DIMENSIONS + 1]) == generator.Get());
        REQUIRE(Approx(POINT_TABLE[8*NUM_DIMENSIONS + 2]) == generator.Get());
        REQUIRE(Approx(POINT_TABLE[8*NUM_DIMENSIONS + 3]) == generator.Get());
        generator.Skip();

        REQUIRE(Approx(POINT_TABLE[9*NUM_DIMENSIONS + 0]) == generator.Get());
        REQUIRE(Approx(POINT_TABLE[9*NUM_DIMENSIONS + 1]) == generator.Get());
        REQUIRE(Approx(POINT_TABLE[9*NUM_DIMENSIONS + 2]) == generator.Get());
        REQUIRE(Approx(POINT_TABLE[9*NUM_DIMENSIONS + 3]) == generator.Get());
        REQUIRE(Approx(POINT_TABLE[9*NUM_DIMENSIONS + 4]) == generator.Get());
        generator.Skip();
    }

    SECTION("Configuration")
    {
        unsigned int const NUM_DIMENSIONS{ 8 };
        Faure generator(NUM_DIMENSIONS);
        REQUIRE(generator.GetConfiguration() == "Faure. dimension = 8");
    };
}