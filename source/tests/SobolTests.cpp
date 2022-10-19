#include "../third-party/catch2/catch.hpp"
#include "../Sobol.hpp"

#include <sstream>

TEST_CASE("Sobol", "[Sobol]")
{
    SECTION("2D")
    {
        unsigned int const NUM_POINTS{ 25 };
        unsigned int const NUM_DIMENSIONS{ 2 };

        double const POINT_TABLE[NUM_POINTS * NUM_DIMENSIONS] = {
            0.5,      0.5,
            0.75,     0.25,
            0.25,     0.75,
            0.375,    0.375,
            0.875,    0.875,
            0.625,    0.125,
            0.125,    0.625,
            0.1875,   0.3125,
            0.6875,   0.8125,
            0.9375,   0.0625,
            0.4375,   0.5625,
            0.3125,   0.1875,
            0.8125,   0.6875,
            0.5625,   0.4375,
            0.0625,   0.9375,
            0.09375,  0.46875,
            0.59375,  0.96875,
            0.84375,  0.21875,
            0.34375,  0.71875,
            0.46875,  0.09375,
            0.96875,  0.59375,
            0.71875,  0.34375,
            0.21875,  0.84375,
            0.15625,  0.15625,
            0.65625,  0.65625
        };

        Sobol generator(NUM_DIMENSIONS);

        for (unsigned int pointId = 0; pointId != NUM_POINTS; ++pointId)
        {
            for (unsigned int dim = 0; dim != NUM_DIMENSIONS; ++dim)
            {
                REQUIRE(POINT_TABLE[pointId*NUM_DIMENSIONS + dim] == generator.Get());
            }
        }
    }

    SECTION("3D")
    {
        unsigned int const NUM_POINTS{ 20 };
        unsigned int const NUM_DIMENSIONS{ 3 };

        double const POINT_TABLE[NUM_POINTS * NUM_DIMENSIONS] = {
            0.5,      0.5,      0.5,
            0.75,     0.25,     0.25,
            0.25,     0.75,     0.75,
            0.375,    0.375,    0.625,
            0.875,    0.875,    0.125,
            0.625,    0.125,    0.875,
            0.125,    0.625,    0.375,
            0.1875,   0.3125,   0.9375,
            0.6875,   0.8125,   0.4375,
            0.9375,   0.0625,   0.6875,
            0.4375,   0.5625,   0.1875,
            0.3125,   0.1875,   0.3125,
            0.8125,   0.6875,   0.8125,
            0.5625,   0.4375,   0.0625,
            0.0625,   0.9375,   0.5625,
            0.09375,  0.46875,  0.46875,
            0.59375,  0.96875,  0.96875,
            0.84375,  0.21875,  0.21875,
            0.34375,  0.71875,  0.71875,
            0.46875,  0.09375,  0.84375
        };

        Sobol generator(NUM_DIMENSIONS);

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
        Sobol longGenerator(2);
        Sobol firstGenerator(2);

        for (int i=0; i!=10; ++i)
        {
            (void)longGenerator.Get();
            (void)firstGenerator.Get();
        }

        std::stringstream stream;
        firstGenerator.save(stream);

        Sobol secondGenerator(2);
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
            0.5, 	    0.5, 	    0.5, 	    0.5, 	    0.5, 	    0.5,
            0.75, 	    0.25, 	    0.25, 	    0.25, 	    0.75, 	    0.75,
            0.25, 	    0.75, 	    0.75, 	    0.75, 	    0.25, 	    0.25,
            0.375, 	    0.375, 	    0.625, 	    0.875, 	    0.375, 	    0.125,
            0.875, 	    0.875, 	    0.125, 	    0.375, 	    0.875, 	    0.625,
            0.625, 	    0.125, 	    0.875, 	    0.625, 	    0.625, 	    0.875,
            0.125, 	    0.625, 	    0.375, 	    0.125, 	    0.125, 	    0.375,
            0.1875, 	0.3125, 	0.9375, 	0.4375, 	0.5625, 	0.3125,
            0.6875, 	0.8125, 	0.4375, 	0.9375, 	0.0625, 	0.8125,
            0.9375, 	0.0625, 	0.6875, 	0.1875, 	0.3125, 	0.5625
        };

        Sobol generator(NUM_DIMENSIONS);

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
