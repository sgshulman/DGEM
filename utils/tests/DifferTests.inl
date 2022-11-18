#include "../../source/third-party/catch2/catch.hpp"
#include <sstream>

TEST_CASE("Differ", "[differ]")
{
    SECTION("Positive values")
    {
        std::stringstream imageStream1;
        imageStream1 << 0.1 << "\t" << 0.2 << "\t" << 0.3 << "\t" << 0.4 << std::endl;

        std::stringstream imageStream2;
        imageStream2 << 0.1 << "\t" << 0.3 << "\t" << 0.1 << "\t" << 0.5 << std::endl;

        std::stringstream difStream;
        std::stringstream absDifStream;
        std::stringstream relDifStream;

        Image image1(imageStream1);
        Image image2(imageStream2);

        REQUIRE(image1.correct());
        REQUIRE(image2.correct());
        REQUIRE(image1.cols() == 4);
        REQUIRE(image2.cols() == 4);
        REQUIRE(image1.rows() == 1);
        REQUIRE(image2.rows() == 1);

        DifferenceStats stats = computeDifference(image1, image2, difStream, absDifStream, relDifStream);

        REQUIRE(Approx(stats.absDifSum) == 0.4);
        REQUIRE(Approx(stats.relDifSum) == 1.4 + 1 / 4.5);
        REQUIRE(Approx(stats.totalSum) == 1.0);
        REQUIRE(stats.pixelNum == 4);

        Image difImage(difStream);
        REQUIRE(Approx(difImage[0]).margin(1e-14) == 0.0);
        REQUIRE(Approx(difImage[1]) == -0.1);
        REQUIRE(Approx(difImage[2]) ==  0.2);
        REQUIRE(Approx(difImage[3]) == -0.1);

        Image absDifImage(absDifStream);
        REQUIRE(Approx(absDifImage[0]).margin(1e-14) == 0.0);
        REQUIRE(Approx(absDifImage[1]) == 0.1);
        REQUIRE(Approx(absDifImage[2]) == 0.2);
        REQUIRE(Approx(absDifImage[3]) == 0.1);

        Image relDifImage(relDifStream);
        REQUIRE(Approx(relDifImage[0]).margin(1e-14) == 0.0);
        REQUIRE(Approx(relDifImage[1]) == 0.4);
        REQUIRE(Approx(relDifImage[2]) == 1.0);
        REQUIRE(Approx(relDifImage[3]) == 1 / 4.5);

        DifferenceStats stats2 = computeDifference(image1, image2);
        REQUIRE(Approx(stats.absDifSum) == stats2.absDifSum);
        REQUIRE(Approx(stats.relDifSum) == stats2.relDifSum);
        REQUIRE(Approx(stats.totalSum) == stats2.totalSum);
        REQUIRE(stats.pixelNum == stats2.pixelNum);
    }

    SECTION("Negative values")
    {
        std::stringstream imageStream1;
        imageStream1 << 0.0 << "\t" <<  0.1 << std::endl << -0.2 << "\t" << 0.3 << std::endl;

        std::stringstream imageStream2;
        imageStream2 << 0.0 << "\t" << -0.1 << std::endl << -0.1 << "\t" << 0.4 << std::endl;

        std::stringstream difStream;
        std::stringstream absDifStream;
        std::stringstream relDifStream;

        Image image1(imageStream1);
        Image image2(imageStream2);

        REQUIRE(image1.correct());
        REQUIRE(image2.correct());
        REQUIRE(image1.cols() == 2);
        REQUIRE(image2.cols() == 2);
        REQUIRE(image1.rows() == 2);
        REQUIRE(image2.rows() == 2);

        DifferenceStats stats = computeDifference(image1, image2, difStream, absDifStream, relDifStream);

        REQUIRE(Approx(stats.absDifSum) == 0.4);
        REQUIRE(Approx(stats.relDifSum) == 2 + 1 / 1.5 + 1 / 3.5);
        REQUIRE(Approx(stats.totalSum) == 0.6);
        REQUIRE(stats.pixelNum == 3);

        Image difImage(difStream);
        REQUIRE(Approx(difImage[0]).margin(1e-14) == 0.0);
        REQUIRE(Approx(difImage[1]) ==  0.2);
        REQUIRE(Approx(difImage[2]) == -0.1);
        REQUIRE(Approx(difImage[3]) == -0.1);

        Image absDifImage(absDifStream);
        REQUIRE(Approx(absDifImage[0]).margin(1e-14) == 0.0);
        REQUIRE(Approx(absDifImage[1]) == 0.2);
        REQUIRE(Approx(absDifImage[2]) == 0.1);
        REQUIRE(Approx(absDifImage[3]) == 0.1);

        Image relDifImage(relDifStream);
        REQUIRE(Approx(relDifImage[0]).margin(1e-14) == 0.0);
        REQUIRE(Approx(relDifImage[1]) == 2.0);
        REQUIRE(Approx(relDifImage[2]) == 1 / 1.5);
        REQUIRE(Approx(relDifImage[3]) == 1 / 3.5);

        DifferenceStats stats2 = computeDifference(image1, image2);
        REQUIRE(Approx(stats.absDifSum) == stats2.absDifSum);
        REQUIRE(Approx(stats.relDifSum) == stats2.relDifSum);
        REQUIRE(Approx(stats.totalSum) == stats2.totalSum);
        REQUIRE(stats.pixelNum == stats2.pixelNum);
    }
}
