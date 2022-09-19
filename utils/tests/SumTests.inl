#include "../../source/third-party/catch2/catch.hpp"
#include <sstream>

TEST_CASE("Sum", "[sum]")
{
    SECTION("One file")
    {
        std::stringstream imageStream;
        imageStream << 0.1 << "\t" << 0.2 << std::endl << 0.3 << "\t" << 0.4 << std::endl;

        std::vector<Image> images;
        images.emplace_back(imageStream);

        std::stringstream sumStream;

        computeSum(images, sumStream);

        Image sumImage(sumStream);
        REQUIRE(sumImage.correct());
        REQUIRE(sumImage.cols() == 2);
        REQUIRE(sumImage.rows() == 2);

        REQUIRE(Approx(sumImage[0]) == 0.1);
        REQUIRE(Approx(sumImage[1]) == 0.2);
        REQUIRE(Approx(sumImage[2]) == 0.3);
        REQUIRE(Approx(sumImage[3]) == 0.4);
    }

    SECTION("Three files")
    {
        std::stringstream imageStream1;
        imageStream1 << 0.1 << "\t" << 0.2 << "\t" << 0.0 << std::endl << 0.3 << "\t" << 0.4 << "\t" << 0.5 << std::endl;

        std::stringstream imageStream2;
        imageStream2 << 0.2 << "\t" << 0.1 << "\t" << 0.1 << std::endl << 0.1 << "\t" << -0.4 << "\t" << 0.1 << std::endl;

        std::stringstream imageStream3;
        imageStream3 << 0.1 << "\t" << 0.2 << "\t" << 0.2 << std::endl << 0.2 << "\t" << 0.1 << "\t" << -0.5 << std::endl;

        std::vector<Image> images;
        images.emplace_back(imageStream1);
        images.emplace_back(imageStream2);
        images.emplace_back(imageStream3);

        std::stringstream sumStream;

        computeSum(images, sumStream);

        Image sumImage(sumStream);
        REQUIRE(sumImage.correct());
        REQUIRE(sumImage.cols() == 3);
        REQUIRE(sumImage.rows() == 2);

        REQUIRE(Approx(sumImage[0]) == 0.4);
        REQUIRE(Approx(sumImage[1]) == 0.5);
        REQUIRE(Approx(sumImage[2]) == 0.3);

        REQUIRE(Approx(sumImage[3]) == 0.6);
        REQUIRE(Approx(sumImage[4]) == 0.1);
        REQUIRE(Approx(sumImage[5]) == 0.1);
    }
}
