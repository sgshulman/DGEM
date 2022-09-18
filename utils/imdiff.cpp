#include <cmath>
#include <fstream>
#include <iostream>

#include "Image.hpp"

struct DifferenceStats
{
    double absDifSum{ 0.0 };
    double relDifSum{ 0.0 };
    double totalSum{ 0.0 };
    int pixelNum{ 0 };
};

DifferenceStats computeDifference(
    Image const& image1,
    Image const& image2,
    std::ostream& difStream,
    std::ostream& absDifStream,
    std::ostream& relDifStream)
{
    DifferenceStats stats;
    std::size_t i=0;
    
    for (std::int64_t x=0; x!=image1.rows(); ++x)
    {
        for (std::int64_t y=0; y!=image1.cols(); ++y)
        {
            double const dif = image1[i] - image2[i];
            double const absDif = std::abs(dif);
            double const mean = 0.5*(image1[i] + image2[i]);
            difStream << dif << "\t";
            absDifStream << absDif << "\t";

            stats.absDifSum += absDif;
            stats.totalSum += mean;

            if (image1[i] != 0 && image2[i] != 0)
            {
                ++stats.pixelNum;
                double const meanAbs = 0.5*(std::abs(image1[i]) + std::abs(image2[i]));
                stats.relDifSum += absDif / meanAbs;
                relDifStream << absDif / meanAbs << "\t";
            } else {
                relDifStream << "0\t";
            }
            ++i;
        }
        difStream << std::endl;
        absDifStream << std::endl;
        relDifStream << std::endl;
    }
    
    return stats;
}

#ifdef ENABLE_UNIT_TESTS
#include "tests/DifferTests.inl"
#else

int main( int argc, char *argv[] )
{
    if (argc > 1 && (std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help"))
    {
        std::cout << "Usage: " << argv[0] << " [OPTION] FILE1 FILE2\n"
            << "Compare image files FILE1 and FILE2 and creates difference files.\n"
            << "Options:\n"
            << "\t-h, --help\tdisplay this help and exit.\n";

        return 0;
    }

    if (argc != 3)
    {
        std::cerr << argv[0] << " requires two image filenames as arguments." << std::endl;
        return 0;
    }

    Image image1(argv[1]);
    Image image2(argv[2]);

    if (!image1.correct() || !image2.correct())
    {
        return 0;
    }

    if (image1.rows() != image2.rows() || image1.cols() != image2.cols())
    {
        std::cerr << argv[1] << " (" << image1.rows() << "x" << image1.cols() << ") and " << argv[2] << " ("
                << image2.rows() <<"x" << image2.cols() << ") should have the same shape." << std::endl;

        return 0;
    }

    std::ofstream difStream("dif.dat");
    std::ofstream absDifStream("dif_abs.dat");
    std::ofstream relDifStream("dif_rel.dat");

    difStream.precision(14);
    absDifStream.precision(14);
    relDifStream.precision(14);

    DifferenceStats stats = computeDifference(image1, image2, difStream, absDifStream, relDifStream);

    std::cout << "dif= " << stats.absDifSum << " meandif= " << stats.absDifSum/stats.pixelNum << " meanreldif= "
            << stats.relDifSum/stats.pixelNum << " n= " << stats.pixelNum << std::endl << "refsum= "
            << stats.absDifSum/stats.totalSum << std::endl;
}

#endif
