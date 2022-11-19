#include <fstream>
#include <iostream>
#include <vector>
#include "Image.hpp"

namespace
{
    void computeMean(
        std::vector<Image> const& images,
        std::ostream& sumStream)
    {
        std::int64_t const rows = images[0].rows();
        std::int64_t const cols = images[0].cols();
        std::size_t i=0;
        double weight = 1.0 / images.size();

        for (std::int64_t x=0; x!=rows; ++x)
        {
            for (std::int64_t y=0; y!=cols; ++y)
            {
                double sum = 0.0;

                for (const auto& image : images)
                {
                    sum += image[i];
                }
                sum *= weight;
                sumStream << sum << "\t";

                ++i;
            }
            sumStream << std::endl;
        }
    }
}

int main(int argc, char *argv[])
{
    if (argc > 1 && (std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help"))
    {
        std::cout << "Usage: " << argv[0] << " [OPTION] FILES -o MEANFILE\n"
            << "Compute mean of images from FILES and store the resulting image into MEANFILE.\n"
            << "Options:\n"
            << "\t-h, --help\tdisplay this help and exit.\n";

        return 0;
    }

    if (argc < 4)
    {
        std::cerr << argv[0] << " requires at least one input image filename and -o MEANFILE as arguments." << std::endl;
        return 0;
    }

    if (std::string(argv[argc-2]) != "-o")
    {
        std::cerr << argv[0] << ": -o MEANFILE are required last arguments." << std::endl;
        return 0;
    }

    std::vector<Image> images;

    if (!readImages(images, argc - 3, argv + 1))
    {
        return 0;
    }

    std::ofstream sumStream(argv[argc-1]);
    sumStream.precision(14);

    computeMean(images, sumStream);
    return 0;
}

