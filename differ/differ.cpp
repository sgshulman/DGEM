#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

class Image
{
    public:
        Image(const char* filename)
        {
            std::ifstream stream(filename);

            if (stream.is_open())
            {
                readFromStream(stream);
            } else {
                std::cerr << "File " << filename << " should exist." << std::endl;
            }
        }

        Image(std::istream& stream)
        {
            readFromStream(stream);
        }

        std::int64_t rows() const
        {
            return nRows_;
        }

        std::int64_t cols() const
        {
            return nCols_;
        }

        bool correct() const
        {
            return correct_;
        }

        double operator[](std::size_t i)
        {
            return i < data_.size() ? data_[i] : 0.0;
        }

    private:
        void readFromStream(std::istream& stream)
        {
            std::string line;

            while (std::getline(stream, line))
            {
                std::int64_t col{ 0 };
                std::istringstream iss(line);

                double v;
                while (iss >> v)
                {
                    data_.push_back(v);
                    ++col;
                }

                if (col > 0)
                {
                    if (nCols_ == -1)
                    {
                        nCols_ = col;
                    } else if (nCols_ != col) {
                        std::cerr << "Different number of values in rows." << std::endl;
                        return;
                    }

                    ++nRows_;
                }
            }
            correct_ = true;
        }

        std::int64_t nRows_{ 0 };
        std::int64_t nCols_{ -1 };
        bool correct_{ false };
        std::vector<double> data_;
};


int main( int argc, char *argv[] )
{
	if (argc != 3)
	{
	    std::cerr << argv[0] << " requires 2 arguments." << std::endl;
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

    std::ofstream difstream("dif.dat");
	std::ofstream absdifstream("dif_abs.dat");
	std::ofstream reldifstream("dif_rel.dat");

    difstream.precision(14);
    absdifstream.precision(14);
    reldifstream.precision(14);

	double absdifsum=0.0, reldifsum=0.0, totsum=0;

	int num=0;
	std::size_t i=0;

	for (std::int64_t x=0; x!=image1.rows(); ++x)
	{
		for (std::int64_t y=0; y!=image1.cols(); ++y)
		{
			double const dif = image1[i] - image2[i];
			double const absdif = std::abs(dif);
            double const mean = 0.5*(image1[i] + image2[i]);
			difstream << dif << "\t";
			absdifstream << absdif << "\t";

			if (mean < 0.07)
			{
				absdifsum += absdif;
				totsum += mean;
			}

			if (image1[i] > 0 && image2[i] > 0)
			{
				++num;
				reldifsum += absdif/mean;
				reldifstream << absdif/mean << "\t";
			} else {
				reldifstream << "0\t";
			}
			++i;
		}
        difstream << std::endl;
		absdifstream << std::endl;
		reldifstream << std::endl;
	}
	std::cout << "dif= " << absdifsum << " meandif= " << absdifsum/num << " meanreldif= "
	        << reldifsum/num << " n= " << num << std::endl << "refsum= " << absdifsum/totsum << std::endl;
}
