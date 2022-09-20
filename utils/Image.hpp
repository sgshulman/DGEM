#ifndef UTILS_IMAGE_HPP_
#define UTILS_IMAGE_HPP_

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

        double operator[](std::size_t const i) const
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


inline bool readImages(std::vector<Image>& images, int n, char *filenames[])
{
    if (n == 0)
    {
        return true;
    }

    images.reserve(n);

    for (int i=0; i != n; ++i)
    {
        images.emplace_back(filenames[i]);
        if (!images.back().correct())
        {
            return false;
        }
    }

    std::int64_t const rows = images[0].rows();
    std::int64_t const cols = images[0].cols();

    for (std::size_t i=1; i != images.size(); ++i)
    {
        if (images[i].rows() != rows || images[i].cols() != cols)
        {
            std::cerr << filenames[0] << " (" << rows << "x" << cols << ") and " << filenames[i] << " ("
                    << images[i].rows() <<"x" << images[i].cols() << ") should have the same shape." << std::endl;

            return false;
        }
    }

    return true;
}

#endif
