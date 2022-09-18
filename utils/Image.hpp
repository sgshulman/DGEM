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
