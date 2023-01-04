#ifndef MATTER_ARRAY_HPP_
#define MATTER_ARRAY_HPP_

#include <vector>
#include "IMatter.hpp"
#include "Predefines.hpp"

class MatterArray final: public IMatter
{
public:
    using AccumulatorType = double (*)(double, double);

    static double sum(double left, double right)
    {
        return left + right;
    }

    static double max(double left, double right)
    {
        return left > right ? left : right;
    }

    explicit MatterArray(
        std::vector<IMatterCPtr> matterArray,
        AccumulatorType accumulator)
        : matterArray_{ std::move(matterArray) }
        , accumulator_{ accumulator }
    {}

    ~MatterArray() override = default;

    double density(Vector3d const& position) const override;

private:
    std::vector<IMatterCPtr> matterArray_;
    AccumulatorType accumulator_;
};

#endif
