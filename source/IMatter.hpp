#ifndef I_MATTER_HPP_
#define I_MATTER_HPP_

class IMatter
{
    public:
        virtual ~IMatter() = default;
        virtual double density(double x, double y, double z) const = 0;
};

#endif
