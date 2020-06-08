#ifndef DUST_HPP_
#define DUST_HPP_

class Dust
{
public:
    Dust(
        double const albedo,
        double const hgg,
        double const pl,
        double const pc,
        double const sc)
    : albedo_{ albedo }
    , hgg_{ hgg }
    , hgg2_{ hgg_ * hgg_ }
    , pl_{ pl }
    , pc_{ pc }
    , sc_{ sc }
    {}

    void scatteringMatrixElements(
        double &p1, double &p2, double &p3, double &p4, double cosTheta) const;

    double fraction(double cosTheta) const;
    double cosRandomTheta(double v) const;

    inline double albedo() const
    {	return albedo_;	}

private:
    double const albedo_;
    double const hgg_;
    double const hgg2_;
    double const pl_;
    double const pc_;
    double const sc_;
};

#endif