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

    inline double albedo() const
    {	return albedo_;	}
    inline double hgg() const
    {	return hgg_;	}
    inline double hgg2() const
    {	return hgg2_;	}
    inline double pl() const
    {	return pl_;		}
    inline double pc() const
    {	return pc_;		}
    inline double sc() const
    {	return sc_;		}

private:
    double const albedo_;
    double const hgg_;
    double const hgg2_;
    double const pl_;
    double const pc_;
    double const sc_;
};

#endif