#include <cmath>
#include <fstream>

#include "MieDust.hpp"
#include "MathUtils.hpp"

MieDust::MieDust(double const albedo, std::string const& tableFile)
    : albedo_{ albedo }
{
    std::ifstream table(tableFile);

    Data entry{};
    double theta, p1, p2;

    while (table >> theta >> p1 >> p2 >> entry.p3 >> entry.p4)
    {
        entry.cosTheta = cos(radians(theta));
        entry.p1 = 0.5 * (p1 + p2);
        entry.p2 = 0.5 * (p1 - p2);
        table_.push_back(entry);
    }

    // Normalize coefficients
    double lastH = 1.;
    double normFactor = 0.;

    for (std::size_t i = 0; i != table_.size()-1; ++i)
    {
        double const newH = 0.5 * (table_.at(i + 1).cosTheta + table_.at(i).cosTheta);
        normFactor += 2 * PI * (lastH - newH) * table_.at(i).p1;
        lastH = newH;
    }

    normFactor += 2 * PI * (lastH + 1.0);

    for (std::size_t i = 0; i != table_.size(); ++i)
    {
        table_.at(i).p1 /= normFactor;
        table_.at(i).p2 /= normFactor;
        table_.at(i).p3 /= normFactor;
        table_.at(i).p4 /= normFactor;
    }
}


void MieDust::scatteringMatrixElements(
    double &p1, double &p2, double &p3, double &p4, double cosTheta) const
{
    auto const it = std::lower_bound(table_.begin(), table_.end(), cosTheta);

    if (it == table_.begin())
    {
        p1 = it->p1;
        p2 = it->p2;
        p3 = it->p3;
        p4 = it->p4;
    }
    else
    {
        auto const prev = it - 1;
        double const delta = prev->cosTheta - it->cosTheta;
        double const wIt = (cosTheta - prev->cosTheta) / delta;
        double const wPrev = (it->cosTheta - cosTheta) / delta;

        p1 = it->p1 * wIt + prev->p1 * wPrev;
        p2 = it->p2 * wIt + prev->p2 * wPrev;
        p3 = it->p3 * wIt + prev->p3 * wPrev;
        p4 = it->p4 * wIt + prev->p4 * wPrev;
    }
}


double MieDust::fraction(double const cosTheta) const
{
    auto const it = std::lower_bound(table_.begin(), table_.end(), cosTheta);

    if (it == table_.begin())
    {
        return it->p1;
    }

    auto const prev = it - 1;
    double const delta = prev->cosTheta - it->cosTheta;
    double const wIt = (cosTheta - prev->cosTheta) / delta;
    double const wPrev = (it->cosTheta - cosTheta) / delta;

    return it->p1 * wIt + prev->p1 * wPrev;
}


double MieDust::cosRandomTheta(double const v) const
{
    (void)v;
    return 0.0;
}

