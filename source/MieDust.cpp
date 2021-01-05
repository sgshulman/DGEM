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
    double normFactor = 0.;

    for (std::size_t i = 0; i != table_.size()-1; ++i)
    {
        double const prevH = table_.at(i).cosTheta;
        double const nextH = table_.at(i + 1).cosTheta;
        normFactor += PI * (prevH - nextH) * (table_.at(i + 1).p1 + table_.at(i).p1);
    }

    for (std::size_t i = 0; i != table_.size(); ++i)
    {
        table_.at(i).p1 /= normFactor;
        table_.at(i).p2 /= normFactor;
        table_.at(i).p3 /= normFactor;
        table_.at(i).p4 /= normFactor;
    }

    // accumulated queue
    accumulated_.at(0) = 0.;
    for (std::size_t i = 0; i != table_.size() - 1; ++i)
    {
        double const prevH = table_.at(i).cosTheta;
        double const nextH = table_.at(i + 1).cosTheta;

        accumulated_.at(i + 1) =
            accumulated_.at(i) + PI * (prevH - nextH) * (table_.at(i + 1).p1 + table_.at(i).p1);
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
    auto const it = std::lower_bound(accumulated_.begin(), accumulated_.end(), 2*v - 1.);

    if (it == accumulated_.begin())
    {
        return 1.0;
    }

    if (it == accumulated_.end() || it == accumulated_.end() - 1)
    {
        return -1.0;
    }

    auto const idx = std::distance(it, accumulated_.begin());
    double const delta = accumulated_.at(idx) - accumulated_.at(idx + 1);
    double const wPrev = (accumulated_.at(idx+1) - v) / delta;
    double const wNext = (v - accumulated_.at(idx)) / delta;

    return table_.at(idx).cosTheta * wPrev + table_.at(idx + 1).cosTheta * wNext;
}

