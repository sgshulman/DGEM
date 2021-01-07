#include <algorithm>
#include <cmath>
#include <fstream>

#include "MieDust.hpp"
#include "DebugUtils.hpp"
#include "MathUtils.hpp"

MieDust::MieDust(double const albedo, std::string const& tableFile)
    : MieDust(albedo, std::ifstream(tableFile))
{}


MieDust::MieDust(double albedo, std::istream&& stream)
    : albedo_{ albedo }
{
    Data entry{};
    double theta, p1, p2;

    while (stream >> theta >> p1 >> p2 >> entry.p3 >> entry.p4)
    {
        entry.cosTheta = std::cos(radians(theta));
        entry.p1 = 0.5 * (p1 + p2);
        entry.p2 = 0.5 * (p1 - p2);
        table_.push_back(entry);
    }

    DATA_ASSERT(
        std::abs(table_.front().cosTheta - 1.) < std::numeric_limits<float>::epsilon(),
        "First cosine in dust table should be equal to 1.");

    DATA_ASSERT(
        std::abs(table_.back().cosTheta + 1.) < std::numeric_limits<float>::epsilon(),
        "Last cosine in dust table should be equal to -1.");

    DATA_ASSERT(
        std::is_sorted(
            table_.begin(),
            table_.end(),
            [](Data const& left, Data const& right){ return left.cosTheta > right.cosTheta;}),
        "Cosines of dust table entries should be sorted");

    normalize();
    computeAccumulatedFractions();
}


void MieDust::normalize()
{
    double normFactor = 0.;

    for (std::size_t i = 0; i != table_.size() - 1; ++i)
    {
        double const prevH = table_.at(i).cosTheta;
        double const nextH = table_.at(i + 1).cosTheta;
        normFactor += 0.25 * (prevH - nextH) * (table_.at(i + 1).p1 + table_.at(i).p1);
    }

    for (std::size_t i = 0; i != table_.size(); ++i)
    {
        table_.at(i).p1 /= normFactor;
        table_.at(i).p2 /= normFactor;
        table_.at(i).p3 /= normFactor;
        table_.at(i).p4 /= normFactor;
    }
}


void MieDust::computeAccumulatedFractions()
{
    accumulated_.resize(table_.size());
    accumulated_.at(0) = 0.;

    for (std::size_t i = 0; i != table_.size() - 1; ++i)
    {
        double const prevH = table_.at(i).cosTheta;
        double const nextH = table_.at(i + 1).cosTheta;

        accumulated_.at(i + 1) =
            accumulated_.at(i) + 0.25 * (prevH - nextH) * (table_.at(i + 1).p1 + table_.at(i).p1);
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
    } else if (it == table_.end()) {
        auto const& back = table_.back();
        p1 = back.p1;
        p2 = back.p2;
        p3 = back.p3;
        p4 = back.p4;
    } else {
        auto const prev = it - 1;
        double const delta = prev->cosTheta - it->cosTheta;
        double const wIt = (prev->cosTheta - cosTheta) / delta;
        double const wPrev = (cosTheta - it->cosTheta) / delta;

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

    if (it == table_.end())
    {
        return table_.back().p1;
    }

    auto const prev = it - 1;
    double const delta = prev->cosTheta - it->cosTheta;
    double const wIt = (prev->cosTheta - cosTheta) / delta;
    double const wPrev = (cosTheta - it->cosTheta) / delta;

    return it->p1 * wIt + prev->p1 * wPrev;
}


double MieDust::cosRandomTheta(double const v) const
{
    auto const it = std::lower_bound(accumulated_.begin(), accumulated_.end(), 1. - v);

    if (it == accumulated_.begin())
    {
        return 1.0;
    }

    if (it == accumulated_.end())
    {
        return -1.0;
    }

    auto const idx = std::distance(accumulated_.begin(), it);
    double const delta = accumulated_.at(idx) - accumulated_.at(idx - 1);
    double const wPrev = (accumulated_.at(idx) - 1. + v) / delta;
    double const wCurr = (1. - v - accumulated_.at(idx - 1)) / delta;

    return table_.at(idx).cosTheta * wCurr + table_.at(idx - 1).cosTheta * wPrev;
}
