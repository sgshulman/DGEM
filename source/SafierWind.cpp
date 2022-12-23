#include "SafierWind.hpp"
#include "DebugUtils.hpp"
#include "IDiskHump.hpp"
#include <cmath>
#include <string>

struct SafierWindModel
{
    double dxi0;
    double a_1, a_2, a_3, a_4, a_5, a_6, a_7, a_8;
    double b_0, b_1, b_2, b_3, b_4, b_5, b_6, b_7, b_8, b_9;
};

namespace
{
    constexpr SafierWindModel modelB{1.73,
        -0.30, 0.54, 1.14, 0.26, 1.76, 0.92, 0.09, 0.27,
        0.035, 7.7, 1.01, 0.31, 1.1, 0.40, -0.10, 0.0, 1.25, 0.60};

    constexpr SafierWindModel modelC{1.73,
        0.21, 0.30, 1.23, 0.21, 1.27, 0.92, 0.04, 0.28,
        0.035, 12.16, 0.65, 0.33, 1.0, 0.40, 0.50, 0.0, 1.25, 0.60};

    constexpr SafierWindModel modelD {1.73,
        0.49, 0.17, 1.27, 0.11, 0.89, 0.97, 0.02, 0.31,
        0.035, 20.08, 0.42, 0.34, 1.0, 0.40, 0.90, 0.0, 1.00, 0.60};

    constexpr SafierWindModel modelI{1.73,
        0.18, 0.14, 1.86, 0.26, 1.30, 0.92, 0.02, 0.24,
        0.035, 33.52, 0.55, 0.26, 2.5, 0.50, 1.00, 0.0, 0.75, 1.00};

    constexpr SafierWindModel modelE{3.73,
        -0.51, 0.62, 1.02, 0.03, 2.46, 0.91, 0.01, 0.32,
        0.01, 6.95, 0.53, 0.30, 1.0, 0.40, 0.90, 0.0, 1.00, 0.60};

    constexpr SafierWindModel modelF{3.73,
        -0.22, 0.50, 0.98, 0.09, 2.78, 0.92, 0.01, 0.27,
        0.01, 11.62, 0.67, 0.20, 1.8, 0.40, 1.00, 0.0, 0.85, 0.55};

    constexpr SafierWindModel modelG {3.73,
        0.07, 0.05, 0.42, 0.05, 2.95, 0.94, 0.008, 0.26,
        0.01, 18.70, 0.78, 0.14, 3.3, 0.55, 1.00, 0.0, 0.30, 1.10};

    SafierWindModel const* getModel(char const key)
    {
        if (key == 'B')
        {
            return &modelB;
        } else if (key == 'C') {
            return &modelC;
        } else if (key == 'D') {
            return &modelD;
        } else if (key == 'I') {
            return &modelI;
        } else if (key == 'E') {
            return &modelE;
        } else if (key == 'F') {
            return &modelF;
        } else if (key == 'G') {
            return &modelG;
        }

        return nullptr;
    }

    double rho0(const SafierWindModel* const model, double mOut, double mStar, double h0, double rRatio)
    {
        return 1.064e-15 * (mOut/1e-7) * std::pow(mStar/0.5, -0.5) / std::log(rRatio) /
            model->b_0/(1.0-h0*model->dxi0);
    }

    double xi(const SafierWindModel* const model, double h)
    {
        double xi1 = (1+model->a_1*h + model->a_2*pow(h, model->a_3)) * exp(-model->a_4*h);
        double xi2 = model->a_5 * pow(h, model->a_6)  * exp(-4*model->a_7*pow(h, model->a_8));
        return xi1 + xi2;
    }

    double psi(const SafierWindModel* const model, double h)
    {
        double psi1 = (model->b_0 + model->b_6*h + model->b_7*h*h) * exp(-model->b_8 * pow(h, model->b_9));
        double psi2 = model->b_1 * exp(-1.0/(model->b_2*pow(h, model->b_3)))  * exp(-model->b_4/(h*model->b_5));
        return psi1 + psi2;
    }

    double dXidChi(const SafierWindModel* const model, double h)
    {
        double dXidChi1 = (model->a_1 + model->a_2*model->a_3*pow(h, model->a_3 - 1.0) - model->a_4 - model->a_1*model->a_4*h - model->a_2*model->a_4*pow(h, model->a_3))
                            *exp(-model->a_4*h);
        double dXidChi2 = (model->a_5 * model->a_6 * pow(h, model->a_6-1.0) -4*model->a_5*model->a_7*model->a_8*pow(h,model->a_6)*pow(h, model->a_8 -1.0) )
                            * exp(-4*model->a_7*pow(h, model->a_8));

        return dXidChi1 + dXidChi2;
    }

    double eta(const SafierWindModel* const model, double chi, double h0)
    {
        double const h = fabs(chi) - h0;
        double const xi = ::xi(model, h);
        double const psi = ::psi(model, h);
        double const dxi = dXidChi(model, h);
        return model->b_0 * (1. - h0*model->dxi0) / (xi*psi*(xi-chi*dxi));
    }
} // namespace


SafierWind::SafierWind(
        char const model,
        double const mOut,
        double const mStar,
        double const h0,
        double const rMin,
        double const rMax,
        IDiskHumpCPtr hump)
    : model_{ getModel(model) }
    , h0_{ h0 }
    , hump_{ std::move(hump) }
{
    DATA_ASSERT(
        std::string("BCDIEFG").find(model) != std::string::npos,
        "Safier Wind model must be B, C, D, I, E, F, or G");

    DATA_ASSERT(mOut > 0., "mOut (the mass outflow rate of the wind) must be positive.");
    DATA_ASSERT(mStar > 0., "mStar (the stellar mass) must be positive.");
    DATA_ASSERT(h0 > 0., "h0 (density of the sphere envelope) must be positive.");
    DATA_ASSERT(rMin > 0., "rMin (inner radius of the wind formation region) must be positive.");
    DATA_ASSERT(rMax > 0., "rMax (outer radius of the wind formation region) must be positive.");
    DATA_ASSERT(rMax > rMin, "rMax (outer radius) must be greater than rMin (inner radius of the wind formation region).");

    if (model_ != nullptr)
    {
        rho0_ = rho0(model_, mOut, mStar, h0, rMax / rMin);
    }
}


double SafierWind::density(Vector3d const& position) const
{
    double const r = std::sqrt(position.x()*position.x() + position.y()*position.y());
    double const chi = std::abs(position.z()/r);

    if (chi < h0_ || rho0_ == 0.0)
    {
        return 0.0;
    }

    double const rho = rho0_ * std::pow(r, -1.5) * eta(model_, chi, h0_);

    if (hump_)
    {
        return hump_->hump(rho, position);
    }

    return rho;
}
