#include "SafierWind.hpp"
#include <cmath>

struct Model
{
    double dxi0;
    double a_1, a_2, a_3, a_4, a_5, a_6, a_7, a_8;
    double b_0, b_1, b_2, b_3, b_4, b_5, b_6, b_7, b_8, b_9;
};

namespace
{
    constexpr Model modelB{1.73,
                           -0.30, 0.54, 1.14, 0.26, 1.76, 0.92, 0.09, 0.27,
                           0.035, 7.7, 1.01, 0.31, 1.1, 0.40, -0.10, 0.0, 1.25, 0.60};

    Model const* getModel(char const key)
    {
        if (key == 'B')
        {
            return &modelB;
        }

        return nullptr;
    }

    double rho0(const Model* const model, double mOut, double mStar, double h0, double rRatio)
    {
        return 1.064e-15 * (mOut/1e-7) * std::pow(mStar/0.5, -0.5) / std::log(rRatio) /
            model->b_0/(1.0-h0*model->dxi0);
    }

    double xi(const Model* const model, double h)
    {
        double xi1 = (1+model->a_1*h + model->a_2*pow(h, model->a_3)) * exp(-model->a_4*h);
        double xi2 = model->a_5 * pow(h, model->a_6)  * exp(-4*model->a_7*pow(h, model->a_8));
        return xi1 + xi2;
    }

    double psi(const Model* const model, double h)
    {
        double psi1 = (model->b_0 + model->b_6*h + model->b_7*h*h) * exp(-model->b_8 * pow(h, model->b_9));
        double psi2 = model->b_1 * exp(-1.0/(model->b_2*pow(h, model->b_3)))  * exp(-model->b_4/(h*model->b_5));
        return psi1 + psi2;
    }

    double dXidChi(const Model* const model, double h)
    {
        double dXidChi1 = (model->a_1 + model->a_2*model->a_3*pow(h, model->a_3 - 1.0) - model->a_4 - model->a_1*model->a_4*h - model->a_2*model->a_4*pow(h, model->a_3))
                            *exp(-model->a_4*h);
        double dXidChi2 = (model->a_5 * model->a_6 * pow(h, model->a_6-1.0) -4*model->a_5*model->a_7*model->a_8*pow(h,model->a_6)*pow(h, model->a_8 -1.0) )
                            * exp(-4*model->a_7*pow(h, model->a_8));

        return dXidChi1 + dXidChi2;
    }

    double eta(const Model* const model, double chi, double h0)
    {
        double const h = fabs(chi) - h0;
        double const xi = ::xi(model, h);
        double const psi = ::psi(model, h);
        double const dxi = dXidChi(model, h);
        return model->b_0 * (1. - h0*model->dxi0) / (xi*psi*(xi-chi*dxi));
    }
}


SafierWind::SafierWind(
        char const model,
        double const mOut,
        double const mStar,
        double const h0,
        double const rMin,
        double const rMax)
    : model_{ getModel(model) }
    , rho0_{ rho0(model_, mOut, mStar, h0, rMax / rMin) }
    , h0_{ h0 }
{}


double SafierWind::density(double x, double y, double z) const
{
    double r = sqrt(x*x + y*y);
    double chi = z/r;

    if (chi < h0_ || rho0_ == 0.0)
    {
        return 0.0;
    }

    return rho0_ * std::pow(r, -1.5) * eta(model_, chi, h0_);
}