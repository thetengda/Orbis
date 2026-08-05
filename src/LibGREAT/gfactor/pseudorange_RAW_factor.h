#ifndef P_RAW_FACTOR_H
#define P_RAW_FACTOR_H

#include "gfgo/gutility.h"
#include "gfactor/raw_equation.h"
#include "gfgo/gprecisebiasFGO.h"

using namespace great;
using namespace gnut;

namespace gfgo
{
    class LibGREAT_LIBRARY_EXPORT PseudorangeRAWFactor : public ceres::SizedCostFunction<1, 3, 1, 1, 1>
    {
    public:
        PseudorangeRAWFactor(const RAWEquMsg &message, const t_gallpar &params, t_gprecisebiasFGO *bias_model);
        bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const override;

    private:
        RAWEquMsg _message;
        t_gallpar _params;
        t_gprecisebiasFGO *_bias_model = nullptr;
    };

    class LibGREAT_LIBRARY_EXPORT MultiPseudorangeRAWFactor : public ceres::SizedCostFunction<1, 3, 1, 1, 1, 1>
    {
    public:
        MultiPseudorangeRAWFactor(const RAWEquMsg &message, const t_gallpar &params, t_gprecisebiasFGO *bias_model);
        bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const override;

    private:
        RAWEquMsg _message;
        t_gallpar _params;
        t_gprecisebiasFGO *_bias_model = nullptr;
    };
}

#endif
