#ifndef L_RAW_FACTOR_H
#define L_RAW_FACTOR_H

#include "gfgo/gutility.h"
#include "gfactor/raw_factor_common.h"

using namespace great;
using namespace gnut;

namespace gfgo
{
    class LibGREAT_LIBRARY_EXPORT CarrierphaseRAWFactor : public ceres::SizedCostFunction<1, 3, 1, 1, 1, 1>, public raw_factor_detail::RawPreparableFactor
    {
    public:
        CarrierphaseRAWFactor(const RAWEquMsg &message, const t_gallpar &params, t_gprecisebiasFGO *bias_model);
        bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const override;
        bool prepare(const std::vector<double *> &blocks) override;

    private:
        RAWEquMsg _message;
        mutable t_gallpar _params;
        t_gprecisebiasFGO *_bias_model = nullptr;
        mutable raw_factor_detail::RawEvaluationCache _cache;
    };

    class LibGREAT_LIBRARY_EXPORT MultiCarrierphaseRAWFactor : public ceres::SizedCostFunction<1, 3, 1, 1, 1, 1, 1>, public raw_factor_detail::RawPreparableFactor
    {
    public:
        MultiCarrierphaseRAWFactor(const RAWEquMsg &message, const t_gallpar &params, t_gprecisebiasFGO *bias_model);
        bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const override;
        bool prepare(const std::vector<double *> &blocks) override;

    private:
        RAWEquMsg _message;
        mutable t_gallpar _params;
        t_gprecisebiasFGO *_bias_model = nullptr;
        mutable raw_factor_detail::RawEvaluationCache _cache;
    };
}

#endif
