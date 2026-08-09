#ifndef P_RAW_FACTOR_H
#define P_RAW_FACTOR_H

#include "gfgo/gutility.h"
#include "gfactor/raw_factor_common.h"

using namespace great;
using namespace gnut;

namespace gfgo
{
    class LibGREAT_LIBRARY_EXPORT PseudorangeRAWFactor : public ceres::SizedCostFunction<1, 3, 1, 1, 1>, public raw_factor_detail::RawPreparableFactor
    {
    public:
        PseudorangeRAWFactor(const RAWEquMsg &message, const t_gallpar &params, t_gprecisebiasFGO *bias_model);
        bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const override;
        bool prepare(const std::vector<double *> &blocks) override;

    private:
        RAWEquMsg _message;
        mutable t_gallpar _params;
        t_gprecisebiasFGO *_bias_model = nullptr;
        mutable raw_factor_detail::RawEvaluationCache _cache;
    };

    class LibGREAT_LIBRARY_EXPORT MultiPseudorangeRAWFactor : public ceres::SizedCostFunction<1, 3, 1, 1, 1, 1>, public raw_factor_detail::RawPreparableFactor
    {
    public:
        MultiPseudorangeRAWFactor(const RAWEquMsg &message, const t_gallpar &params, t_gprecisebiasFGO *bias_model);
        bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const override;
        bool prepare(const std::vector<double *> &blocks) override;

    private:
        RAWEquMsg _message;
        mutable t_gallpar _params;
        t_gprecisebiasFGO *_bias_model = nullptr;
        mutable raw_factor_detail::RawEvaluationCache _cache;
    };

    class LibGREAT_LIBRARY_EXPORT PseudorangeRAWIFBFactor : public ceres::SizedCostFunction<1, 3, 1, 1, 1, 1>, public raw_factor_detail::RawPreparableFactor
    {
    public:
        PseudorangeRAWIFBFactor(const RAWEquMsg &message, const t_gallpar &params, t_gprecisebiasFGO *bias_model);
        bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const override;
        bool prepare(const std::vector<double *> &blocks) override;

    private:
        RAWEquMsg _message;
        mutable t_gallpar _params;
        t_gprecisebiasFGO *_bias_model = nullptr;
        mutable raw_factor_detail::RawEvaluationCache _cache;
    };

    class LibGREAT_LIBRARY_EXPORT MultiPseudorangeRAWIFBFactor : public ceres::SizedCostFunction<1, 3, 1, 1, 1, 1, 1>, public raw_factor_detail::RawPreparableFactor
    {
    public:
        MultiPseudorangeRAWIFBFactor(const RAWEquMsg &message, const t_gallpar &params, t_gprecisebiasFGO *bias_model);
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
