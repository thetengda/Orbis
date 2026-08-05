#include "pseudorange_RAW_factor.h"
#include "raw_factor_common.h"

namespace gfgo
{
    PseudorangeRAWFactor::PseudorangeRAWFactor(const RAWEquMsg &message, const t_gallpar &params, t_gprecisebiasFGO *bias_model)
        : _message(message), _params(params), _bias_model(bias_model)
    {
    }

    bool PseudorangeRAWFactor::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
    {
        const double crd[3] = {parameters[0][0], parameters[0][1], parameters[0][2]};
        raw_factor_detail::RawLinearization linear;
        if (!raw_factor_detail::linearize(_message, _params, _bias_model, crd,
                                          parameters[1][0], parameters[2][0], parameters[3][0],
                                          0.0, 0.0, false, false, linear))
        {
            return false;
        }

        residuals[0] = -linear.sqrt_info * linear.residual;
        if (!jacobians)
            return true;

        if (jacobians[0])
        {
            for (int i = 0; i < 3; ++i)
                jacobians[0][i] = linear.sqrt_info * linear.coefficient[i];
        }
        if (jacobians[1])
            jacobians[1][0] = linear.sqrt_info * linear.coefficient[3];
        if (jacobians[2])
            jacobians[2][0] = linear.sqrt_info * linear.coefficient[4];
        if (jacobians[3])
            jacobians[3][0] = linear.sqrt_info * linear.coefficient[5];
        return true;
    }

    MultiPseudorangeRAWFactor::MultiPseudorangeRAWFactor(const RAWEquMsg &message, const t_gallpar &params, t_gprecisebiasFGO *bias_model)
        : _message(message), _params(params), _bias_model(bias_model)
    {
    }

    bool MultiPseudorangeRAWFactor::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
    {
        const double crd[3] = {parameters[0][0], parameters[0][1], parameters[0][2]};
        raw_factor_detail::RawLinearization linear;
        if (!raw_factor_detail::linearize(_message, _params, _bias_model, crd,
                                          parameters[1][0], parameters[2][0], parameters[3][0],
                                          parameters[4][0], 0.0, true, false, linear))
        {
            return false;
        }

        residuals[0] = -linear.sqrt_info * linear.residual;
        if (!jacobians)
            return true;

        if (jacobians[0])
        {
            for (int i = 0; i < 3; ++i)
                jacobians[0][i] = linear.sqrt_info * linear.coefficient[i];
        }
        if (jacobians[1])
            jacobians[1][0] = linear.sqrt_info * linear.coefficient[3];
        if (jacobians[2])
            jacobians[2][0] = linear.sqrt_info * linear.coefficient[4];
        if (jacobians[3])
            jacobians[3][0] = linear.sqrt_info * linear.coefficient[5];
        if (jacobians[4])
            jacobians[4][0] = linear.sqrt_info * linear.isb_coefficient;
        return true;
    }
}
