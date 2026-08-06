#ifndef FIXED_AMBIGUITY_FACTOR_H
#define FIXED_AMBIGUITY_FACTOR_H

#include <ceres/ceres.h>

namespace gfgo
{
    /** Absolute ambiguity equation ca * a + cb * b = target. */
    class FixedAmbiguityFactor : public ceres::SizedCostFunction<1, 1, 1>
    {
    public:
        FixedAmbiguityFactor(double coefficient_a, double coefficient_b,
                             double target, double sqrt_information)
            : _coefficient_a(coefficient_a),
              _coefficient_b(coefficient_b),
              _target(target),
              _sqrt_information(sqrt_information)
        {
        }

        bool Evaluate(double const *const *parameters, double *residuals,
                      double **jacobians) const override
        {
            if (!parameters || !parameters[0] || !parameters[1] || !residuals)
                return false;

            residuals[0] = _sqrt_information *
                           (_coefficient_a * parameters[0][0] +
                            _coefficient_b * parameters[1][0] - _target);
            if (jacobians)
            {
                if (jacobians[0])
                    jacobians[0][0] = _sqrt_information * _coefficient_a;
                if (jacobians[1])
                    jacobians[1][0] = _sqrt_information * _coefficient_b;
            }
            return true;
        }

    private:
        double _coefficient_a;
        double _coefficient_b;
        double _target;
        double _sqrt_information;
    };
}

#endif
