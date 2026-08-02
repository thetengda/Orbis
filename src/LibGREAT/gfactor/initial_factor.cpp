#include"initial_factor.h"

gfgo::InitialFactor::InitialFactor(const double & value, const double & sqrt_info) :
	_value(value), _sqrt_info(sqrt_info)
{}

bool gfgo::InitialFactor::Evaluate(double const * const * parameters, double * residuals, double ** jacobians) const
{
	double X = parameters[0][0];
	double o = _value;
	double c = X;
	double cmo = c - o;

	residuals[0] = _sqrt_info * cmo;

	if (jacobians)
	{
		if (jacobians[0])
		{
			Eigen::Map<Eigen::Matrix<double, 1, 1, Eigen::RowMajor>> jacobian_X(jacobians[0]);
			jacobian_X = _sqrt_info * Eigen::Matrix<double, 1, 1>::Identity();
		}
	}
	return true;
}