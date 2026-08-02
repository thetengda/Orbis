#include"random_walk_factor.h"

gfgo::RandomWalkFactor::RandomWalkFactor(const double & sqrt_info) :
	_sqrt_info(sqrt_info)
{}

bool gfgo::RandomWalkFactor::Evaluate(double const * const * parameters, double * residuals, double ** jacobians) const
{
	double X0 = parameters[0][0];
	double X1 = parameters[1][0];
	double o = 0.0;
	double c = -X0 + X1;
	double cmo = c - o;

	residuals[0] = _sqrt_info * cmo;

	if (jacobians)
	{
		if (jacobians[0])
		{
			Eigen::Map<Eigen::Matrix<double, 1, 1, Eigen::RowMajor>> jacobian_X0(jacobians[0]);
			jacobian_X0 = -_sqrt_info * Eigen::Matrix<double, 1, 1>::Identity();
		}
		if (jacobians[1])
		{
			Eigen::Map<Eigen::Matrix<double, 1, 1, Eigen::RowMajor>> jacobian_X1(jacobians[1]);
			jacobian_X1 = _sqrt_info * Eigen::Matrix<double, 1, 1>::Identity();
		}
	}
	return true;
}