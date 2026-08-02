#ifndef GALIGNMENT_GNSS_PRIOR_FACTOR_H
#define GALIGNMENT_GNSS_PRIOR_FACTOR_H

#include "gfgo/gutility.h"

namespace gfgo
{
	class AlignmentGNSSPriorFactor : public ceres::CostFunction
	{
	public:
		AlignmentGNSSPriorFactor(
			const Eigen::VectorXd& mean,
			const Eigen::MatrixXd& sqrt_info,
			int scalar_block_count)
			: _mean(mean), _sqrt_info(sqrt_info)
		{
			set_num_residuals(static_cast<int>(_mean.size()));
			mutable_parameter_block_sizes()->push_back(7);
			for (int i = 0; i < scalar_block_count; ++i)
				mutable_parameter_block_sizes()->push_back(1);
		}

		bool Evaluate(
			double const* const* parameters,
			double* residuals,
			double** jacobians) const override
		{
			const int state_size = static_cast<int>(_mean.size());
			Eigen::VectorXd delta(state_size);
			delta.head<3>() = Eigen::Vector3d(
				parameters[0][0], parameters[0][1], parameters[0][2]) -
				_mean.head<3>();

			for (int i = 3; i < state_size; ++i)
				delta(i) = parameters[i - 2][0] - _mean(i);

			Eigen::Map<Eigen::VectorXd> residual_map(residuals, state_size);
			residual_map = _sqrt_info * delta;

			if (jacobians != nullptr)
			{
				if (jacobians[0] != nullptr)
				{
					Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 7, Eigen::RowMajor>>
						pose_jacobian(jacobians[0], state_size, 7);
					pose_jacobian.setZero();
					pose_jacobian.leftCols<3>() = _sqrt_info.leftCols<3>();
				}

				for (int i = 3; i < state_size; ++i)
				{
					if (jacobians[i - 2] == nullptr)
						continue;
					Eigen::Map<Eigen::VectorXd> scalar_jacobian(
						jacobians[i - 2], state_size);
					scalar_jacobian = _sqrt_info.col(i);
				}
			}
			return true;
		}

	private:
		Eigen::VectorXd _mean;
		Eigen::MatrixXd _sqrt_info;
	};
}

#endif
