/**
 * @file         gpose_local_parameterization_factor.cpp
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        Pose update for ceres.
 * @version      1.0
 * @date         2025-11-04
 *
 * @copyright Copyright (c) 2025, Wuhan University. All rights reserved.
 *
 */

#include "gpose_local_parameterization.h"

namespace gfgo
{
	bool PoseLocalParameterization::Plus(const double *x, const double *delta, double *x_plus_delta) const
	{
		Eigen::Map<const Eigen::Vector3d> _p(x);
		Eigen::Map<const Eigen::Quaterniond> _q(x + 3);

		Eigen::Map<const Eigen::Vector3d> dp(delta);

		Eigen::Quaterniond dq = t_gfgo_utility::deltaQ(Eigen::Map<const Eigen::Vector3d>(delta + 3));

		Eigen::Map<Eigen::Vector3d> p(x_plus_delta);
		Eigen::Map<Eigen::Quaterniond> q(x_plus_delta + 3);

		p = _p + dp;
		q = (_q * dq).normalized();

		return true;
	}
	bool PoseLocalParameterization::ComputeJacobian(const double *x, double *jacobian) const
	{
		Eigen::Map<Eigen::Matrix<double, 7, 6, Eigen::RowMajor>> j(jacobian);
		j.topRows<6>().setIdentity();
		j.bottomRows<1>().setZero();

		return true;

		// Eigen::Map<Eigen::Matrix<double, 7, 6, Eigen::RowMajor>> J(jacobian);
		// J.setZero();
		//
		// // position: p_new = p + delta_p
		// J.block<3, 3>(0, 0).setIdentity();
		//
		// // Quaternion is stored in memory as [qx, qy, qz, qw]
		// const Eigen::Vector3d qv(x[3], x[4], x[5]);
		// const double qw = x[6];
		//
		// Eigen::Matrix3d S;
		// S <<     0.0, -qv.z(),  qv.y(),
		// 	  qv.z(),     0.0, -qv.x(),
		// 	 -qv.y(),  qv.x(),    0.0;
		//
		// // q_new = q * deltaQ(delta_theta)
		// J.block<3, 3>(3, 3) = 0.5 * (qw * Eigen::Matrix3d::Identity() + S);
		// J.block<1, 3>(6, 3) = -0.5 * qv.transpose();

		return true;
	}


	

}