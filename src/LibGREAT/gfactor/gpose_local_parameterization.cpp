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
	bool PoseLocalParameterization::PlusJacobian(const double *x, double *jacobian) const
	{
		Eigen::Map<Eigen::Matrix<double, 7, 6, Eigen::RowMajor>> j(jacobian);
		j.topRows<6>().setIdentity();
		j.bottomRows<1>().setZero();

		return true;

		// A fuller Jacobian (still 7x6, tangent->ambient); kept for reference.
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

	bool PoseLocalParameterization::Minus(const double *y, const double *x, double *y_minus_x) const
	{
		Eigen::Map<const Eigen::Vector3d> yp(y), xp(x);
		Eigen::Map<const Eigen::Quaterniond> yq(y + 3), xq(x + 3);

		Eigen::Map<Eigen::Vector3d> dp(y_minus_x);
		Eigen::Map<Eigen::Vector3d> dtheta(y_minus_x + 3);

		// Plus(x, delta) = y  =>  delta_p = y_p - x_p ;  dtheta = 2*Im(xq^-1 * yq)
		dp = yp - xp;
		Eigen::Quaterniond dq = (xq.conjugate() * yq).normalized();
		dtheta = 2.0 * dq.vec();

		return true;
	}

	bool PoseLocalParameterization::MinusJacobian(const double *x, double *jacobian) const
	{
		// d(Minus)/dy | (y=x), size 6x7. Position part: [I3 0]; rotation part from
		// d(2*Im(xq^-1*yq))/dy. Quaternion stored xyzw at x+3.
		Eigen::Map<Eigen::Matrix<double, 6, 7, Eigen::RowMajor>> j(jacobian);
		j.setZero();
		j.topLeftCorner<3, 3>().setIdentity();

		const double xw = x[6];
		const Eigen::Vector3d xv(x[3], x[4], x[5]);
		Eigen::Matrix3d S;   // skew of xv
		S <<  0.0, -xv.z(),  xv.y(),
		      xv.z(),   0.0, -xv.x(),
		     -xv.y(),  xv.x(),   0.0;
		j.block<3, 3>(3, 3) = 2.0 * (xw * Eigen::Matrix3d::Identity() - S);
		j.block<3, 1>(3, 6) = -2.0 * xv;

		return true;
	}


	

}