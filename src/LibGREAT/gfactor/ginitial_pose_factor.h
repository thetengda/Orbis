#ifndef GINITIAL_POSE_FACTOR_H
#define GINITIAL_POSE_FACTOR_H
/**
 * @file         ginitial_pose_factor.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        Initialization of pose info factor.
 * @version      1.0
 * @date         2025-11-04
 *
 * @copyright Copyright (c) 2025, Wuhan University. All rights reserved.
 *
 */

#include"gfgo/gutility.h"

namespace gfgo
{
	/**
	*@brief InitialPoseFactor Class for initializing pose factors
	*/
	class InitialPoseFactor : public ceres::SizedCostFunction<6, 7>
	{
	public:

		/**
		 * @brief Initial pose prior factor constructor
		 * Sets initial position and attitude constraints with information matrix
		 */
		InitialPoseFactor(const Eigen::Vector3d &_P, const Eigen::Quaterniond &_Q)
		{
			init_P = _P;
			init_Q = _Q;
			sqrt_info = 1000 * Eigen::Matrix<double, 6, 6>::Identity();
		}

		/**
		 * @brief Evaluate initial pose prior factor
		 *
		 * Computes residuals and Jacobians for position and attitude constraints.
		 * Handles quaternion manifold for attitude errors with proper tangent space mapping.
		 */
		virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
		{
			Eigen::Vector3d P(parameters[0][0], parameters[0][1], parameters[0][2]);
			Eigen::Quaterniond Q(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);

			Eigen::Map<Eigen::Matrix<double, 6, 1>> residual(residuals);
			residual.block<3, 1>(0, 0) = P - init_P;
			residual.block<3, 1>(3, 0) = 2 * (init_Q.inverse() * Q).vec();
			residual = sqrt_info * residual;
			
			if (jacobians)
			{
				if (jacobians[0])
				{
					Eigen::Map<Eigen::Matrix<double, 6, 7, Eigen::RowMajor>> jacobian_pose(jacobians[0]);
					jacobian_pose.setZero();
					jacobian_pose.block<3, 3>(0, 0) = Eigen::Matrix3d::Identity();
					jacobian_pose.block<3, 3>(3, 3) = t_gfgo_utility::Qleft(init_Q.inverse() * Q).bottomRightCorner<3, 3>();
					jacobian_pose = sqrt_info * jacobian_pose;
					//std::cout << jacobian_pose << endl;
				}

			}
			return true;

		}

		/**
		 * @brief Debug function for pose factor verification
		 *
		 * Tests residual and Jacobian computation with numerical differentiation.
		 * Validates analytical derivatives against finite differences for quality assurance.
		 */
		void check(double **parameters)
		{
			double *res = new double[6];
			double **jaco = new double *[1];
			jaco[0] = new double[6 * 7];
			Evaluate(parameters, res, jaco);
			puts("check begins");

			puts("my");

			std::cout << Eigen::Map<Eigen::Matrix<double, 6, 1>>(res).transpose() << std::endl
				<< std::endl;
			std::cout << Eigen::Map<Eigen::Matrix<double, 6, 7, Eigen::RowMajor>>(jaco[0]) << std::endl
				<< std::endl;

			Eigen::Vector3d P(parameters[0][0], parameters[0][1], parameters[0][2]);
			Eigen::Quaterniond Q(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);

			Eigen::Matrix<double, 6, 1> residual;
			residual.block<3, 1>(0, 0) = P - init_P;
			residual.block<3, 1>(3, 0) = 2 * (init_Q.inverse()* Q).vec();
			residual = sqrt_info * residual;

			puts("num");
			std::cout << residual.transpose() << std::endl;

			const double eps = 1e-6;
			Eigen::Matrix<double, 6, 6> num_jacobian;
			for (int k = 0; k < 6; k++)
			{
				Eigen::Vector3d P(parameters[0][0], parameters[0][1], parameters[0][2]);
				Eigen::Quaterniond Q(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);

				int a = k / 3, b = k % 3;
				Eigen::Vector3d delta = Eigen::Vector3d(b == 0, b == 1, b == 2) * eps;

				if (a == 0)
					P += delta;
				else if (a == 1)
					Q = Q * t_gfgo_utility::deltaQ(delta);


				Eigen::Matrix<double, 6, 1> tmp_residual;
				tmp_residual.block<3, 1>(0, 0) = P - init_P;
				tmp_residual.block<3, 1>(3, 0) = 2 * (init_Q.inverse()* Q).vec();
				tmp_residual = sqrt_info * tmp_residual;

				num_jacobian.col(k) = (tmp_residual - residual) / eps;
			}
			std::cout << num_jacobian << std::endl;

		}

		Eigen::Vector3d init_P;
		Eigen::Quaterniond init_Q;
		Eigen::Matrix<double, 6, 6> sqrt_info;
	};

	class InitialPosFactor : public ceres::SizedCostFunction<3, 3>
	{
	public:

		/**
		 * @brief Initial position prior factor constructor
		 * Sets position constraints with diagonal information matrix weights
		 */
		InitialPosFactor(const Eigen::Vector3d& _P)
		{
			init_P = _P;
			sqrt_info =  Eigen::Matrix<double, 3, 3>::Identity();
			for (int i = 0; i < 3; i++)
				sqrt_info(i, i) = 0.033;
		}

		/**
		 * @brief Evaluate initial position prior factor
		 *
		 * Computes residuals and Jacobians for position constraints.
		 * Simple linear factor with identity Jacobian for coordinate priors.
		 */
		virtual bool Evaluate(double const* const* parameters, double* residuals, double** jacobians) const
		{
			Eigen::Vector3d P(parameters[0][0], parameters[0][1], parameters[0][2]);

			Eigen::Map<Eigen::Matrix<double, 3, 1>> residual(residuals);
			residual.block<3, 1>(0, 0) = P - init_P;
			residual = sqrt_info * residual;

			if (jacobians)
			{
				if (jacobians[0])
				{
					Eigen::Map<Eigen::Matrix<double, 3, 3, Eigen::RowMajor>> jacobian_pose(jacobians[0]);
					jacobian_pose.setZero();
					jacobian_pose.block<3, 3>(0, 0) = Eigen::Matrix3d::Identity();
					jacobian_pose = sqrt_info * jacobian_pose;
				}

			}
			return true;

		}

		/**
		 * @brief Debug function for position factor verification
		 *
		 * Tests residual and Jacobian computation with numerical differentiation.
		 * Validates analytical derivatives against finite differences for position priors.
		 */
		void check(double** parameters)
		{
			double* res = new double[3];
			double** jaco = new double* [1];
			jaco[0] = new double[3 * 3];
			Evaluate(parameters, res, jaco);
			puts("check begins");

			puts("my");

			std::cout << Eigen::Map<Eigen::Matrix<double, 3, 1>>(res).transpose() << std::endl
				<< std::endl;
			std::cout << Eigen::Map<Eigen::Matrix<double, 3, 3, Eigen::RowMajor>>(jaco[0]) << std::endl
				<< std::endl;

			Eigen::Vector3d P(parameters[0][0], parameters[0][1], parameters[0][2]);

			Eigen::Matrix<double, 3, 1> residual;
			residual.block<3, 1>(0, 0) = P - init_P;
			residual = sqrt_info * residual;

			puts("num");
			std::cout << residual.transpose() << std::endl;

			const double eps = 1e-6;
			Eigen::Matrix<double, 3, 3> num_jacobian;
			for (int k = 0; k < 3; k++)
			{
				Eigen::Vector3d P(parameters[0][0], parameters[0][1], parameters[0][2]);

				int a = k / 3, b = k % 3;
				Eigen::Vector3d delta = Eigen::Vector3d(b == 0, b == 1, b == 2) * eps;

				if (a == 0)
					P += delta;

				Eigen::Matrix<double, 3, 1> tmp_residual;
				tmp_residual.block<3, 1>(0, 0) = P - init_P;
				tmp_residual = sqrt_info * tmp_residual;

				num_jacobian.col(k) = (tmp_residual - residual) / eps;
			}
			std::cout << num_jacobian << std::endl;

		}

		Eigen::Vector3d init_P;
		Eigen::Matrix<double, 3, 3> sqrt_info;
	};


	/**
	*@brief InitialPoseFactor Class for initializing pose factors
	*/
	class InitialGnssCRD : public ceres::SizedCostFunction<3, 3>
	{
	public:
		InitialGnssCRD(const Eigen::Vector3d & gnss_ard)
		{
			init_P = gnss_ard;
			sqrt_info.setIdentity();
			for (int i = 0; i < 3; i++)
				sqrt_info(i, i) = 1.0 / 10000.0;
		}
		InitialGnssCRD(const Eigen::Vector3d & gnss_ard, const double sqrt)
		{
			init_P = gnss_ard;
			sqrt_info.setIdentity();
			for (int i = 0; i < 3; i++)
				sqrt_info(i, i) = sqrt;
		}
		InitialGnssCRD(const Eigen::Vector3d & gnss_ard, const Eigen::Vector3d & sqrt)
		{
			init_P = gnss_ard;
			sqrt_info.setIdentity();
			for (int i = 0; i < 3; i++)
				sqrt_info(i, i) = sqrt(i);
		}
		virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
		{
			Eigen::Vector3d P(parameters[0][0], parameters[0][1], parameters[0][2]);


			Eigen::Map<Eigen::Matrix<double, 3, 1>> residual(residuals);
			residual.block<3, 1>(0, 0) = P - init_P;
			residual = sqrt_info * residual;

			if (jacobians)
			{
				if (jacobians[0])
				{
					Eigen::Map<Eigen::Matrix<double, 3, 3, Eigen::RowMajor>> jacobian_pos(jacobians[0]);
					jacobian_pos.setZero();
					jacobian_pos.block<3, 3>(0, 0) = Eigen::Matrix3d::Identity();
					jacobian_pos = sqrt_info * jacobian_pos;
				}

			}
			return true;

		}
		Eigen::Vector3d init_P;
		Eigen::Matrix3d sqrt_info;
	};
}
#endif