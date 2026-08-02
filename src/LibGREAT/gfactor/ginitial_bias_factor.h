#ifndef GINITIAL_BIAS_FACTOR_H
#define GINITIAL_BIAS_FACTOR_H
/**
 * @file         ginitial_bias_factor.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        Initialization of accelerometer bias factor and gyroscope bias factor.
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
	*@brief InitialBiasFactor Class for initializing bias factors
	*/
	class InitialBiasFactor : public ceres::SizedCostFunction<6, 9>
	{
	public:

		/**
		 * @brief Initial bias prior factor constructor
		 * Sets accelerometer and gyroscope bias constraints with uncertainty weighting
		 */
		InitialBiasFactor(const Eigen::Vector3d &_Ba, const Eigen::Vector3d &_Bg)
		{
			init_Ba = _Ba;
			init_Bg = _Bg;
			sqrt_info = 1.0 / (0.001) * Eigen::Matrix<double, 6, 6>::Identity();
		}

		/**
		 * @brief Evaluate initial bias prior factor
		 *
		 * Computes residuals and Jacobians for IMU bias constraints.
		 * Handles accelerometer and gyroscope bias parameters with proper indexing.
		 */
		virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
		{
			Eigen::Vector3d Ba(parameters[0][3], parameters[0][4], parameters[0][5]);
			Eigen::Vector3d Bg(parameters[0][6], parameters[0][7], parameters[0][8]);

			Eigen::Map<Eigen::Matrix<double, 6, 1>> residual(residuals);
			residual.block<3, 1>(0, 0) = Ba - init_Ba;
			residual.block<3, 1>(3, 0) = Bg - init_Bg;
			residual = sqrt_info * residual;

			if (jacobians)
			{
				if (jacobians[0])
				{
					Eigen::Map<Eigen::Matrix<double, 6, 9, Eigen::RowMajor>> jacobian_bias(jacobians[0]);
					jacobian_bias.setZero();
					jacobian_bias.block<3, 3>(0, 3) = Eigen::Matrix3d::Identity();
					jacobian_bias.block<3, 3>(3, 6) = Eigen::Matrix3d::Identity();
					jacobian_bias = sqrt_info * jacobian_bias;
				}
			}
			return true;
		}

		Eigen::Vector3d init_Ba, init_Bg;				///< bias error of accelerometer,bias error of gyroscope
		Eigen::Matrix<double, 6, 6> sqrt_info;			///< info weight
	};

/**
*@brief InitialBiasFactor Class for initializing vel and bias factors
*/
	class InitialVelBiasFactor : public ceres::SizedCostFunction<9, 9>
	{
	public:

		/**
		 * @brief Initial velocity and bias prior factor constructor
		 * Sets velocity, accelerometer and gyroscope bias constraints with uncertainty weighting
		 */
		InitialVelBiasFactor(const Eigen::Vector3d& _Vel, const Eigen::Vector3d& _Ba, const Eigen::Vector3d& _Bg)
		{
			init_Vel = _Vel;
			init_Ba = _Ba;
			init_Bg = _Bg;
			sqrt_info = 1.0 / (0.001) * Eigen::Matrix<double, 9, 9>::Identity();
		}

		/**
		 * @brief Evaluate initial velocity and bias prior factor
		 *
		 * Computes residuals and Jacobians for velocity and IMU bias constraints.
		 * Handles velocity, accelerometer and gyroscope bias parameters with identity Jacobians.
		 */
		virtual bool Evaluate(double const* const* parameters, double* residuals, double** jacobians) const
		{
			Eigen::Vector3d Vel(parameters[0][0], parameters[0][1], parameters[0][2]);
			Eigen::Vector3d Ba(parameters[0][3], parameters[0][4], parameters[0][5]);
			Eigen::Vector3d Bg(parameters[0][6], parameters[0][7], parameters[0][8]);

			Eigen::Map<Eigen::Matrix<double, 9, 1>> residual(residuals);
			residual.block<3, 1>(0, 0) = Vel - init_Vel;
			residual.block<3, 1>(3, 0) = Ba - init_Ba;
			residual.block<3, 1>(6, 0) = Bg - init_Bg;
			residual = sqrt_info * residual;
			
			if (jacobians)
			{
				if (jacobians[0])
				{
					Eigen::Map<Eigen::Matrix<double, 9, 9, Eigen::RowMajor>> jacobian_bias(jacobians[0]);
					jacobian_bias.setZero();
					jacobian_bias.block<3, 3>(0, 0) = Eigen::Matrix3d::Identity();
					jacobian_bias.block<3, 3>(3, 3) = Eigen::Matrix3d::Identity();
					jacobian_bias.block<3, 3>(6, 6) = Eigen::Matrix3d::Identity();
					jacobian_bias = sqrt_info * jacobian_bias;
				}
			}
			return true;
		}

		Eigen::Vector3d init_Vel, init_Ba, init_Bg;
		Eigen::Matrix<double, 9, 9> sqrt_info;	 
	};

}
#endif