#ifndef GINITIAL_AMB_FACTOR_H
#define GINITIAL_AMB_FACTOR_H
/**
 * @file         ginitial_gnss_amb_factor.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        Initialization of gnss ambiguity info factor.
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
	class InitialGnssAMB : public ceres::SizedCostFunction<1, 1>
	{
	public:

		/**
		 * @brief GNSS ambiguity prior factor constructor
		 * Sets initial ambiguity constraint with uncertainty weighting
		 */
		InitialGnssAMB(const double & gnss_ambiguity)
		{
			initial_amb = gnss_ambiguity;
			sqrt_info = 0.033;
		}
		// to control the strength of the ambiguity constraint
		InitialGnssAMB(const double & gnss_ambiguity, const double & set_sqrt_info)
		{
			initial_amb = gnss_ambiguity;
			sqrt_info = set_sqrt_info;
		}
		/**
		 * @brief Evaluate GNSS ambiguity prior factor
		 *
		 * Computes residual and Jacobian for ambiguity parameter constraint.
		 * Simple linear factor for float ambiguity initialization.
		 */
		virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
		{
			double c_amb = parameters[0][0];
			
			double res = c_amb - initial_amb;
			
			residuals[0] = sqrt_info * res;
			if (jacobians)
			{
				if (jacobians[0])
				{
					Eigen::Map<Eigen::Matrix<double, 1, 1, Eigen::RowMajor>> jacobian_ambiguity1(jacobians[0]);
					jacobian_ambiguity1 = sqrt_info* Eigen::Matrix<double, 1, 1>::Identity();
				}

			}
			return true;

		}				
		double initial_amb;
		double  sqrt_info;
	};

}
#endif