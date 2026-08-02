#ifndef GPOSE_LOCAL_PARAMETERIZATION_FACTOR_H
#define GPOSE_LOCAL_PARAMETERIZATION_FACTOR_H
/**
 * @file         gpose_local_parameterization_factor.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        Pose update for ceres.
 * @version      1.0
 * @date         2025-11-04
 *
 * @copyright Copyright (c) 2025, Wuhan University. All rights reserved.
 *
 */

#include"gfgo/gutility.h"
#include "gexport/ExportLibGREAT.h"
namespace gfgo
{
	/**
	*@brief  Class for pose parameterization
	*/
	class  LibGREAT_LIBRARY_EXPORT PoseLocalParameterization : public ceres::LocalParameterization
	{
		public:
		/**
		 * @brief Apply pose perturbation on manifold
		 * Updates position and quaternion using tangent space increments
		 * Maintains quaternion normalization for SO(3) manifold
		 */
		virtual bool Plus(const double *x, const double *delta, double *x_plus_delta) const;

		/**
		 * @brief Compute pose parameterization Jacobian
		 * Returns identity mapping for position and quaternion tangent space
		 */
		virtual bool ComputeJacobian(const double *x, double *jacobian) const;


		virtual int GlobalSize() const { return 7; };


		virtual int LocalSize() const { return 6; };
	};

}
#endif