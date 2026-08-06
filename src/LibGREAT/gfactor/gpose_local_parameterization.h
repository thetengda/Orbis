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
	*@brief  Class for pose parameterization (Ceres 2.x Manifold, quaternion stored xyzw)
	*/
	class  LibGREAT_LIBRARY_EXPORT PoseLocalParameterization : public ceres::Manifold
	{
		public:
		/**
		 * @brief Apply pose perturbation on manifold
		 * Updates position and quaternion using tangent space increments
		 * Maintains quaternion normalization for SO(3) manifold
		 */
		virtual bool Plus(const double *x, const double *delta, double *x_plus_delta) const override;

		/**
		 * @brief Compute pose parameterization Jacobian (tangent -> ambient)
		 * Returns identity mapping for position and quaternion tangent space
		 */
		virtual bool PlusJacobian(const double *x, double *jacobian) const override;

		/**
		 * @brief Inverse of Plus: ambient difference -> tangent vector
		 */
		virtual bool Minus(const double *y, const double *x, double *y_minus_x) const override;

		/**
		 * @brief Jacobian of Minus w.r.t. y evaluated at y = x (ambient -> tangent)
		 */
		virtual bool MinusJacobian(const double *x, double *jacobian) const override;


		virtual int AmbientSize() const override { return 7; };


		virtual int TangentSize() const override { return 6; };
	};

}
#endif