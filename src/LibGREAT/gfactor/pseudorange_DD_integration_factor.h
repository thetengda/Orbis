#ifndef P_DD_ING_FACTOR
#define P_DD_ING_FACTOR
/**
 * @file         pseudorange_DD_integration_factor.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        construct  pseudo-range factor for GNSS/INS integration
 * @version      1.0
 * @date         2025-11-04
 *
 * @copyright Copyright (c) 2025, Wuhan University. All rights reserved.
 *
 */


#include"gfgo/gutility.h"
#include "gexport/ExportLibGREAT.h"
#include "gfgo/gprecisebiasFGO.h"

using namespace great;

namespace gfgo
{
	class LibGREAT_LIBRARY_EXPORT PseudorangeDDINGFactor : public ceres::SizedCostFunction<1, 7>
	{
	public:

		/**
		 * @brief Pseudorange DD with INS/GNSS lever arm factor constructor
		 * Initializes DD factor including IMU-GNSS antenna lever arm offset
		 */
		PseudorangeDDINGFactor(const t_gtime &cur_time, const pair<string, string> &base_rover_site, const t_gallpar  &params, const vector<pair<t_gsatdata, t_gsatdata>> &DD_sat_data, t_gprecisebiasFGO *bias_model, const pair<FREQ_SEQ, GOBSBAND> &freq_band, const Eigen::Vector3d &lever_arm);
		
		/**
		 * @brief Update rover coordinates in parameter set for INS/GNSS
		 * Replaces rover position parameters considering lever arm transformation
		 */
		void updatePara(t_gallpar &params_tmp, const Eigen::Vector3d &Pi, const Eigen::Vector3d &Vi = Eigen::Vector3d::Identity()) const;
		
		/**
		 * @brief Transform design matrices to Eigen format for INS/GNSS pseudorange
		 * Converts sparse design matrix, weight matrix and residuals to dense Eigen matrices
		 * Handles coordinate parameters considering lever arm effects
		 */
		void trans2Eigen(const vector<vector<pair<int, double>>> &B, const vector<double> &P, const vector<double> &l, Eigen::Matrix<double, 2, 3> &B_new, Eigen::Matrix<double, 2, 2> &P_new, Eigen::Matrix<double, 2, 1> &l_new) const;
		
		/**
		 * @brief Evaluate INS/GNSS pseudorange DD factor with lever arm compensation
		 *
		 * Computes residuals and Jacobians considering IMU-GNSS antenna offset.
		 * Handles pose parameters and lever arm transformation for pseudorange observations.
		 * Transforms INS pose to GNSS antenna frame and constructs DD observations with proper covariance scaling for tight coupling integration.
		 */
		virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
		
		/**
		 * @brief Debug function for factor verification
		 * Tests residual and Jacobian computation with current parameters
		 * Outputs numerical results for validation and debugging purposes
		 */
		void check(double **parameters);
	protected:
		t_gtime _cur_time;
		pair<string, string> _base_rover_site;
		t_gallpar _params;
		//shared_ptr<t_gbiasmodel> _gprecise_bias_model = nullptr;
		t_gprecisebiasFGO *_gprecise_bias_model = nullptr;
		pair<FREQ_SEQ, GOBSBAND> _freq_band;
		vector<pair<t_gsatdata, t_gsatdata>> _DD_sat_data;
		Eigen::Vector3d _lever_arm;
	};
}
#endif