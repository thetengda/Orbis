#ifndef L_DD_ING_FACTOR
#define L_DD_ING_FACTOR
/**
 * @file         carrierphase_DD_integration_factor.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        construct  carrier-phase factor for GNSS/INS integration
 * @version      1.0
 * @date         2025-11-04
 *
 * @copyright Copyright (c) 2025, Wuhan University. All rights reserved.
 *
 */

#include"gfgo/gutility.h"
#include "gfgo/gprecisebiasFGO.h"

using namespace great;
namespace gfgo
{
	class LibGREAT_LIBRARY_EXPORT CarrierphaseDDINGFactor : public ceres::SizedCostFunction<1, 7, 1, 1>
	{
	public:

		/**
		 * @brief Carrier phase DD with INS/GNSS lever arm factor constructor
		 * Initializes DD factor including IMU-GNSS antenna lever arm offset
		 */
		CarrierphaseDDINGFactor(const t_gtime &cur_time, const pair<string, string> &base_rover_site, const t_gallpar  &params, const vector<pair<t_gsatdata, t_gsatdata>> &DD_sat_data, t_gprecisebiasFGO *bias_model, const pair<FREQ_SEQ, GOBSBAND> &freq_band, const Eigen::Vector3d &lever_arm);
		
		/**
		 * @brief Update parameters for INS/GNSS carrier phase DD factor
		 * Updates rover coordinates and satellite ambiguities including lever arm effects
		 * Maintains parameter consistency for INS and GNSS integration
		 */
		void updatePara(t_gallpar & params_tmp, const double &ref_sd_amb, const double &nonref_sd_amb, const Eigen::Vector3d &Pi, const Eigen::Vector3d &Vi = Eigen::Vector3d::Identity()) const;
		
		/**
		 * @brief Transform design matrices to Eigen format for INS/GNSS
		 * Converts sparse design matrix, weight matrix and residuals to dense Eigen matrices
		 */
		void trans2Eigen(const vector<vector<pair<int, double>>> &B, const vector<double> &P, const vector<double> &l, Eigen::Matrix<double, 2, 5> &B_new, Eigen::Matrix<double, 2, 2> &P_new, Eigen::Matrix<double, 2, 1> &l_new) const;
		
		/**
		 * @brief Evaluate INS/GNSS carrier phase DD factor with lever arm compensation
		 *
		 * Computes residuals and Jacobians considering IMU-GNSS antenna offset.
		 * Handles pose parameters, lever arm transformation, and ambiguity states.
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

