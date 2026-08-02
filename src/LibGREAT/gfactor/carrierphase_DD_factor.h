#ifndef L_DD_FACTOR
#define L_DD_FACTOR
/**
 * @file         carrierphase_DD_factor.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        Class for factor graph-based GNSS/SINS integrated solution
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
	class LibGREAT_LIBRARY_EXPORT CarrierphaseDDFactor : public ceres::SizedCostFunction<1, 3, 1, 1>
	{
	public:
		
		/**
		 * @brief Carrier phase DD factor constructor
		 * Initializes with time, site pairs, parameters, satellite data and frequency band
		 */
		CarrierphaseDDFactor(const t_gtime &cur_time, const pair<string, string> &base_rover_site, const t_gallpar  &params, const vector<pair<t_gsatdata, t_gsatdata>> &DD_sat_data, t_gprecisebiasFGO *bias_model, const pair<FREQ_SEQ, GOBSBAND> &freq_band);
		
		/**
		 * @brief Update parameters for carrier phase DD factor
		 * Updates rover coordinates, reference and non-reference satellite ambiguities in parameter set for linearization point adjustment
		 */
		void updatePara(t_gallpar & params_tmp, const double &ref_sd_amb, const double &nonref_sd_amb, const Eigen::Vector3d &Pi, const Eigen::Vector3d &Vi = Eigen::Vector3d::Identity()) const;
		
		/**
		 * @brief Transform design matrices to Eigen format
		 * Converts sparse design matrix, weight matrix and residuals to dense Eigen matrices for efficient linear algebra operations in optimization
		 */
		void trans2Eigen(const vector<vector<pair<int, double>>> &B, const vector<double> &P, const vector<double> &l, Eigen::Matrix<double, 2, 5> &B_new, Eigen::Matrix<double, 2, 2> &P_new, Eigen::Matrix<double, 2, 1> &l_new) const;		
		
		/**
		 * @brief Evaluate carrier phase double-difference factor for optimization
		 *
		 * Computes residuals and Jacobians for carrier phase DD observations in factor graph.
		 * Handles coordinate updates, ambiguity parameters, and measurement weighting.
		 * Constructs DD equations from single-difference observations and transforms to optimization-friendly format with proper covariance scaling.
		 */
		virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;

	/*	void ResetJacobianCache() const { _jac_cached = false; }
		void EnableFreezeJacobian(bool on) const { _freeze_jac = on; if (!on) _jac_cached = false; }*/
		
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

	//private:
	//	mutable bool _freeze_jac = true;
	//	mutable bool _jac_cached = false;

	//	// 缓存 Jacobian（whitened 后，与你最终写入 jacobians[] 的数值一致）
	//	mutable double _J_p[3] = { 0.0, 0.0, 0.0 };  // 对 parameters[0] (1x3)
	//	mutable double _J_a1 = 0.0;              // 对 parameters[1] (1x1)
	//	mutable double _J_a2 = 0.0;              // 对 parameters[2] (1x1)

	//	// 可选：缓存线性化点，用于“变化大则重算”保险（强烈建议加）
	//	mutable Eigen::Vector3d _lin_P = Eigen::Vector3d::Zero();
	//	mutable double _relinearize_thresh_m = 0.02; // 2cm，可自己调

	};
}
#endif