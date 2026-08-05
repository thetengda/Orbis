#ifndef GNSS_INFO_FACTOR
#define GNSS_INFO_FACTOR
/**
 * @file         gnss_info_factor.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        construct Double-Different carrier-phase factor 
 * @version      1.0
 * @date         2025-11-04
 *
 * @copyright Copyright (c) 2025, Wuhan University. All rights reserved.
 *
 */

#include"gfgo/gutility.h"
#include "gexport/ExportLibGREAT.h"
#include <cstdint>

namespace gfgo
{
	using ParameterBlockKey = std::uintptr_t;

	/**
		*@brief class for construting residual blocks (for single observation equation)
		*/
	struct LibGREAT_LIBRARY_EXPORT GNSSResidualBlockInfo
	{
		GNSSResidualBlockInfo(ceres::CostFunction *_cost_function, ceres::LossFunction *_loss_function, const std::vector<double *> &_parameter_blocks)
			: cost_function(_cost_function), loss_function(_loss_function), parameter_blocks(_parameter_blocks) {}

		GNSSResidualBlockInfo(ceres::CostFunction *_cost_function, ceres::LossFunction *_loss_function, const std::vector<double *> &_parameter_blocks, const std::vector<int> &_drop_set)
			: cost_function(_cost_function), loss_function(_loss_function), parameter_blocks(_parameter_blocks), drop_set(_drop_set) {}
		
		/**
		 * @brief Evaluate GNSS residual block with loss function
		 *
		 * Computes residuals and Jacobians for GNSS factors, applies robust loss functions, and scales residuals/Jacobians for outlier rejection in optimization.
		 */
		void Evaluate();
		ceres::CostFunction *cost_function;
		ceres::LossFunction *loss_function;
		std::vector<double *> parameter_blocks;
		std::vector<int> drop_set;
		double **raw_jacobians = nullptr;
		std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> jacobians;  //sqrt(P)*A  
		Eigen::VectorXd residuals; //sqrt(P)*l	
	};	

	class LibGREAT_LIBRARY_EXPORT GNSSInfo
	{

	public:
		GNSSInfo() { valid = false; };
		~GNSSInfo();

		/**
		 * @brief Add GNSS residual block to information structure
		 *
		 * Stores GNSS factors and parameter indexing for covariance computation.
		 */
		void addResidualBlockInfo(GNSSResidualBlockInfo *residual_block_info, const map<ParameterBlockKey, vector<int>> &para_ids);
		
		/**
		 * @brief Add GNSS residual block for marginalization
		 *
		 * Stores GNSS factors and identifies parameter blocks to be marginalized.
		 * Tracks parameter sizes and drop sets for Schur complement computation.
		 */
		void addResidualBlockInfo(GNSSResidualBlockInfo *residual_block_info);
		
		/**
		 * @brief Construct GNSS observation equations from covariance matrix
		 *
		 * Builds linearized observation equations for posterior analysis and outlier detection.
		 * Computes design matrix, residuals, covariance matrices, and statistical metrics including normalized residuals and unit variance for quality assessment.
		 */
		void constructEqu(Eigen::MatrixXd variance);
		
		/**
		 * @brief Construct GNSS equations with covariance type selection
		 *
		 * Builds linearized observation equations with configurable covariance handling.
		 * Supports both prior-included and prior-excluded covariance computation for different stages of GNSS processing and quality control.
		 */
		void constructEqu(Eigen::MatrixXd variance, int var_type); //var_type: 0-> initial variance; 1 posterior variance

		/**
		 * @brief Construct GNSS equations using Ceres covariance results
		 *
		 * Builds linearized observation equations from Ceres solver covariance output.
		 * Integrates optimization results with GNSS measurement model for quality assessment and statistical analysis in post-processing validation.
		 */
		void constructEqu_fromCeres(Eigen::MatrixXd variance);
		
		/**
		 * @brief Precompute residuals and cache parameter data for GNSS marginalization
		 *
		 * Evaluates all GNSS residual factors and stores current parameter values.
		 * Prepares data structures for Schur complement computation in GNSS context.
		 */
		void preMarginalize();
		
		
		std::vector<ParameterBlockKey> ordered_addrs;   // 记录参数块顺序
		
		/**
		 * @brief Perform Schur complement marginalization for GNSS factors
		 *
		 * Executes GNSS-specific marginalization using Schur complement method.
		 * Partitions parameters into marginalized and retained sets, computes information matrices.
		 */
		void marginalize();
		
		
		int localSize(int size) const;

		/**
		 * @brief Get retained parameter blocks after GNSS marginalization
		 *
		 * Extracts GNSS parameter blocks to be kept (not marginalized) and updates their addresses using the provided shift mapping.
		 * Prepares data structures for GNSS prior factor construction in subsequent optimization.
		 */
		std::vector<double *> getParameterBlocks(std::unordered_map<ParameterBlockKey, double *> &addr_shift);
		/// manager the  factors 
		/// Combine all factor's jacobians and residuals into a large matrix
		std::vector<GNSSResidualBlockInfo *> factors;
		//t_gamb_manager* amb_manager = nullptr;  // 指向外部模糊度管理器
		int m = 0, n = 0;			
		//std::unordered_map<ParameterBlockKey, int> parameter_block_size;   /// Store all parameter addresses and sizes	总的变量块
		int sum_block_size = 0;
		//std::unordered_map<ParameterBlockKey, int> parameter_block_idx;    /// parameter index in Jacobians要被丢弃的变量块
		//std::unordered_map<ParameterBlockKey, double *> parameter_block_data;
		std::map<ParameterBlockKey, int> parameter_block_size;
		std::map<ParameterBlockKey, int> parameter_block_idx;
		std::map<ParameterBlockKey, double*> parameter_block_data;
		std::vector<int> keep_block_size;						/// global size
		std::vector<int> keep_block_idx;						/// local size
		std::vector<double *> keep_block_data;
		int n_size = 0;
		int para_size = 0;
		std::vector<std::map<ParameterBlockKey, std::vector<int>>> para_index;       ///for reconstruct jacobian

		// add for debug hwzhang
		Eigen::VectorXd postfit_qv_diag;
		int postfit_nobs = 0;
		int postfit_npar = 0;
		int postfit_dof_raw = 0;
		/// final residuals and jacobians 
		Eigen::MatrixXd linearized_jacobians;
		Eigen::VectorXd linearized_residuals;
		Eigen::MatrixXd weight;
		Eigen::VectorXd v_norm;
		Eigen::MatrixXd Qx;               /// parameter variance after optimization
		Eigen::MatrixXd Qx0;              /// initial variance
		double          sig_unit;		  ///< sig unit
		double          vtpv;
		const double eps = 1e-8;
		bool valid;
	};

	class LibGREAT_LIBRARY_EXPORT MarginalizationGNSSFactor : public ceres::CostFunction
	{
	public:

		/**
		 * @brief Construct GNSS marginalization prior factor
		 * Initializes parameter block sizes from GNSS marginalization info
		 * Sets up residual dimension for GNSS prior constraint
		 */
		explicit MarginalizationGNSSFactor(GNSSInfo* _gnss_info);

		/**
		 * @brief Evaluate GNSS marginalization prior factor
		 *
		 * Computes residuals and Jacobians for marginalized GNSS prior constraint.
		 * Handles linear error propagation from Schur complement marginalization results for GNSS-specific parameters.
		 */
		virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
		GNSSInfo* gnss_info = nullptr;
	};
}
#endif
