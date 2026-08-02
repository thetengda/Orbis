#ifndef GFGO_GINS_H
#define GFGO_GINS_H
#define NOMINMAX
/**
 * @file         gfgo_gins.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        Class for factor graph-based GNSS/SINS integrated solution
 * @version      1.0
 * @date         2025-11-04
 *
 * @copyright Copyright (c) 2025, Wuhan University. All rights reserved.
 *
 */

#include"gall/gallpar.h"
#include "gproc/gpvtflt.h"
#include "gmsf/gintegration.h"
#include "gins/gins.h"
#include "gfgomsf/gpara.h"
#include "gfgognss/gpvtfgo.h"
#include "gexport/ExportLibGREAT.h"
#include <ceres/ceres.h>
#include <Eigen/Dense>
#include <iostream>
#include <set>
#include <vector>
#include "gmsf/gpublish.h"
#include "gfactor/gpose_local_parameterization.h"




#include "gfactor/carrierphase_IF_integration_factor.h"
#include "gfactor/pseudorange_IF_integration_factor.h"
#include "gfactor/multi_carrierphase_IF_integration_factor.h"
#include "gfactor/multi_pseudorange_IF_integration_factor.h"


namespace gfgomsf
{
	using namespace gnut;
	using namespace gfgo;
	/**
	*@brief   Class for factor graph optimization based multi-sensor fusion operation
	*/
	class t_gnss_node
	{
	public:
		t_gnss_node() {};
		double t = 0.0;
		bool is_key_node;	///<Time		
		Eigen::Matrix3d R;														///<Rotation
		Eigen::Vector3d T;														///<Translation
		IntegrationBase* pre_integration = nullptr;								///<Parameter integration
		///<Whether is key frame
	};

	class LibGREAT_LIBRARY_EXPORT t_gfgo_gins : virtual public t_gsinskf,
		virtual public t_gpvtfgo,
		virtual public t_gfgo,
		virtual public t_gpara
		//virtual public t_gfgo_slam
	{
	public:

		/**
		 * @brief GNSS/INS tight coupling integration constructor
		 * Initializes multi-inheritance navigation system with GNSS PPP, FGO, INS KF components
		 * Sets up optimization states, IMU pre-integration buffers, and result publishing
		 */
		t_gfgo_gins(string site, string site_base, t_gsetbase* gset, std::shared_ptr<spdlog::logger> spdlog, t_gallproc* allproc);/// initialization of parameters for factor graph based multi_sensor fusion solution

		virtual void init() {}				/// the same name function of t_gsinskf
		virtual void set_Hk() {}			/// the same name function of t_gsinskf

		/**
		* @brief Get imu data
		*/
		void get_imu(great::t_gimudata* _imu);

		/**
		 * @brief Execute GNSS/INS integrated processing in batch mode
		 *
		 * Main processing loop for tightly-coupled GNSS/INS integration. Handles:
		 * - IMU and GNSS data synchronization and time alignment
		 * - Initial coarse alignment using GNSS position/velocity
		 * - IMU mechanization and gravity updates
		 * - GNSS measurement processing and factor graph optimization
		 * - Real-time state estimation and result publishing
		 *
		 * @param beg Start time of processing batch
		 * @param end End time of processing batch
		 * @param beg_end Flag indicating batch start/end
		 * @return Processing status code
		 */
		int ProcessBatch(const t_gtime& beg, const t_gtime& end, bool beg_end); 
		

	protected:

		/**
		 * @brief Initialize GNSS/INS integration processing
		 * Sets up time ranges, validates IMU data availability, and initializes INS filter
		 * Synchronizes GNSS and IMU data streams for integrated processing
		 */
		virtual int _init();/// initialize imu 		
		/**
		 * @brief Validate and synchronize GNSS observation time
		 *
		 * Checks time alignment between IMU and GNSS data, loads observations,
		 * applies DCB corrections, and validates data availability for both
		 * rover and base stations. Handles cycle slip detection and time synchronization.
		 */
		bool _time_valid(t_gtime gt, t_gtime imut);
		/**
		 * @brief Virtual function, determine available measurement types for current epoch
		 *
		 * Identifies valid GNSS measurements based on time alignment and motion state. Updates motion type mapping and clears yaw when needed.
		 */
		virtual int _getMeas();								/// obtain all observation types 

		/**
		 * @brief Get GNSS position solution for alignment
		 * Processes single GNSS epoch, retrieves position/velocity solution,
		 * and updates INS states for initial alignment. Validates solution quality.
		 */
		MEAS_TYPE _getPOS(t_gposdata::data_pos& pos);		/// obtain gnss position measurement type	

		/**
		 * @brief Process GNSS measurements for tight coupling
		 * Retrieves GNSS observations and prepares factor graph constraints
		 * Sets GNSS factor availability flag for integration
		 */
		int _gnss_processing();								/// process gnss data 

		virtual void _double_to_vector() {};					/// unfished

		virtual void _vector_to_double() {};					/// unfished	

		/**
		 * @brief Update gravity vector in navigation frame
		 * Computes gravity components in either ECEF or NED frame, default ECEF(_E_F)
		 */
		void _gravity_update();

		/**
		 * @brief Prepare GNSS measurements for tight coupling
		 *
		 * Sets up GNSS data processing including coordinate initialization,
		 * double-difference combination, and ambiguity resolution preparation.
		 * Validates satellite visibility and data quality for integration.
		 */
		MEAS_TYPE _get_gnss_measurements(const t_gtime& runEpoch);

		/**
		 * @brief Initialize GNSS state parameters for new epoch
		 * Performs ambiguity synchronization, coordinate and ambiguity prediction
		 * Sets up initial state vector for GNSS factor graph optimization
		 */
		virtual void _get_initial_value(const t_gtime& runEpoch) override;

		/**
		 * @brief Initialize GNSS state for new epoch with slip detection
		 *
		 * Manages sliding window, detects cycle slips and satellite losses,
		 * updates ambiguity states, and maintains satellite tracking consistency
		 * across epochs for GNSS/INS tight coupling
		 */
		virtual void _set_initial_value(const t_gtime& runEpoch) override;

		/**
		 * @brief Store current INS state in sliding window frame
		 *
		 * Saves current position, velocity, attitude and bias states into the factor graph optimization window. Handles initial GNSS position alignment
		 */
		void _set_frame_pose(const int& framecount);

		/**
		 * @brief Perform IMU pre-integration between two time epochs
		 *
		 * Integrates IMU measurements between GNSS epochs to compute relative motion
		 * Manages IMU buffer, handles bias states, and builds pre-integration factors
		 */
		bool _imu_pre_integration(double t1, double t2);

		/**
		 * @brief Execute GNSS/INS tight coupling optimization
		 *
		 * Manages the complete factor graph optimization pipeline including:
		 * - GNSS factor integration and pre-integration management
		 * - Initialization phase handling and solver state transitions
		 * - Nonlinear optimization, ambiguity resolution, and marginalization
		 * - Sliding window maintenance for real-time processing
		 */
		int _gins_processing();

		/**
		 * @brief Execute GNSS/INS tight coupling optimization
		 *
		 * Performs iterative nonlinear optimization with outlier rejection for GNSS/INS integration.
		 * Sets up factor graph with IMU pre-integration, GNSS DD observations, ambiguity constraints,
		 * and marginalization priors. Handles robust estimation and quality validation.
		 */
		void _gins_optimization();

		void _gins_optimization_PPP();
		void _capture_alignment_gnss_prior();
		bool _make_alignment_gnss_prior(
			ceres::CostFunction*& factor,
			vector<double*>& parameter_blocks,
			vector<int>& drop_set,
			bool marginalizing) const;

		/**
		 * @brief Perform posterior covariance analysis for GNSS/INS solution
		 *
		 * Computes parameter covariance matrices and validates solution quality after optimization.
		 * Handles GNSS double-difference observations, constructs covariance blocks, and
		 * performs statistical analysis for solution validation and outlier detection.
		 */
		void _gins_posteriori_test(ceres::Problem& problem);
		void _gins_posteriori_test_PPP(ceres::Problem& problem);


		/**
		 * @brief Perform marginalization for GNSS/INS sliding window
		 *
		 * Manages the marginalization process when the sliding window reaches maximum size.
		 * Handles IMU factors, GNSS double-difference observations, static constraints,
		 * and ambiguity parameters. Maintains consistency while removing oldest states.
		 */
		void _gins_marginalization();
		void _gins_marginalization_PPP();
		/**
		 * @brief Convert optimized parameters back to state vectors
		 * Updates rotation matrices, positions, velocities, biases and ambiguities
		 * from double array optimization results to Eigen vector states
		 */
		void _gins_double_to_vector();

		/**
		 * @brief Convert state vectors to double arrays for optimization
		 * Prepares pose, velocity, bias and ambiguity parameters for ceres solver
		 * Transforms Eigen vectors and quaternions to double array format
		 */
		void _gins_vector_to_double();

		/**
		 * @brief Slide GNSS/INS processing window forward
		 *
		 * Manages sliding window mechanism by shifting states, IMU data, and GNSS observations.
		 * Removes oldest epoch data while maintaining window size constraints and buffer consistency.
		 */
		void _slide_gins();

		/**
		 * @brief Buffer IMU measurements for pre-integration
		 * Stores accelerometer and gyroscope data in time-synchronized buffers
		 */
		void _input_imu_data(const double& t, const vector<Eigen::Vector3d>& wm, const vector<Eigen::Vector3d>& vm);

		/**
		 * @brief Extract IMU measurements between two time epochs
		 * Retrieves and buffers accelerometer and gyroscope data for pre-integration
		 * between specified time intervals, maintaining measurement synchronization
		 */
		bool _get_imu_interval(double t0, double t1, vector<pair<double, Eigen::Vector3d>>& accVector,
			vector<pair<double, Eigen::Vector3d>>& gyrVector);

		/**
		 * @brief Execute GNSS ambiguity resolution
		 */
		void _gnss_amb_resolution();

		/**
		 * @brief Update INS states with optimization results
		 *
		 * Feeds back optimized pose, velocity, and bias states to INS mechanization
		 * Handles robust position fixing and maintains attitude consistency between optimization framework and inertial navigation system
		 */
		void _gins_feedback();

		// Undo only the newest PPP/INS GNSS node.  This is the PPP counterpart of
		// the RTK/INS FGO slideNew() failure path: the INS pre-integration is kept
		// so that it can span the failed GNSS epoch.
		void _rollback_ppp_node(bool gnss_node_created);

		// Cycle-slip tracking is persistent across rejected PPP candidates.  A
		// pending slip is consumed only after the corresponding new ambiguity arc
		// is retained by a successfully accepted graph node.
		void _record_ppp_cycle_slips(const std::vector<t_gsatdata>& epoch_data);
		void _commit_ppp_cycle_slips();

		/**
		 * @brief Clear all GNSS processing states and data
		 * Resets ambiguity management, satellite tracking, and optimization parameters
		 * Prepares for new GNSS processing session after time gaps or failures
		 */
		void clearGNSSmsg();

		/**
		 * @brief Write navigation solution results to output files
		 *
		 * Outputs INS states, GNSS quality metrics, ambiguity status, and statistical
		 * information to log files. Formats position, velocity, DOP values, and
		 * ambiguity resolution status for analysis and visualization.
		 */
		virtual void write(string info = "");

		/**
		 * @brief Get robust fixed position solution
		 * Returns fixed ambiguity position if solution meets quality thresholds
		 * Validates satellite count and ambiguity resolution ratio
		 */
		bool _getRobustFixedPosition(Eigen::Vector3d&);
		
		/**
		* @brief write trajectory to KML file (xml need to set)
		*/
		int prt_ins_kml();

		


	protected:

		enum SolverFlag
		{
			INITIAL,
			NON_LINEAR
		};
		t_gtime _gnss_beg;									///< begin time of gnss
		t_gtime _gnss_end;									///< end time of gnss
		t_gtime _gnss_crt;
		MEAS_TYPE _msf_flag;
		t_gtime _imu_beg;										///< begin time of imu
		t_gtime _imu_end;///< current time of gnss
		t_gtime _imu_crt, _img_crt, _laser_crt;
		int _fgoflag = 0;										///< coefficient used to control frequency
		bool _c_gnss_factor;								///< sign of whether the gnss factor is constructed
		int _opt_valid ;
		std::set<std::string> _pending_ppp_slips;
		std::set<std::string> _candidate_ppp_slips;
		std::map<std::string, t_gtime> _last_ppp_phase_epoch;
		int _ppp_candidate_global_sat_id = -1;
		int _ppp_candidate_global_amb_id = -1;
	protected:
		double _cur_node_time = -1;
		//int _rover_count = -1;
		double _last_pre_integration_time = -1;             ///last pre integration  time
		bool _opt_flag = false;
		Eigen::Vector3d _initial_pos;
		bool _initial_gnss_pos = true;

		Eigen::Vector3d        _Ps[(WINDOW_SIZE + 1)];							///< position of frame
		Eigen::Vector3d        _Vs[(WINDOW_SIZE + 1)];							///< velocity of frame
		Eigen::Matrix3d        _Rs[(WINDOW_SIZE + 1)];							///< rotation of frame
		Eigen::MatrixXd _global_variance;
		Eigen::Vector3d        _Bas[(WINDOW_SIZE + 1)];						///< bias of acceleration
		Eigen::Vector3d        _Bgs[(WINDOW_SIZE + 1)];
		IntegrationBase* _pre_integrations[(WINDOW_SIZE + 1)];
		vector<double> _dt_buf[(WINDOW_SIZE + 1)];
		vector<Eigen::Vector3d> _linear_acceleration_buf[(WINDOW_SIZE + 1)];	///< buffer to store linear_acceleration info
		vector<Eigen::Vector3d> _angular_velocity_buf[(WINDOW_SIZE + 1)];
		IntegrationBase* _tmp_pre_integration=nullptr;
		map<double, t_gnss_node> _all_gnss_node;
		SolverFlag _solver_flag;
		Eigen::Vector3d _acc_0, _gyr_0;
		//double          _headers[GWINDOW_SIZE + 1];
		double         _headers[(GWINDOW_SIZE + 1)];

		// debug hwzhang
		double _clk_prior_mean[(GWINDOW_SIZE + 1)] = {};

		// debug hwzhang
		// Fixed prior of the very first INS state.  The mean must never follow the
		// current estimate; its information is transferred by the first marginalization.
		bool _fixed_initial_state_prior_valid = false;
		Eigen::Vector3d _fixed_initial_P = Eigen::Vector3d::Zero();
		Eigen::Quaterniond _fixed_initial_Q = Eigen::Quaterniond::Identity();
		Eigen::Vector3d _fixed_initial_V = Eigen::Vector3d::Zero();
		Eigen::Vector3d _fixed_initial_Ba = Eigen::Vector3d::Zero();
		Eigen::Vector3d _fixed_initial_Bg = Eigen::Vector3d::Zero();
		Eigen::Matrix<double, 6, 6> _fixed_initial_pose_sqrt_info =
			Eigen::Matrix<double, 6, 6>::Identity();
		Eigen::Matrix<double, 9, 9> _fixed_initial_speed_bias_sqrt_info =
			Eigen::Matrix<double, 9, 9>::Identity();
		bool _alignment_gnss_prior_valid = false;
		vector<par_type> _alignment_gnss_prior_types;
		Eigen::VectorXd _alignment_gnss_prior_mean;
		Eigen::MatrixXd _alignment_gnss_prior_covariance;

		bool _first_imu = false;
		int mar_reset;
		map<t_gtime, vector<t_gsatdata>> _map_basedata;
		t_gpublish _publish_result;
		map<t_gtime, t_gallpar> _map_param;
		map<t_gtime, MOTION_TYPE> _map_motion;
		queue<pair<double, Eigen::Vector3d>> _accBuf;					///< buffer to  store acceleration measurement
		queue<pair<double, Eigen::Vector3d>> _gyrBuf;
		map<t_gtime, pair<Eigen::Vector3d, Eigen::Quaterniond>> _map_pose;
		IMUState _imu_state;

	};



}

#endif
