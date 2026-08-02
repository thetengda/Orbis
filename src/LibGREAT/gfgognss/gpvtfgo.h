#ifndef GPVTFGO_H
#define GPVTFGO_H
/**
 * @file         gpvtfgo.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        main code of GNSS RTK besed on factor graph optimization
 * @version      1.0
 * @date         2025-01-01
 *
 * @copyright Copyright (c) 2025, Wuhan University. All rights reserved.
 *
 */

#include "gfgo/gfgo_estimator.h"
#include "gfgo/gutility.h"
#include "gproc/gpvtflt.h"
#include "gutils/gtimesync.h"
#include "gamb_manager.h"
#include "gexport/ExportLibGREAT.h"
#include "gfactor/carrierphase_DD_factor.h"
#include "gfactor/pseudorange_DD_factor.h"
#include "gfactor/gnss_info_factor.h"
#include "gfgo/gprecisebiasFGO.h"


#include "gfgognss/gambIF_manager.h"
#include "gfactor/initial_factor.h"
#include "gfactor/carrierphase_IF_factor.h"
#include "gfactor/pseudorange_IF_factor.h"
#include "gfactor/multi_carrierphase_IF_factor.h"
#include "gfactor/multi_pseudorange_IF_factor.h"
#include "gfactor/random_walk_factor.h"

#include <map>
#include <set>


namespace gfgomsf
{
	using namespace gnut;
	using namespace gfgo;
	
	/**
	*@brief  class for graph optimzation based GNSS pvt, derive from gpvtflt
	*/
	class LibGREAT_LIBRARY_EXPORT t_gpvtfgo: virtual public t_gpvtflt,
		                                   virtual public t_gfgo
	{
	public:	
		
		explicit t_gpvtfgo(string site, string site_base, t_gsetbase* gset, std::shared_ptr<spdlog::logger> spdlog, t_gallproc* allproc);/// initialization of parameters for factor graph based pvt
		/** @brief default destructor. */
		virtual ~t_gpvtfgo();
	public:
		/** @brief processBatch. */
		virtual int processBatch(const t_gtime &beg, const t_gtime &end, bool prtOut) override;
		/**
		* @brief One epoch processing
		* @note virtual function, Wangxb Created for fgopvt
		* @param[in] Now	current time
		* @param[in] data_rover
		* @param[in] data_base
		* @return
			   -1,failed;
				0,float;
				1,fixed.
		*/		
		int processWindow(const t_gtime &now, vector<t_gsatdata> *data_rover = NULL, vector<t_gsatdata> *data_base = NULL);
		
		/**
		* @brief Clear optimization window
		* @Clear all parameters in the window, including _DD_msg, _para_CRD, _para_amb, _amb_manager and so on.
		*/
		void clearWindow();
		virtual void publish_foat();
		
	protected:
		class DDEquMsg
		{
		public:
			t_gtime time;
			t_gsatdata  rover_ref_sat;
			t_gsatdata  rover_nonref_sat;
			t_gsatdata  base_ref_sat;
			t_gsatdata  base_nonref_sat;
			GOBSTYPE  obs_type;
			FREQ_SEQ  freq;
			GOBSBAND band;
			string base_site;
			string rover_site;
			string ref_sat;
			string nonref_sat;
			int ref_sat_global_id = -1;
			int nonref_sat_global_id = -1;
		public:
			explicit DDEquMsg(const t_gsatdata& _ref_sat, const t_gsatdata& _nonref_sat, const GOBSTYPE& _obs_type, const FREQ_SEQ& _freq);

		};		


		//for PPP IF:
		class IFEquMsg
		{
		public:
			t_gtime time;
			t_gsatdata satdata;
			GOBSTYPE obs_type;
			FREQ_SEQ freq1;
			FREQ_SEQ freq2;
			GOBSBAND band1;
			GOBSBAND band2;
			string site;
			int sat_id = -1;		//global sat id
		public:
			explicit IFEquMsg(const t_gsatdata& _satdata, const GOBSTYPE& _obs_type, const FREQ_SEQ& _freq1, const FREQ_SEQ& _freq2, const GOBSBAND& _band1, const GOBSBAND& _band2) :
				satdata(_satdata), obs_type(_obs_type), freq1(_freq1), freq2(_freq2), band1(_band1), band2(_band2) {
			}
		};



		vector<DDEquMsg> _DD_msg;
		vector<vector<DDEquMsg>> _vDD_msg;
		//for PPP IF:
		vector<IFEquMsg> _IF_msg;
		vector<vector<IFEquMsg>> _vIF_msg;

		vector<pair<pair<string, int>, pair<FREQ_SEQ, GOBSTYPE>>> _gnss_obs_index;
		std::vector<PseudorangeDDFactor*> window_pseudo_factors;
		std::vector<CarrierphaseDDFactor*> window_carrierphase_dd_factors;
		int              _rover_count = -1;
		int              _global_sat_id = -1;
		int              _global_amb_id = -1;
		int              _ppp_candidate_global_sat_id = -1;
		int              _ppp_candidate_global_amb_id = -1;
		std::set<std::string> _pending_ppp_slips;
		std::set<std::string> _candidate_ppp_slips;
		std::map<std::string, t_gtime> _last_ppp_phase_epoch;
		double           _loss_func_value = 3;
		double           _gtime_interval = 1.0;
		shared_ptr<t_gamb_manager>  _amb_manager;
		//for PPP IF:
		shared_ptr<t_gambIF_manager>  _ambIF_manager;
		Eigen::Vector3d        _Pos[GWINDOW_SIZE + 1];							///< position of rover
		Eigen::Vector3d        _Vel[GWINDOW_SIZE + 1];


		//for PPP IF:
		double _trp_ini = 0.0;
		double _isb_GAL_ini = 0.0;
		double _isb_BDS_ini = 0.0;
		double _isb_GLO_ini = 0.0;
		double _clk[GWINDOW_SIZE + 1];
		double _trp[GWINDOW_SIZE + 1];
		double _isb_GAL[GWINDOW_SIZE + 1];
		double _isb_BDS[GWINDOW_SIZE + 1];
		double _isb_GLO[GWINDOW_SIZE + 1];
		bool _lost_isb_GAL[GWINDOW_SIZE + 1];
		bool _lost_isb_BDS[GWINDOW_SIZE + 1];
		bool _lost_isb_GLO[GWINDOW_SIZE + 1];

		ofstream _output_float_solution;//add for MultiWindow
		ofstream _output_estimator_info;//add for MultiWindow

		///< velocity of rover
		double          _headers[GWINDOW_SIZE + 1];						    ///< header info of rover	
		t_grover_msg *_rover_window[GWINDOW_SIZE + 1];                         ///< raw data of rover		                 
		vector<vector<t_gsatdata>>  _win_base_data;                         ///<base sat data
		t_gallpar _para_window[GWINDOW_SIZE + 1];
		t_gallpar _all_para_win;
		map<string, int> cur_sat_prn;  // current satellite list
		vector<string> outlier_sats;
		vector<pair<string, int>> _removed_sats;	
		bool _batch_remove = false;
		bool _initial_prior = true;
		std::unordered_map<long, int> amb_idx;	



	protected:
		double _gclk_rover, _gclk_base;		
		int _obs_level = 3;
		t_gtriple _gcrd_base;
		t_gprecisebiasFGO *_gbias_model = nullptr;					                    ///< baise model
		GNSSInfo *_last_gnss_info = nullptr;                             ///save gnss equ
		GNSSInfo *_last_gnss_marginalization_info = nullptr;             ///for gnss marginalization
		vector<double *> _last_gnss_marginalization_para_blocks;
		std::unordered_map<long, int> all_parameter_block; /// block to store info of last marginalized parameters
	protected:	
		void _prepare_equ();
		virtual bool _update_all_equ();
		bool _get_satdata(const int epoch_id, const string & sat, const string & rec, t_gsatdata & satdata);
		bool _set_satdata(const int epoch_id, const string & sat, const string & rec, const t_gsatdata & satdata);
		/**
		 * @brief Get processed GNSS observation data for rover and base stations
		 *
		 * Retrieves satellite observation data, applies DCB corrections, and analyzes
		 * frequency availability for both rover and base stations.
		 * @param data_rover Output parameter for rover station satellite data
		 * @param data_base Output parameter for base station satellite data (if applicable)
		 * @return true if data retrieval and processing successful for required stations
		 * @return false if no observation data found for required stations
		 */
		bool _get_gdata(const t_gtime &now, vector<t_gsatdata> *data_rover = NULL, vector<t_gsatdata> *data_base = NULL);
		/** Set GLONASS FDMA channels from observation headers or navigation data. */
		void _set_glonass_channels(vector<t_gsatdata>& data, const t_gtime& epoch);
		/**
		 * @brief Initialize state parameters for GNSS positioning
		 * Performs ambiguity prediction, coordinate estimation, and state vector initialization
		 */
		virtual void _get_initial_value(const t_gtime& runEpoch);
		/**
		 * @brief Initialize state for current epoch with slip detection
		 * Manages sliding window, detects cycle slips, and updates ambiguity parameters
		 */
		virtual void _set_initial_value(const t_gtime& runEpoch);

		//void _set_initial_value_beta(const t_gtime& runEpoch);//add for MultiWindow
		/**
		 * @brief Create temporary parameter set with base station coords and clock offsets
		 * Adds base station coordinates and receiver clock parameters to parameter set
		 */
		bool _gtemp_params(t_gallpar &params, t_gallpar &params_temp);
		/**
		 * @brief Set receiver coordinates and clock offsets
		 * Stores base station position and receiver clock parameters for processing
		 */
		void _set_rec_info(const t_gtriple &xyz_base, double clk_rover, double clk_base);	
		/**
		 * @brief Generate double-difference observations
		 * Creates double-difference combinations between reference and non-reference satellites
		 */
		int _combine_DD();
		/**
		 * @brief Select reference satellite and form DD combinations
		 * Chooses highest elevation satellite as reference and creates double-difference pairs with other satellites for all GNSS systems and frequencies
		 */
		int _combine_IF();

		void _select_ref_sat();	
		/**
		 * @brief Complete DD message with base station data and frequency information
		 * Populates the double-difference message with corresponding base station observations and frequency band details. Verifies common satellite visibility between rover and base.
		 *
		 * @param dd_msg Double-difference message to be completed (input/output)
		 * @param base_sat_data Base station satellite observation data
		 * @return true if DD data successfully prepared, false if common visibility check fails
		 */
		bool _get_DD_data(DDEquMsg &dd_msg, vector<t_gsatdata> base_sat_data);
		/**
		 * @brief Remove outlier satellite from current processing
		 * Eliminates outlier satellite from DD observations and ambiguity management
		 * @param outlier Pair of satellite name and global ID to be removed
		 * @return true if removal successful, false if insufficient satellites remain
		 */
		bool _remove_outlier_sat(const pair<string, int> & outlier);
		/**
		 * @brief Prepare data for ambiguity resolution
		 * Transfers linearized GNSS observation equations to filter for parameter estimation
		 * Sets up design matrix, weight matrix, residuals and covariance matrices
		 */
		//virtual bool _remove_outlier_sat_beta(const pair<string, int> & outlier);
		bool _pre_amb_resolution();	

	protected:
		/**
		 * @brief Add GNSS prior constraints to optimization problem
		 * Integrates marginalized GNSS information from previous epoch as prior factor
		 * Maintains temporal consistency in sliding window optimization
		 */
		void _prior_factor(ceres::Problem &problem);
		/**
		 * @brief Convert optimized parameters back to state vectors, virtual function
		 * Updates rover positions and ambiguity parameters from optimization results
		 * Transforms double array values to Eigen vectors and ambiguity states
		 */
		virtual void _double_to_vector() override;
		/**
		 * @brief Convert state vectors to double arrays for optimization, virtual function
		 * Prepares rover coordinates and ambiguity parameters for ceres solver
		 * Transforms Eigen vectors and ambiguity states to double arrays
		 */
		virtual void _vector_to_double() override;	
		/**
		 * @brief Perform nonlinear optimization for GNSS positioning and ambiguity resolution
		 *
		 * Executes iterative optimization with outlier detection and removal. Sets up Ceres problem with:
		 * - Prior factors from marginalization
		 * - GNSS double-difference factors (pseudorange and carrier phase)
		 * - Ambiguity constraints and initial position priors
		 * - Automatic outlier detection and iterative re-optimization
		 *
		 * @return Optimization status code
		 */
		virtual int _optimization();

		virtual int _optimization_PPP();
		/**
		 * @brief Roll back a rejected pure-PPP candidate epoch
		 *
		 * Removes only the newest PPP node and its ambiguity associations while
		 * preserving the previously accepted FGO window and marginalization prior.
		 */
		void _rollback_current_ppp_node();

		// Track carrier continuity independently of accepted PPP graph nodes.
		// A pending event is consumed only after a new ambiguity arc is retained
		// by a successfully accepted node.
		void _record_ppp_cycle_slips(const std::vector<t_gsatdata>& epoch_data);
		void _commit_ppp_cycle_slips();
		bool _ppp_candidate_has_phase(const std::string& sat_name) const;
		/**
		 * @brief Slide processing window forward by one epoch
		 * Shifts window data when full, removes oldest epoch, updates ambiguity states and _para_CRD
		 */
		virtual void _slide_window();
		/**
		 * @brief Construct marginalization information for sliding window optimization
		 * Prepares marginalization factors including GNSS observations and ambiguity constraints
		 * Maintains consistency when removing oldest epoch from the processing window
		 */
		void _marginalization();


		void _marginalization_PPP();
		/**
		 * @brief Detect GNSS observation outliers using normalized residuals
		 * Identifies and flags satellites with excessive residuals for removal
		 * @param outlier Output pair of satellite name and ID to be removed
		 * @return Index of detected outlier, -1 if none found
		 */
		int  _gobs_outlier_detection(pair<string, int> & outlier);
		/**
		 * @brief Check if satellite is marked as outlier
		 */
		bool _check_outlier(const string &sat);	

		vector<double*> _parameter_blocks;

		double posterior_std_X, posterior_std_Y, posterior_std_Z;

		/**
		 * @brief Perform posteriori validation after optimization
		 * Computes covariance matrices and validates solution quality
		 * Updates GNSS information with residuals and statistical metrics
		 */
		void _posteriori_test(ceres::Problem &problem);
		//void _posteriori_test_PPP(ceres::Problem& problem);
		void _posteriori_test_PPP(ceres::Problem & problem);
	};
}
#endif
