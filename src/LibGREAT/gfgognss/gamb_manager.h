#ifndef	AMBMANAGER_H
#define AMBMANAGER_H
/**
 * @file         gamb_manager.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        manager of ambiguity parameters
 * @version      1.0
 * @date         2025-01-01
 *
 * @copyright Copyright (c) 2025, Wuhan University. All rights reserved.
 *
 */
#include "gproc/gpvtflt.h"
#include "gexport/ExportLibGREAT.h"
namespace gfgomsf
{
	class t_gsat_map
	{
	public:
		t_gsat_map(GSYS sys, string sat_name, int id) :
			_gnss_system(sys), _sat_name(sat_name),_global_id(id)

		{
			DD_used = false;
		}			
		t_gsat_map(GSYS sys, string sat_name, double obs_time, int id) :
			_gnss_system(sys), _sat_name(sat_name), _start_obs_time(obs_time), _global_id(id)

		{
			DD_used = false;
			time_span.push_back(_start_obs_time);
		}
		GSYS _gnss_system;
		string _sat_name;
		double _start_obs_time;
		const int _global_id;
		vector<int> _amb_ids;   //stored the amb id;
		bool DD_used;
		vector<double> time_span;
		
	};
	class t_grover_msg
	{
	public:
		t_grover_msg()
		{

		}
		t_grover_msg(t_gtime _cur_time, Eigen::Vector3d _crd, Eigen::Vector3d _vel=Eigen::Vector3d::Identity()):
			cur_time(_cur_time),crd(_crd),vel(_vel)
		{
		}
	public:
		t_gtime cur_time;         //current epoch
		Eigen::Vector3d crd;      //position
		Eigen::Vector3d vel;      //velocity
	protected:		
		vector<t_gsatdata> _sat_data; //sat data
	};
	class t_gamb_per_ID
	{
	public:
		t_gamb_per_ID(pair<FREQ_SEQ, GOBSBAND> freq_band, GSYS sys, string prn, int sat_id, int amb_id, int start_rover, double initial_amb) :
			_freq_band(freq_band), _gnss_system(sys), _sat_prn(prn), _sat_global_id(sat_id), _amb_id(amb_id), _start_rover_count(start_rover),
			_initial_amb(initial_amb), _estimated_amb(initial_amb), _fixed_amb(-1.0), _is_fix(false),_is_slip(false), _is_estimated(false)
		{
			_rover_list.push_back(start_rover);
		}

		t_gamb_per_ID(GSYS sys, string prn, int sat_id, int amb_id, int start_rover, double initial_amb) :
	_gnss_system(sys), _sat_prn(prn), _sat_global_id(sat_id), _amb_id(amb_id), _start_rover_count(start_rover),
	_initial_amb(initial_amb), _estimated_amb(initial_amb), _fixed_amb(-1.0), _is_fix(false), _is_slip(false), _is_estimated(false)
		{
			_rover_list.push_back(start_rover);
		}
		int endRover();
		int startRover();
		void setBeg(t_gtime beg);
		void setEnd(t_gtime end);
		void setSlip();
		void addRover(int rover_count);
		double get_initial_value();
		double get_est_value();
		double is_estimated();
		void set_est_value(double est_value);
		void set_fixed_value(double fixed_value);
		void setFix();
		string getPRN();
		pair<FREQ_SEQ, GOBSBAND> getFB();
	public:
		GSYS _gnss_system;
		string _sat_prn;	
		int _sat_global_id = -1;
		pair<FREQ_SEQ, GOBSBAND> _freq_band;
		const  int _amb_id;							
		int _start_rover_count = 0;						
		vector<int> _rover_list;  //stored the rover index in sliding window
		bool _is_slip;
		double _initial_amb = 0.0;
		double _estimated_amb = -1.0;
		double _fixed_amb = -1.0;
		bool _is_estimated = false;
		bool _is_fix = false;
		t_gtime _beg;	  ///< begin time
		t_gtime _end;	  ///< end time	(end time corresponds to a new amb arc start)
	};
	class LibGREAT_LIBRARY_EXPORT t_gamb_manager
	{
	public:
		t_gamb_manager();
		t_gamb_manager(map<GSYS, map<FREQ_SEQ, GOBSBAND>> band_index);
		~t_gamb_manager();
		/**
		 * @brief Clear all ambiguity management states
		 * Resets satellite maps, ambiguity parameters and search indices
		 * Prepares for new processing session
		 */
		void clearState();
		/**
		 * @brief Remove satellite and update ambiguity states
		 * Handles satellite loss by clearing or updating affected ambiguities
		 * Maintains consistency in multi-epoch tracking
		 */
		void removeSat(const int &sat_global_id, const int &rover_index);
		/**
		 * @brief Slide window and update ambiguity tracking
		 * Removes oldest epoch data and shifts rover indices
		 * Maintains satellite-ambiguity consistency in window
		 */
		void slidingWindow();	
		/**
		 * @brief Get total ambiguity parameter count
		 */
		int getAmbCount();
		/**
		 * @brief Get start rover ID for ambiguity
		 * Returns first epoch index where ambiguity appears
		 */
		int getAmbStartRoverID(const int &amb_search_id);
		/**
		 * @brief Get end rover ID for ambiguity
		 * Returns last epoch index where ambiguity appears
		 */
		int getAmbEndRoverID(const int &amb_search_id);
		/**
		 * @brief Get ambiguity parameter vector
		 * Returns all current ambiguity values as vector
		 */
		Eigen::VectorXd getAmbVector();
		/**
		 * @brief Get estimated ambiguity value
		 */
		double getAmb(int amb_id);
		/**
		* @brief Get initial ambiguity value 
		*/
		double getInitialAmb(const int &amb_id);
		/**
		 * @brief Get ambiguities for marginalization
		 * Returns ambiguities ending at oldest epoch
		 */
		vector<int> getMarginAmb();
		/**
		 * @brief Get current window ambiguities
		 * Returns all active ambiguities in window
		 */
		vector<int> getCurWinAmb();
		/**
		 * @brief Update ambiguity parameter value
		 * Sets new estimated float ambiguity
		 */
		void updateAmb(int amb_id, double value);
		/**
		 * @brief Generate ambiguity search index
		 * Creates mapping from satellite-frequency to ambiguity ID
		 */
		void generateAmbSearchIndex();	
		/**
		 * @brief Get current epoch satellites
		 * Updates satellite list for latest epoch
		 */
		int get_last_epoch_sats(double time);
		/**
		 * @brief Get satellite global ID
		 * Returns satellite ID from PRN name
		 */
		int get_sat_id(string prn);
		/**
		 * @brief Set estimated ambiguities from vector
		 * Updates all ambiguity values from solution vector
		 */
		void setEstAmb(const Eigen::VectorXd &x);
		/**
		 * @brief Get ambiguity index by satellite and frequency
		 */
		int getAmbSearchIndex(const pair<int, FREQ_SEQ> &sat_freq);	
		/**
		 * @brief Add new satellite with ambiguities
		 * Initializes satellite tracking and creates ambiguity parameters
		 */
		bool addNewSat(const t_gtime & cur_time, const int &rover_index, const int &sat_index, int &amb_index, const t_gsatdata &sat_data, t_gallpar params);
		/**
		 * @brief Extend satellite tracking to new epoch
		 * Updates all ambiguities for given satellite to include new rover epoch
		 */
		void addRover(string sat_name,const int &rover_index);
		/**
		 * @brief Add rover epoch with timestamp
		 * Extends satellite tracking time span and ambiguity rover lists
		 */
		void addRover(double time, string sat_name, const int &rover_index);
		/**
		 * @brief Create ambiguity parameters for new satellite
		 * Initializes single-difference ambiguities for all frequencies
		 * Updates satellite-ambiguity mapping and global ID tracking
		 */
		bool addAmb(const t_gtime & cur_time, vector<t_gpar> amb_para, const  GSYS & gnss_system, int rover_index, const int &sat_index, int &amb_index);
		/**
		 * @brief Add ambiguity parameter to parameter set
		 */
		void addGpara(t_gallpar &params, const int &idx, double value);
		/**
		 * @brief Mark satellite as used in DD
		 */
		void checkDDSat(const int &sat_id);
		/**
		 * @brief Remove unused satellites and ambiguities
		 * Cleans up satellites not used in double-difference processing
		 * Maintains only actively tracked ambiguities
		 */
		void updateAmb();
		/**
		 * @brief Check if satellite is currently tracked
		 */
		bool is_sat_tracking(string prn);
	public:
		vector<int> ambiguity_ids;			
		vector<pair<string, int>> cur_sats;
		map<string, int> removed_sat;	
		vector<int> removed_ambs;
	protected:
		map<int, shared_ptr<t_gamb_per_ID>> _ambiguity;		// All ambiguity information within the current window (corresponding satellite and index, frequency, ambiguity index, whether cycle slip, whether fixed)
		map<int, shared_ptr<t_gsat_map>> _sat_map;			// All satellite information within the current window (satellite index and corresponding ambiguity index)
		map<GSYS, map<FREQ_SEQ, GOBSBAND>> _band_index;
		map<par_type, FREQ_SEQ> ambtype_list = {
				{par_type::AMB_L1,FREQ_1},
				{par_type::AMB_L2,FREQ_2},
				{par_type::AMB_L3,FREQ_3},
				{par_type::AMB_L4,FREQ_4},
				{par_type::AMB_L5,FREQ_5} };
		map<FREQ_SEQ, par_type> freq_ambtype_list = {
				{FREQ_1, par_type::AMB_L1},
				{FREQ_2, par_type::AMB_L2},
				{FREQ_3, par_type::AMB_L3},
				{FREQ_4, par_type::AMB_L4},
				{FREQ_5, par_type::AMB_L5} };
		map<pair<int,FREQ_SEQ>,int> search_index;		
	};
}
#endif
