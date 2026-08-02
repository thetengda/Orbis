/**
 * @file         gpvtfgo.cpp
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        main code of GNSS RTK besed on factor graph optimization
 * @version      1.0
 * @date         2025-01-01
 *
 * @copyright Copyright (c) 2025, Wuhan University. All rights reserved.
 *
 */

#include "gpvtfgo.h"
#include <cmath>
#include <fstream>
#include "gmodels/gprecisebiasGPP.h"
#include "gfactor/ginitial_pose_factor.h"
#include <gutils/gcommon.cpp>

namespace
{
    double graph_interval_random_walk_q(gnut::t_randomwalk* model, double graph_dt)
    {
        const double stochastic_dt = std::fabs(model->get_dt());
        const double model_q = model->getQ();
        if (!std::isfinite(graph_dt) || graph_dt <= 0.0 ||
            !std::isfinite(stochastic_dt) || stochastic_dt <= 1e-9)
        {
            return model_q;
        }
        return model_q * graph_dt / stochastic_dt;
    }
}


gfgomsf::t_gpvtfgo::t_gpvtfgo(string site, string site_base, t_gsetbase * gset, std::shared_ptr<spdlog::logger> spdlog, t_gallproc * allproc):
t_gspp(site, gset, spdlog),
t_gpvtflt(site, site_base, gset, spdlog, allproc),
t_gfgo(gset),
t_gfgo_para(gset) {
	/*t_gbiasmodel *precise_bias(new t_gprecisebiasGPP(_allproc, _spdlog, gset));
	_gbias_model = precise_bias;*/
	_gbias_model = new t_gprecisebiasFGO(_allproc, _spdlog, gset);
	for (int i = 0; i <= gwindow_size; i++)
	{
		_Pos[i].setZero();
		_Vel[i].setZero();
		//PPP
		if (!_isBase)
		{
			_clk[i] = 0.0;
			_trp[i] = 0.0;
			_isb_GAL[i] = 0.0;
			_isb_BDS[i] = 0.0;
			_isb_GLO[i] = 0.0;
			_lost_isb_GAL[i] = false;
			_lost_isb_BDS[i] = false;
			_lost_isb_GLO[i] = false;
		}
		_headers[i]=0.0;
		_rover_window[i] = nullptr;
		_para_window->delAllParam();
	}
	_win_base_data.resize(gwindow_size+1);
	cur_sat_prn.clear();
	_all_para_win.delAllParam();
	_DD_msg.clear();
	_vDD_msg.clear();
	//add zhang
	if (!_isBase)//for PPP IF
	{
		_IF_msg.clear();
		_vIF_msg.clear();
		shared_ptr<t_gambIF_manager> ambIF_m(new t_gambIF_manager(_band_index));
		_ambIF_manager = ambIF_m;
	}
	else
	{
		shared_ptr<t_gamb_manager> amb_m(new t_gamb_manager(_band_index));
		_amb_manager = amb_m;
	}

	_pos_constrain = false;
	_cntrep = 0;
	int sign = 1;
	if (!_isBase)
	{
		// Keep a pure-PPP window across at most ten missing graph intervals.
		// For 1 Hz data this is 10 seconds; RTK keeps its original behavior.
		if (_sampling > 1)
			_gtime_interval=int(sign * _sampling * 10);
		else
			_gtime_interval = sign * _sampling*10;
	}
	else {
		if (_sampling > 1)
			_gtime_interval=int(sign * _sampling);
		else
			_gtime_interval = sign * _sampling;
	}
	//_gtime_interval = 20;//?????????????????????????




	//调试信息
	string _output_path = ".";
	string o_path = _output_path + "/GREAT_FGO_FLOAT.txt";
	_output_float_solution.open(o_path, ofstream::out);
	if (_output_float_solution.good()) {
		_output_float_solution << "#FACTOR GRAPH OPTIMIZATION BASED GNSS SOLUTION" << endl;
		//_output_float_solution << "#" << "ambiguity propogation: " << dynamic_cast<t_gsetfgo*>(gset)->_amb_propagation() << endl;
		_output_float_solution << "#" << "FGO WINDOW LENGTH: " << dynamic_cast<t_gsetfgo*>(gset)->gwindow_size() << endl;
		set<string> sys = dynamic_cast<t_gsetgen*>(gset)->sys();
		_output_float_solution << "#" << "GNSS system: ";
		for (auto it = sys.begin(); it != sys.end(); it++)
		{
			_output_float_solution << *it << " ";
		}
		_output_float_solution << endl;
		if (_isBase)
		{
			t_gtriple xyz_base = _gallobj->obj(_site_base)->crd_arp(_epoch);
			_output_float_solution << "#Base coordinate: " << fixed << setprecision(4) << xyz_base[0] << "," << xyz_base[1] << "," << xyz_base[2] << endl;
			_output_float_solution << endl;
		}
		_output_float_solution << "#" << setw(15) << "Seconds of Week"
			<< setw(12) << "X-ECEF " << // [m]
			" " << setw(15) << "Y-ECEF" <<      // [m]
			" " << setw(15) << "Z-ECEF" <<      // [m]
			" " << setw(10) << "Vx-ECEF" <<      // [m/s]
			" " << setw(10) << "Vy-ECEF" <<      // [m/s]
			" " << setw(10) << "Vz-ECEF" <<      // [m/s]
			" " << setw(9) << "X-RMS"
			<< " " << setw(9) << "Y-RMS"
			<< " " << setw(9) << "Z-RMS"
			<< " " << setw(9) << "Vx-RMS"
			<< " " << setw(9) << "Vy-RMS"
			<< " " << setw(9) << "Vz-RMS" <<
			" " << setw(5) << "NSat"
			<< " " << setw(5) << "PDOP"
			<< " " << setw(8) << "sigma0"
			<< " " << setw(10) << "AmbStatus";

		_output_float_solution << endl;
		_output_float_solution << "#" << setw(15) << "(s)" <<
			" " << setw(12) << "(m)" << // [m]
			" " << setw(15) << "(m)" <<      // [m]
			" " << setw(15) << "(m)" <<      // [m]
			" " << setw(10) << "(m/s)" << // [m/s]
			" " << setw(10) << "(m/s)" << // [m/s]
			" " << setw(10) << "(m/s)" << // [m/s]
			" " << setw(9) << "(m)" <<      // [m]
			" " << setw(9) << "(m)" <<      // [m]
			" " << setw(9) << "(m)" <<
			" " << setw(9) << "(m/s)" <<  // [m/s]
			" " << setw(9) << "(m/s)" <<  // [m/s]
			" " << setw(9) << "(m/s)" <<      // [m/s]
			" " << setw(5) << "(#)" <<
			" " << setw(5) << "(#)" <<
			" " << setw(8) << "(m)" <<
			" " << setw(10) << " ";

		_output_float_solution << endl;
	}
}

gfgomsf::t_gpvtfgo::~t_gpvtfgo()
{
	 delete _last_gnss_info;
	 delete _last_gnss_marginalization_info;
	 delete _gbias_model;
	 for (int i = 0; i<=_rover_count; i++)
	 {
		 if (_rover_window[i] != nullptr)
		 {
			 delete _rover_window[i];
		 }

	 }

}

int gfgomsf::t_gpvtfgo::processBatch(const t_gtime &beg_r, const t_gtime &end_r, bool prtOut)
{
#ifdef BMUTEX
	boost::mutex::scoped_lock lock(_mutex);
#endif
	_gmutex.lock();

	if (_grec == nullptr)
	{
		ostringstream os;
		os << "ERROR: No object found (" << _site << "). Processing terminated!!! " << beg_r.str_ymdhms() << " -> " << end_r.str_ymdhms() << endl;
		if (_spdlog) SPDLOG_LOGGER_ERROR(_spdlog, string("t_gpvtfgo ") + os.str());
		_gmutex.unlock();
		return -1;
	}

	int sign = 1;


	double subint = 0.1;

	if (!_beg_end)
		sign = -1;

	InitProc(beg_r, end_r, &subint);

	t_gtime now(_beg_time);

	std::cerr << _site << ": Start GNSS Processing: " << now.str_ymdhms() << " " << _end_time.str_ymdhms() << endl;

	bool time_loop = true;
	while (time_loop)
	{

		if (_beg_end && (now < _end_time || now == _end_time))
			time_loop = true;
		else if (_beg_end && now > _end_time)
		{
			time_loop = false;
			break;
		}
		// synchronization
		if (now != _end_time)
		{
			if (!time_sync(now, _sampling, _scale, _spdlog))
			{									   // now.sod()%(int)_sampling != 0 ){
				now.add_dsec(sign * subint / 100); // add_dsec used for synchronization!

				continue;
			}
			if (_sampling > 1)
				now.reset_dsec();
		}

		bool irc_slip = _slip_detect(now);


		int irc_epo = 0;
		irc_epo = processWindow(now); //sliding window model
		if (irc_epo < 0)
		{
			if (_sampling > 1)
				now.add_secs(int(sign * _sampling)); // =<1Hz data
			else
				now.add_dsec(sign * _sampling); //  >1Hz data

			continue;
		}
		else
			_success = true;

		if (_spdlog) SPDLOG_LOGGER_ERROR(_spdlog, string("t_gpvtfgo ") ,_site + now.str_ymdhms(" processing epoch: "));
		double percent = now.diff(_beg_time) / _end_time.diff(_beg_time) * 100.0;
		std::cerr << "\r" << now.str_ymdhms() << setw(5) << " Q = " << (_amb_state ? 1 : 2) << fixed << setprecision(1) << setw(6) << percent << "%";

		if (_sampling > 1)
			now.add_secs(int(sign * _sampling)); // =<1Hz data
		else
			now.add_dsec(sign * _sampling); //  >1Hz data




	}
	_gmutex.unlock();
	return 1;
}

int gfgomsf::t_gpvtfgo::processWindow(const t_gtime & now, vector<t_gsatdata>* data_rover, vector<t_gsatdata>* data_base)
{
	if (!_get_gdata(now, data_rover, data_base))
		return -1;
	if (!_isBase)
		_record_ppp_cycle_slips(_data);
	t_gtime runEpoch = _data.begin()->epoch();
	_timeUpdate(runEpoch);
	_epoch = runEpoch;


	if (_reset_par > 0)
	{
		if (now.sod() % _reset_par == 0)
		{
			_reset_param();
		}
	}
	// save apriory coordinates
	if (_crd_est != CONSTRPAR::FIX)
		_saveApr(runEpoch, _param, _Qx);
	if (!_crd_xml_valid())
		_sig_init_crd = 100.0;
	// select obs for tb log (rover)
	// add
	if (!_isBase)
			_cntrep = 0; // number of iterations caused by outliers

	if (_prepareData() < 0) return -1;
	if (_isBase)
	{
		_set_rec_info(_gallobj->obj(_site_base)->crd_arp(_epoch), _vBanc(4), _vBanc_base(4));
	}
	// _prepareData() filters its internal GNSS information, but the PPP/FGO
	// factor path below still consumes _data.  Remove observations below the
	// configured elevation mask here as well, before ambiguity creation and
	// factor construction.  The small guard band prevents the nominal filter
	// elevation and the factor's recomputed elevation from landing on opposite
	// sides of the hard cutoff because of state/rounding differences.
	if (!_isBase)
	{
		for (auto it = _data.begin(); it != _data.end();)
		{
			if (it->ele_deg() < _minElev + 0.1)
				it = _data.erase(it);
			else
				++it;
		}
	}
	if (_data.size() < _minsat)
	{
		if (_spdlog) SPDLOG_LOGGER_ERROR(_spdlog, string("gpvtfgo "), ("Not enough visible satellites!"));
		return -1;
	}
	_get_initial_value(runEpoch); //for current epoch
	// ??????  necessary
	if (_data.size() < 6)
	{
		cout << "nsat < 6!!" << endl;
	}

	// add hwzhang
	//RTK FGO
	if (_isBase && _combine_DD() < 0)
	{
		std::cout << "Epoch: " << runEpoch.sow() << " combine DD wrong!!!" << endl;
		clearWindow();
		return -1;
	}
	if (!_isBase && _combine_IF() < 0)
	{
		std::cout << "Epoch: " << runEpoch.sow() << " combine IF wrong!!!" << endl;
		_rollback_current_ppp_node();
		return -1;
	}

	if (_isBase) {
		_prepare_equ();
	}


	if (_isBase)
	{
		_optimization();
	}
	else
	{
		if (_optimization_PPP() < 0)
		{
			std::cout << "Epoch: " << runEpoch.sow()
				<< " rejected: outlier remains at minimum satellite count ("
				<< _minsat << ")" << endl;
			_rollback_current_ppp_node();
			return -1;
		}
	}


	// ambiguity resolution


	if (_last_gnss_info->valid)
	{
		publish_foat();
		//if(_isBase) {
			_pre_amb_resolution();
			_amb_resolution();
		//}
	}
	else
	{
		std::cout << "Epoch: " << runEpoch.sow() << "solving faild" << endl;
		if (_isBase)
			clearWindow();
		else
			_rollback_current_ppp_node();
		return -1;
	}

	if (_isBase)
	{
		_marginalization();
	}
	else
	{
		_marginalization_PPP();//modify internally for MultiWindow
	}
	if (!_isBase)
		_commit_ppp_cycle_slips();
	//if (_rover_count == 1)
	//	_initial_prior = false;

	_slide_window();
	return _amb_state ? 1 : 0;
}


void gfgomsf::t_gpvtfgo::clearWindow()
{
    for (int i = 0; i <= _rover_count; i++)
    {
    	_headers[i]=0.0;
    	_Pos[i].setZero();
    	if (!_isBase)
    	{
    		_clk[i] = 0.0;
    		_trp[i] = 0.0;
    		_isb_GAL[i] = 0.0;
    		_isb_BDS[i] = 0.0;
			_isb_GLO[i] = 0.0;
    		_lost_isb_GAL[i] = false;
    		_lost_isb_BDS[i] = false;
			_lost_isb_GLO[i] = false;
    	}
		if (_rover_window[i] != nullptr)
			delete _rover_window[i];
		_rover_window[i] = nullptr;
		_para_window[i].delAllParam();
    }
	if (_last_gnss_info != nullptr)
		delete _last_gnss_info;
	_last_gnss_info = nullptr;
	if (_last_gnss_marginalization_info != nullptr)
		delete _last_gnss_marginalization_info;
	_last_gnss_marginalization_info = nullptr;
	if (_isBase) {
		_win_base_data.clear();
		_win_base_data.resize(gwindow_size + 1);
	}
	cur_sat_prn.clear();
	if (_isBase){
	_DD_msg.clear();
	_vDD_msg.clear();
		}
	if (!_isBase)
	{
		_IF_msg.clear();
		_vIF_msg.clear();
	}
	if (_isBase)
	{
		_amb_manager->clearState();
	}
	else
	{
		_ambIF_manager->clearState();
		_pending_ppp_slips.clear();
		_candidate_ppp_slips.clear();
		_last_ppp_phase_epoch.clear();
	}
	_rover_count = -1;
    _global_sat_id = -1;
	_global_amb_id = -1;
	_ppp_candidate_global_sat_id = -1;
	_ppp_candidate_global_amb_id = -1;
	_initial_prior = true;
	_amb_state = false;
	memset(_para_amb, 0, sizeof(_para_amb));
	memset(_para_CRD, 0, sizeof(_para_CRD));
	if (!_isBase)
	{
		memset(_para_CLK, 0, sizeof(_para_CLK));
		memset(_para_TRP, 0, sizeof(_para_TRP));
		memset(_para_ISB_GAL, 0, sizeof(_para_ISB_GAL));
		memset(_para_ISB_BDS, 0, sizeof(_para_ISB_BDS));
		memset(_para_ISB_GLO, 0, sizeof(_para_ISB_GLO));
		memset(_para_AMB_IF, 0, sizeof(_para_AMB_IF));
	}
}


void gfgomsf::t_gpvtfgo::publish_foat()
{
    // get CRD params
    t_gtriple xyz, ell;
    _all_para_win.getCrdParam(_site, xyz, _epoch, _epoch);
    xyz2ell(xyz, ell, false);

    // CRD using eccentricities
    t_gtriple xyz_ecc = xyz - _grec->eccxyz(_epoch); // MARKER + ECC = ARP

    double Xrms = 0.0, Yrms = 0.0, Zrms = 0.0,
        Vxrms = 0.0, Vyrms = 0.0, Vzrms = 0.0;
    Xrms = sqrt(_last_gnss_info->Qx(0, 0));
    Yrms = sqrt(_last_gnss_info->Qx(1, 1));
    Zrms = sqrt(_last_gnss_info->Qx(2, 2));

    t_gtriple crd_rms(Xrms, Yrms, Zrms);
    t_gtriple vRec(0, 0, 0);
    vRec = t_gtriple(0.0, 0.0, 0.0);
    double pdop = sqrt(_last_gnss_info->Qx(0, 0) + _last_gnss_info->Qx(1, 1) + _last_gnss_info->Qx(2, 2));

    set<string> ambs = _all_para_win.amb_prns();
    int nsat = ambs.size();

    // get amb status
    string amb = "Float";

    t_gtriple blh;
    xyz2ell(xyz_ecc, blh, true);

    string str_dsec = dbl2str(_epoch.dsec());

    double bl = 0;
    if (_isBase)
    {
        auto crd_base = _gallobj->obj(_site_base)->crd_arp(_epoch);
        t_gtriple tmpell, tmpdxyz, tmpneu;
        xyz2ell(crd_base, tmpell, false);
        tmpdxyz = xyz - crd_base;
        xyz2neu(tmpell, tmpdxyz, tmpneu);
        bl = tmpneu.norm();
    }

    _output_float_solution << fixed << setprecision(4) << " "
        // << epoch.str_ymdhms() << str_dsec.substr(2) << setprecision(4)
        << " " << _epoch.sow() + _epoch.dsec()
        << fixed << setprecision(4)
        << " " << setw(10) << xyz_ecc[0] // [m]
        << " " << setw(10) << xyz_ecc[1] // [m]
        << " " << setw(10) << xyz_ecc[2] // [m]
        << " " << setw(10) << vRec[0]    // [m/s]
        << " " << setw(10) << vRec[1]    // [m/s]
        << " " << setw(10) << vRec[2]    // [m/s]
        << " " << setw(10) << Xrms       // [m]
        << " " << setw(10) << Yrms       // [m]
        << " " << setw(10) << Zrms       // [m]
        << " " << setw(10) << Vxrms      // [m/s]
        << " " << setw(10) << Vyrms      // [m/s]
        << " " << setw(10) << Vzrms      // [m/s]
        << fixed << setprecision(0)
        << " " << setw(5) << nsat // nsat
        << fixed << setprecision(2)
        << " " << setw(5) << pdop // pdop
        << fixed << setprecision(2)
        << " " << setw(5) << _last_gnss_info->sig_unit // pdop
        << fixed << setprecision(2)
        << " " << setw(8) << amb
        << endl;
}

// only use for RTK
void gfgomsf::t_gpvtfgo::_prepare_equ()
{
	if (!_DD_msg.empty())
	{
		for (int i = 0; i < _vDD_msg.size(); i++)
		{
			for (int j = 0; j < _vDD_msg[i].size(); j++)
			{
				_vDD_msg[i][j].base_nonref_sat.is_process = false;
				_vDD_msg[i][j].base_ref_sat.is_process = false;
				_vDD_msg[i][j].rover_nonref_sat.is_process = false;
				_vDD_msg[i][j].rover_ref_sat.is_process = false;
			}
		}
		_update_all_equ();
	}
}
// for RTK
bool gfgomsf::t_gpvtfgo::_update_all_equ()
{
	t_gprecisebiasFGO gbias_model = *_gbias_model;
	for (int i = 0; i <= _rover_count; i++)//epoch loop
	{
		//sat-rec-band-all_obstype
		map<string, map<string, map<GOBSBAND, set<GOBSTYPE>>>> obs_info_map;
		t_gallpar params_temp(_para_window[i]);
		//scan all obs
		for (auto &dd_iter : _vDD_msg[i])
		{
			string sat_non = dd_iter.base_nonref_sat.sat();
			string sat_ref = dd_iter.base_ref_sat.sat();
			string base = dd_iter.base_site;
			string rover = dd_iter.rover_site;
			GOBSBAND band = dd_iter.band;
			GOBSTYPE obs_type = dd_iter.obs_type;

			// 1 DDEquMsg => 4 rec-sat pair
			obs_info_map[sat_non][base][band].insert(obs_type);
			obs_info_map[sat_non][rover][band].insert(obs_type);
			obs_info_map[sat_ref][base][band].insert(obs_type);
			obs_info_map[sat_ref][rover][band].insert(obs_type);
		}

		//cmb all equ, only want: _prepare_obs_GPP
		for (auto i_sat = obs_info_map.begin(); i_sat != obs_info_map.end(); i_sat++)//sat loop
		{
			string sat = i_sat->first;
			for (auto i_rec = obs_info_map[sat].begin(); i_rec != obs_info_map[sat].end(); i_rec++)//rec loop
			{
				string rec = i_rec->first;
				t_gsatdata satdata;
				t_gbaseEquation tempP, tempL;
				//get gsatdata according to sat-rec pair
				if (!_get_satdata(i, sat, rec, satdata))
				{
					continue;
				}
				for (auto i_band = obs_info_map[sat][rec].begin(); i_band != obs_info_map[sat][rec].end(); i_band++)//band loop
				{
					GOBSBAND band = i_band->first;
					t_gobs obsP(satdata.select_range(band));
					t_gobs obsL(satdata.select_phase(band));
					t_gtime crt = satdata.epoch();
					if (!gbias_model.cmb_equ(true, false, crt, params_temp, satdata, obsP, tempP))
					{
						continue;
					}
					if (!gbias_model.cmb_equ(true, false, crt, params_temp, satdata, obsL, tempL))
					{
						continue;
					}
					_set_satdata(i, sat, rec, satdata);
				}//end band loop
			}//end rec loop
		}//end sat loop

	}//end epoch loop

	return true;
}
// for RTK
bool gfgomsf::t_gpvtfgo::_get_satdata(const int epoch_id, const string & sat, const string & rec, t_gsatdata & satdata)
{
	if (rec == _vDD_msg[epoch_id][0].base_site)//base rec
	{
		for (int i = 0; i < _vDD_msg[epoch_id].size(); i++)
		{
			if (_vDD_msg[epoch_id][i].base_nonref_sat.sat() == sat)
			{
				satdata = _vDD_msg[epoch_id][i].base_nonref_sat;
				return true;
			}
			if (_vDD_msg[epoch_id][i].base_ref_sat.sat() == sat)
			{
				satdata = _vDD_msg[epoch_id][i].base_ref_sat;
				return true;
			}
		}
	}
	else if (rec == _vDD_msg[epoch_id][0].rover_site)//rover rec
	{
		for (int i = 0; i < _vDD_msg[epoch_id].size(); i++)
		{
			if (_vDD_msg[epoch_id][i].rover_nonref_sat.sat() == sat)
			{
				satdata = _vDD_msg[epoch_id][i].rover_nonref_sat;
				return true;
			}
			if (_vDD_msg[epoch_id][i].rover_ref_sat.sat() == sat)
			{
				satdata = _vDD_msg[epoch_id][i].rover_ref_sat;
				return true;
			}
		}
	}
	else//no rec
	{
		return false;
	}
	return false;
}

bool gfgomsf::t_gpvtfgo::_set_satdata(const int epoch_id, const string & sat, const string & rec, const t_gsatdata & satdata)
{
	bool updated = false;
	for (int i = 0; i < _vDD_msg[epoch_id].size(); i++)
	{
		// 1 DDEquMsg => 4 rec-sat pair
		pair<string, string> rec_sat = make_pair(rec, sat);
		pair<string, string> rec_sat1 = make_pair(_vDD_msg[epoch_id][i].base_site, _vDD_msg[epoch_id][i].base_nonref_sat.sat());
		pair<string, string> rec_sat2 = make_pair(_vDD_msg[epoch_id][i].base_site, _vDD_msg[epoch_id][i].base_ref_sat.sat());
		pair<string, string> rec_sat3 = make_pair(_vDD_msg[epoch_id][i].rover_site, _vDD_msg[epoch_id][i].rover_nonref_sat.sat());
		pair<string, string> rec_sat4 = make_pair(_vDD_msg[epoch_id][i].rover_site, _vDD_msg[epoch_id][i].rover_ref_sat.sat());
		if (rec_sat == rec_sat1)
		{
			_vDD_msg[epoch_id][i].base_nonref_sat = satdata;
			updated = true;
		}
		else if (rec_sat == rec_sat2)
		{
			_vDD_msg[epoch_id][i].base_ref_sat = satdata;
			updated = true;
		}
		else if (rec_sat == rec_sat3)
		{
			_vDD_msg[epoch_id][i].rover_nonref_sat = satdata;
			updated = true;
		}
		else if (rec_sat == rec_sat4)
		{
			_vDD_msg[epoch_id][i].rover_ref_sat = satdata;
			updated = true;
		}
	}
	return updated;
}

void gfgomsf::t_gpvtfgo::_set_glonass_channels(vector<t_gsatdata>& data, const t_gtime& epoch)
{
	map<string, int> channels;
	if (_gobs)
		channels = _gobs->glo_freq_num();
	if (_gnav)
	{
		for (const auto& item : _gnav->glo_freq_num())
			if (channels.find(item.first) == channels.end())
				channels[item.first] = item.second;
	}

	for (auto it = data.begin(); it != data.end();)
	{
		if (it->gsys() != GSYS::GLO)
		{
			++it;
			continue;
		}

		int channel = it->channel();
		if (channel < -7 || channel > 13)
		{
			auto found = channels.find(it->sat());
			if (found != channels.end())
				channel = found->second;
		}

		if (channel < -7 || channel > 13)
		{
			if (_spdlog)
				_spdlog->warn("PPP FGO: remove GLONASS satellite {} at {} because its FDMA channel is unavailable",
					it->sat(), epoch.str_ymdhms(""));
			it = data.erase(it);
			continue;
		}

		it->channel(channel);
		++it;
	}
}

bool gfgomsf::t_gpvtfgo::_get_gdata(const t_gtime& now, vector<t_gsatdata>* data_rover, vector<t_gsatdata>* data_base)
{
	if (/*_data.size()*/ _getData(now, data_rover, false) == 0)
	{
		if (_spdlog)
		{
			SPDLOG_LOGGER_DEBUG(_spdlog, _site + now.str_ymdhms(" no observation found at epoch: "));
		}
		return false;
	}
	// apply dcb
	if (_gallbias)
	{
		for (auto& itdata : _data)
		{
			itdata.apply_bias(_gallbias);
		}
	}
	if (!_isBase)
		_set_glonass_channels(_data, now);

	vector<t_gsatdata>::iterator it = _data.begin();
	string double_freq = "";
	string single_freq = "";


	_sat_freqs.clear();
	while (it != _data.end())
	{
		GOBSBAND b1 = _band_index[it->gsys()][FREQ_1];
		GOBSBAND b2 = _band_index[it->gsys()][FREQ_2];

		auto obsL1 = it->select_phase(b1);
		auto obsL2 = it->select_phase(b2);

		if (obsL1 == GOBS::X && obsL2 != GOBS::X || obsL1 != GOBS::X && obsL2 == GOBS::X)
		{
			single_freq += "  " + it->sat();
			_sat_freqs[it->sat()] = "1";
		}

		if (obsL1 != GOBS::X && obsL2 != GOBS::X)
		{
			double_freq += "  " + it->sat();
			_sat_freqs[it->sat()] = "2";
		}

		++it;
	}


	if (_isBase)
	{

		if (/*_data_base.size()*/_getData(now, data_base, true) == 0)
		{
			if (_spdlog)
			{
				SPDLOG_LOGGER_DEBUG(_spdlog, _site_base + now.str_ymdhms(" no observation found at epoch: "));
			}
			return false;
		}
		// apply dcb
		if (_gallbias)
		{
			for (auto& itdata_base : _data_base)
			{
				itdata_base.apply_bias(_gallbias);
			}
		}
	}


	if (_gallobj != nullptr) {
		auto it_data = _data.begin();
		while (it_data != _data.end()) {
			string sat_id = it_data->sat();
			shared_ptr<t_gobj> sat_obj = _gallobj->obj(sat_id);


			if (sat_obj == nullptr) {
				if (_spdlog) {
					SPDLOG_LOGGER_DEBUG(_spdlog, "remove satellite " + sat_id + " due to missing object");
				}
				it_data = _data.erase(it_data);
			}
			else {
				shared_ptr<t_gpcv> sat_pcv = sat_obj->pcv(now);
				if (sat_pcv == nullptr) {
					if (_spdlog) {
						SPDLOG_LOGGER_DEBUG(_spdlog, "remove satellite " + sat_id + " due to missing PCV data");
					}
					it_data = _data.erase(it_data);
				}
				else {
					++it_data;
				}
			}
		}


		if (_isBase && data_base != nullptr) {
			auto it_base = data_base->begin();
			while (it_base != data_base->end()) {
				string sat_id = it_base->sat();
				shared_ptr<t_gobj> sat_obj = _gallobj->obj(sat_id);

				if (sat_obj == nullptr) {
					it_base = data_base->erase(it_base);
				}
				else {
					shared_ptr<t_gpcv> sat_pcv = sat_obj->pcv(now);
					if (sat_pcv == nullptr) {
						it_base = data_base->erase(it_base);
					}
					else {
						++it_base;
					}
				}
			}
		}
	}


	return true;
}


void gfgomsf::t_gpvtfgo::_get_initial_value(const t_gtime& runEpoch)
{


	// Predict ambiguity
	// Add or remove ambiguity parameter and appropriate rows/columns covar. matrix
	if (!_isBase)
	{
		_cntrep++;
		_syncSys();
	}
	if (_phase)
	{
		_syncAmb();
	}
	_predictCrd();
	if (!_isBase)
	{
		_predictClk();
		_predictBias();
		_predictTropo();
	}
	_predictAmb();
	//_set_initial_value(runEpoch);//delete for MultiWindow
	_set_initial_value(runEpoch);//add for MultiWindow
	_initialized = true;
}


void gfgomsf::t_gpvtfgo::_set_initial_value(const t_gtime& runEpoch)
{
    _rover_count++;
    if (_rover_count >= 1)
    {
        //detect time gap
        if (runEpoch.sow() - _headers[_rover_count - 1] > _gtime_interval)
        {
            cout << "[" << _headers[_rover_count - 1] << "]--" << "[" << runEpoch.sow() << "]" << ": time gap !!!]" << endl;
            clearWindow();
            _rover_count++;
        }
    }
    t_gtriple xyz;
    _param.getCrdParam(_site, xyz);  //!!!need to check the solution!!!
    _Pos[_rover_count] = Eigen::Vector3d(xyz.crd(0), xyz.crd(1), xyz.crd(2));
    if (!_isBase)
    {
        int id = 0;
        id = _param.getParam(_site, par_type::CLK, "");
        if (id >= 0)
        {
            _clk[_rover_count] = _param[id].value();
        }

        if (_rover_count > 0)
        {
            _trp[_rover_count] = _trp[_rover_count - 1];
        }
        else
        {
            id = _param.getParam(_site, par_type::TRP, "");
            if (id >= 0)
            {
                _trp[_rover_count] = _param[id].value();
            }
        }

        id = _param.getParam(_site, par_type::GAL_ISB, "");
        if (id >= 0)
        {
            _lost_isb_GAL[_rover_count] = false;
            if (_rover_count > 0)
            {
                _isb_GAL[_rover_count] = _isb_GAL[_rover_count - 1];
            }
            else
            {
                _isb_GAL[_rover_count] = _param[id].value();
            }
        }
        else
        {
            _lost_isb_GAL[_rover_count] = true;
            _isb_GAL[_rover_count] = 0.0;
        }

        id = _param.getParam(_site, par_type::BDS_ISB, "");
        if (id >= 0)
        {
            _lost_isb_BDS[_rover_count] = false;
            if (_rover_count > 0)
            {
                _isb_BDS[_rover_count] = _isb_BDS[_rover_count - 1];
            }
            else
            {
                _isb_BDS[_rover_count] = _param[id].value();
            }
        }
        else
        {
            _lost_isb_BDS[_rover_count] = true;
            _isb_BDS[_rover_count] = 0.0;
        }

		id = _param.getParam(_site, par_type::GLO_ISB, "");
		if (id >= 0)
		{
			_lost_isb_GLO[_rover_count] = false;
			if (_rover_count > 0)
				_isb_GLO[_rover_count] = _isb_GLO[_rover_count - 1];
			else
				_isb_GLO[_rover_count] = _param[id].value();
		}
		else
		{
			_lost_isb_GLO[_rover_count] = true;
			_isb_GLO[_rover_count] = 0.0;
		}

        //save the ini value
        double epsilon = 1e-10;
        if (fabs(_trp_ini) < epsilon)
        {
            _trp_ini = _trp[_rover_count];
        }
        if (fabs(_isb_GAL_ini) < epsilon && !_lost_isb_GAL[_rover_count])
        {
            _isb_GAL_ini = _isb_GAL[_rover_count];
        }
        if (fabs(_isb_BDS_ini) < epsilon && !_lost_isb_BDS[_rover_count])
        {
            _isb_BDS_ini = _isb_BDS[_rover_count];
        }
		if (fabs(_isb_GLO_ini) < epsilon && !_lost_isb_GLO[_rover_count])
		{
			_isb_GLO_ini = _isb_GLO[_rover_count];
		}
    }
    _headers[_rover_count] = runEpoch.sow();
	// new add
	_rover_window[_rover_count] =new t_grover_msg(runEpoch, _Pos[_rover_count]);
    if (_isBase)
    {	// new add
    	_win_base_data[_rover_count] = _data_base;
        t_gallpar params_add;
        _gtemp_params(_param, params_add);
        _para_window[_rover_count] = params_add;
    }
    else
    {
        _para_window[_rover_count] = _param;
		// Snapshot the ID generators before this candidate epoch creates any
		// satellite/ambiguity arcs.  A rejected epoch must not consume IDs.
		_ppp_candidate_global_sat_id = _global_sat_id;
		_ppp_candidate_global_amb_id = _global_amb_id;
		_candidate_ppp_slips.clear();
    }
    if (_isBase)
    {
        for (auto it : _data)
        {
            string sat_name = it.sat();
            if (!_amb_manager->is_sat_tracking(sat_name))
            {
                //new sat
                _global_sat_id++;
                _amb_manager->addNewSat(runEpoch, _rover_count, _global_sat_id, _global_amb_id, it, _param);
            }
            else
            {
                bool is_slip = false;
                vector<GOBSBAND> band = dynamic_cast<t_gsetgnss*>(_set)->band(it.gsys());
                GSYS gs = it.gsys();
                int nf = 5;
                if (band.size())
                    nf = band.size();
                for (FREQ_SEQ f = FREQ_1; f <= nf; f = (FREQ_SEQ)(f + 1))
                {
                    GOBSBAND b;
                    if (f > _frequency)
                        continue; // modified by lvhb in 20201211
                    if (band.size() >= f)
                    { // automatic dual band selection -> for Anubis purpose
                        b = band[f - 1];
                    }
                    else
                        b = t_gsys::band_priority(gs, f);

                    t_gobs gobs(it.select_phase(b));
                    if (true && it.getlli(gobs.gobs()) >= 1)
                    {
                        // slip occured
                        _global_sat_id++;
                        _amb_manager->addNewSat(runEpoch, _rover_count, _global_sat_id, _global_amb_id, it, _param);
                        is_slip = true;
                        break;
                    }
                }
                if (!is_slip)
                {
                    //_amb_manager->addRover(sat_name, _rover_count);
                    _amb_manager->addRover(_epoch.sow(), sat_name, _rover_count);
                }
            }
        }
        _amb_manager->get_last_epoch_sats(_epoch.sow());
        assert(_amb_manager->cur_sats.size() == _data.size());
    }
    else
    {
		for (const auto& sat_data : _data)
		{
			const string sat_name = sat_data.sat();
			const bool was_tracking =
				_ambIF_manager->is_sat_tracking(sat_name);
			const bool force_new_arc =
				_pending_ppp_slips.find(sat_name) !=
				_pending_ppp_slips.end();

			if (!was_tracking || force_new_arc)
			{
				_global_sat_id++;
				_ambIF_manager->addNewSat(
					runEpoch, _rover_count, _global_sat_id,
					_global_amb_id, sat_data, _param);

				if (force_new_arc)
					_candidate_ppp_slips.insert(sat_name);
			}
			else
			{
				_ambIF_manager->addRover(
					_epoch.sow(), sat_name, _rover_count);
			}
		}
        _ambIF_manager->get_last_epoch_sats(_epoch.sow());
        assert(_ambIF_manager->cur_sats.size() == _data.size());
    }
}
bool gfgomsf::t_gpvtfgo::_gtemp_params(t_gallpar & params, t_gallpar & params_temp)
{
	params_temp = params;
	t_gpar par_x_base; par_x_base.site = _site_base; par_x_base.parType = par_type::CRD_X; par_x_base.value(_gcrd_base[0]);
	t_gpar par_y_base; par_y_base.site = _site_base; par_y_base.parType = par_type::CRD_Y; par_y_base.value(_gcrd_base[1]);
	t_gpar par_z_base; par_z_base.site = _site_base; par_z_base.parType = par_type::CRD_Z; par_z_base.value(_gcrd_base[2]);
	t_gpar par_clk_rover; par_clk_rover.site = _site; par_clk_rover.parType = par_type::CLK; par_clk_rover.value(_gclk_rover);
	t_gpar par_clk_base; par_clk_base.site = _site_base; par_clk_base.parType = par_type::CLK; par_clk_base.value(_gclk_base);
	params_temp.addParam(par_x_base);
	params_temp.addParam(par_y_base);
	params_temp.addParam(par_z_base);
	params_temp.addParam(par_clk_rover);
	params_temp.addParam(par_clk_base);
	params_temp.reIndex();
	return true;
}

void gfgomsf::t_gpvtfgo::_set_rec_info(const t_gtriple & xyz_base, double clk_rover, double clk_base)
{
	_gcrd_base = xyz_base;
	_gclk_rover = clk_rover;
	_gclk_base = clk_base;
}

int gfgomsf::t_gpvtfgo::_combine_DD()
{	 if (!_isBase)
	{
	return -1;
	}
	int flag = -1;
	_select_ref_sat();
	auto it_dd = _DD_msg.begin();
	while (it_dd != _DD_msg.end())
	{
		if (!_get_DD_data(*it_dd, _data_base))
		{
			cout << "[" << it_dd->rover_ref_sat.sat() << "-" << it_dd->rover_nonref_sat.sat() << it_dd->obs_type << "," << it_dd->freq << "]" << "DD msg wrong!" << endl;
			it_dd = _DD_msg.erase(it_dd);
			continue;
		}

		it_dd++;
	}
	if (!_DD_msg.empty())
	{
		_vDD_msg.push_back(_DD_msg);
		flag = 1;
	}
	return flag;
}

int gfgomsf::t_gpvtfgo::_combine_IF()
{
    if (_isBase)
    {
        return -1;
    }
    _IF_msg.clear();
    vector<t_gsatdata>::iterator it = _data.begin();
    for (it = _data.begin(); it != _data.end();)
    {
        const GOBSBAND& b1 = _band_index[it->gsys()][FREQ_1];
        const GOBSBAND& b2 = _band_index[it->gsys()][FREQ_2];

        it->apply_bias(_gallbias);

        if (b1 == BAND || b2 == BAND)
        {	it++;
            continue;
        }
        t_gobs gobsP1(it->select_range(b1));
        t_gobs gobsP2(it->select_range(b2));
        t_gobs gobsL1(it->select_phase(b1));
        t_gobs gobsL2(it->select_phase(b2));

        double P1_value = it->getobs(gobsP1.gobs());
        double P2_value = it->getobs(gobsP2.gobs());
        double L1_value = it->getobs(gobsL1.gobs());
        double L2_value = it->getobs(gobsL2.gobs());
        auto gsys = it->gsys();

        if (!double_eq(P1_value, 0.0) && !double_eq(P2_value, 0.0))
        {
            IFEquMsg if_msg(*it, GOBSTYPE::TYPE_C, _freq_index[gsys][b1], _freq_index[gsys][b2], b1, b2);
            if_msg.sat_id = _ambIF_manager->get_sat_id(it->sat());//add for MultiWindow
            if_msg.time = it->epoch();
            if_msg.site = _site;
            _IF_msg.push_back(if_msg);
        }
        if (!double_eq(L1_value, 0.0) && !double_eq(L2_value, 0.0))
        {
            IFEquMsg if_msg(*it, GOBSTYPE::TYPE_L, _freq_index[gsys][b1], _freq_index[gsys][b2], b1, b2);
            //if_msg.sat_id = cur_sat_prn[it->sat()];//delete for MultiWindow
            if_msg.sat_id = _ambIF_manager->get_sat_id(it->sat());//add for MultiWindow
            if_msg.time = it->epoch();
            if_msg.site = _site;
            _IF_msg.push_back(if_msg);
        }
        _combineMW(*it);//for AMB FIX

        it++;
    }

    if (!_IF_msg.empty())
    {
        _vIF_msg.push_back(_IF_msg);
        return 1;
    }
    else
    {
        return -1;
    }
}

void gfgomsf::t_gpvtfgo::_prior_factor(ceres::Problem & problem)
{
	//prior
	if (_last_gnss_marginalization_info && _last_gnss_marginalization_info->valid)
	{
		MarginalizationGNSSFactor *marginalization_gnss_factor = new MarginalizationGNSSFactor(_last_gnss_marginalization_info);
		problem.AddResidualBlock(marginalization_gnss_factor, NULL,
			_last_gnss_marginalization_para_blocks);
	}
}

void gfgomsf::t_gpvtfgo::_double_to_vector()
{	// new add yfj
	all_parameter_block.clear();
	for (int i = 0; i <= _rover_count; i++)
	{
		_Pos[i] = Eigen::Vector3d(_para_CRD[i][0], _para_CRD[i][1], _para_CRD[i][2]);

	}
	if (_isBase) {
		for (int i = 0; i < _amb_manager->ambiguity_ids.size(); i++)
		{
			int amb_id = _amb_manager->ambiguity_ids[i];
			_amb_manager->updateAmb(amb_id, _para_amb[amb_id][0]);

			//cout << amb_id <<": "<< _para_amb[amb_id][0] << "   " ;
		}
	}
	else
	{
		for (int i = 0; i < _ambIF_manager->ambiguity_ids.size(); i++)
		{
			int amb_id = _ambIF_manager->ambiguity_ids[i];
			_ambIF_manager->updateAmb(amb_id, _para_AMB_IF[amb_id][0]);
			//cout << amb_id <<": "<< _para_amb[amb_id][0] << "   " ;
		}
	}
	// PPP
	if (!_isBase)
	{
		for (int i = 0; i <= _rover_count; i++)
		{
			_clk[i] = _para_CLK[i][0];
			_trp[i] = _para_TRP[i][0];
			if (!_lost_isb_GAL[i])
				_isb_GAL[i] = _para_ISB_GAL[i][0];
			if (!_lost_isb_BDS[i])
				_isb_BDS[i] = _para_ISB_BDS[i][0];
			if (!_lost_isb_GLO[i])
				_isb_GLO[i] = _para_ISB_GLO[i][0];
		}
	}

}

void gfgomsf::t_gpvtfgo::_vector_to_double()
{
	for (int i = 0; i <= _rover_count; i++)
	{
		_para_CRD[i][0] = _Pos[i].x();
		_para_CRD[i][1] = _Pos[i].y();
		_para_CRD[i][2] = _Pos[i].z();

	}
	if (_isBase) {
		for (int i = 0; i < _amb_manager->ambiguity_ids.size(); i++)
		{
			int amb_id = _amb_manager->ambiguity_ids[i];
			_para_amb[amb_id][0] = _amb_manager->getAmb(amb_id);
		}
		_amb_manager->generateAmbSearchIndex();
	}
	else
	{
		for (int i = 0; i < _ambIF_manager->ambiguity_ids.size(); i++)
		{
			int amb_id = _ambIF_manager->ambiguity_ids[i];
			_para_AMB_IF[amb_id][0] = _ambIF_manager->getAmb(amb_id);
		}
		_ambIF_manager->generateAmbSearchIndex();
	}
	//for PPP IF:
	if (!_isBase)
	{
		for (int i = 0; i <= _rover_count; i++)
		{
			_para_CLK[i][0] = _clk[i];
			_para_TRP[i][0] = _trp[i];
			if (!_lost_isb_GAL[i])
				_para_ISB_GAL[i][0] = _isb_GAL[i];
			if (!_lost_isb_BDS[i])
				_para_ISB_BDS[i][0] = _isb_BDS[i];
			if (!_lost_isb_GLO[i])
				_para_ISB_GLO[i][0] = _isb_GLO[i];
		}
	}
}

void gfgomsf::t_gpvtfgo::_select_ref_sat()
{
	_DD_msg.clear();
	//_obs_index.clear();
	_sat_ref.clear();
	set<string> sysall = dynamic_cast<t_gsetgen *>(_set)->sys();
	bool isSetRefSat = dynamic_cast<t_gsetamb *>(_set)->isSetRefSat();
	bool isPhaseProcess = true;
	for (int obslevel = _obs_level; obslevel <= 3; obslevel++)
	{
		for (auto sys_iter = sysall.begin(); sys_iter != sysall.end(); sys_iter++)
		{
			enum GSYS sys = t_gsys::str2gsys(*sys_iter);
			vector<GOBSBAND> band = dynamic_cast<t_gsetgnss *>(_set)->band(sys);
			int nf = 5;
			if (band.size())
				nf = band.size();
			if (_observ == OBSCOMBIN::IONO_FREE)
				nf = 1;
			FREQ_SEQ f;
			string sat_ref;
			for (FREQ_SEQ freq = FREQ_1; freq <= 2 * nf; freq = (FREQ_SEQ)(freq + 1))
			{
				if (freq <= nf)
				{ //phase equations
					isPhaseProcess = true;
					f = freq;
				}
				else
				{ //code equations
					isPhaseProcess = false;
					f = (FREQ_SEQ)(freq - nf);
				}
				if (f > _frequency)
					continue;
				string sat;
				t_gsatdata obs_sat_ref;
				enum GSYS gs;
				if (!isSetRefSat || (_observ == OBSCOMBIN::RAW_MIX && !isPhaseProcess))
					sat_ref.clear();
				//sat_ref.empty();
				if (sat_ref.empty())
				{
					for (auto it = _data.begin(); it != _data.end(); ++it)
					{
						std::string sat_id = it->sat();   // Obtain sat ID
						double elevation = it->ele();       // Obtain the elevation angle

					}
					for (auto it = _data.begin(); it != _data.end(); it++)
					{
						string sat = it->sat();

						gs = it->gsys();
						if (gs == QZS)
							gs = GPS;
						if (gs != sys)
							continue;
						if ((gs == BDS) && t_gsys::bds_geo(sat))
							continue;
						if (/*(_observ == RAW_ALL||_observ == IONO_FREE||_observ == RAW_DOUBLE)
						//	 &&*/
							!_reset_amb && !_reset_par && it->islip())
							continue;

						GOBSBAND b = _band_index[gs][f];
						if (isPhaseProcess)
						{
							if (!it->band_avail(true).count(b))
							{
								continue;
							}
						}
						else
						{
							if (!it->band_avail(true).count(b) || !it->band_avail(false).count(b))
							{
								continue;
							}
						}
						int base_flag = 0;
						for (auto it_base = _data_base.begin(); it_base != _data_base.end(); it_base++)
						{
							if (it_base->sat() != sat)
							{
								continue;
							}

							if (isPhaseProcess)
							{
								if (it_base->band_avail(true).count(b))
								{
									base_flag = 1;
								}
							}
							else
							{
								if (it_base->band_avail(true).count(b) && it_base->band_avail(false).count(b))
								{
									base_flag = 1;
								}
							}
						}
						if (!base_flag)
						{
							continue;
						}


						if (sat_ref.empty())
						{
							sat_ref = sat;
							obs_sat_ref = *it;
							continue;
						}

						double e = it->ele_deg();
						double e2 = obs_sat_ref.ele_deg();
						if (e>=e2)
						{
							sat_ref = sat;
							sat_ref = it->sat();
							obs_sat_ref = *it;
						}

					} //end select sat_ref
				}
				if (sat_ref.empty())
					continue;
				if (_ipSatRep[sys] != "" && sat_ref != _ipSatRep[sys])
				{
					cerr << "refsat bug" << endl;
					continue;
				}
				if (freq == FREQ_1)
					_sat_ref.insert(sat_ref);
				t_gsatdata ref_sat_data;
				auto it = find_if(_data.begin(), _data.end(), [sat_ref](t_gsatdata it)
				{
					return it.sat() == sat_ref;
				});
				assert(it != _data.end());
				ref_sat_data = *it;
				for (auto it = _data.begin(); it != _data.end(); it++)
				{
					sat = it->sat();
					if (sat == sat_ref)
					{
						continue;
					}
					gs = it->gsys();
					if (gs == QZS)
						gs = GPS;
					if (gs != sys)
						continue;
					GOBSTYPE obstype = TYPE_C;
					if (isPhaseProcess)
						obstype = GOBSTYPE::TYPE_L;
					//_obs_index.push_back(make_pair(it->sat(), make_pair(f, obstype)));
					// delete unrecorded observations, added by hyChang
					GOBSBAND b = _band_index[gs][f];
					if (isPhaseProcess)
					{
						if (!it->band_avail(true).count(b))
						{
							continue;
						}
					}
					else
					{
						if (!it->band_avail(true).count(b) || !it->band_avail(false).count(b))
						{
							continue;
						}
					}
					int base_flag = 0;
					for (auto it_base = _data_base.begin(); it_base != _data_base.end(); it_base++)
					{
						if (it_base->sat() != sat)
						{
							continue;
						}

						if (isPhaseProcess)
						{
							if (it_base->band_avail(true).count(b))
							{
								base_flag = 1;
							}
						}
						else
						{
							if (it_base->band_avail(true).count(b) && it_base->band_avail(false).count(b))
							{
								base_flag = 1;
							}
						}
					}
					if (!base_flag)
					{
						continue;
					}

					DDEquMsg dd_msg(ref_sat_data, *it, obstype, f);
					//if (gins_window_size>0)
					//{
					//	dd_msg.ref_sat_global_id = cur_sat_prn[ref_sat_data.sat()];
					//	dd_msg.nonref_sat_global_id = cur_sat_prn[it->sat()];
					//}
					//else
					//{
						dd_msg.ref_sat_global_id = _amb_manager->get_sat_id(ref_sat_data.sat());
						//dd_msg.ref_sat_global_id = cur_sat_prn[ref_sat_data.sat()];
						dd_msg.nonref_sat_global_id = _amb_manager->get_sat_id(it->sat());
					//}

					assert(dd_msg.ref_sat_global_id != -1);
					assert(dd_msg.nonref_sat_global_id != -1);

					_DD_msg.push_back(dd_msg);

				} //end sat
			}      //end f
		}          //end sys
	}
}



bool gfgomsf::t_gpvtfgo::_get_DD_data(DDEquMsg & dd_msg,  vector<t_gsatdata> base_sat_data)
{
	t_gsatdata rover_ref_sat = dd_msg.rover_ref_sat;
	t_gsatdata rover_nonref_sat = dd_msg.rover_nonref_sat;
	FREQ_SEQ   freq = dd_msg.freq;
	GOBSTYPE   obstype = dd_msg.obs_type;
	t_gsatdata base_ref_sat, base_nonref_sat;
	if (base_sat_data.empty()) return false;
	int ref_i = 0;
	int nonref_i = 0;

	for (auto it : base_sat_data)
	{
		if (it.sat() == rover_ref_sat.sat())
		{
			ref_i = 1;
			base_ref_sat = it;
		}
		if (it.sat() == rover_nonref_sat.sat())
		{
			nonref_i = 1;
			base_nonref_sat = it;
		}
	}
	if (!ref_i || !nonref_i)
	{
		cout << "no common view between base and rover!!!" << endl;
		return false;
	}
	map<FREQ_SEQ, GOBSBAND> crt_bands = _band_index[rover_ref_sat.gsys()];
	if (crt_bands.empty()) return false;
	if (freq > _frequency) return false;
	dd_msg.band = crt_bands[freq];
	dd_msg.base_ref_sat = base_ref_sat;
	dd_msg.base_nonref_sat = base_nonref_sat;
	dd_msg.time = rover_ref_sat.epoch();
	dd_msg.base_site = _site_base;
	dd_msg.rover_site = _site;
	return true;
}


bool gfgomsf::t_gpvtfgo::_remove_outlier_sat(const pair<string, int>& outlier)
{
    //attention: the newest sat observation will be removed
    if (_isBase && outlier.first != " ")
    {
    	if (_amb_manager->cur_sats.size()  > _minsat)
    	{
    		pair<string, int> sat_id = outlier;
    		auto it_DD = _vDD_msg[_rover_count].begin();
    		int sat_global_id = outlier.second;
    		/*
			for (int i = _rover_count; i <= _rover_count; i++)
			{
				for (vector<DDEquMsg>::iterator it_DD = _vDD_msg[i].begin(); it_DD != _vDD_msg[i].end(); )
				{
					if (it_DD->nonref_sat_global_id == sat_id.second)
					{
						it_DD = _vDD_msg[i].erase(it_DD);
					}
					else
					{
						++it_DD;
					}
				}
			}
			*/

    		while (it_DD != _vDD_msg[_rover_count].end())
    		{
    			if (it_DD->ref_sat_global_id == sat_global_id || it_DD->nonref_sat_global_id == sat_global_id)
    			{
    				_vDD_msg[_rover_count].erase(it_DD);
    				continue;
    			}

    			it_DD++;

    		}


    		_amb_manager->removeSat(sat_global_id, _rover_count);
    		_amb_manager->get_last_epoch_sats(_epoch.sow());//update cur sat map

    		//update cur sat map
    		//auto it = cur_sat_prn.find(sat_id.first);
    		//if (it != cur_sat_prn.end())
    		//	cur_sat_prn.erase(sat_id.first);
    		_last_gnss_info->valid = false;
    		auto it = cur_sat_prn.find(sat_id.first);
    		if (it != cur_sat_prn.end())
    			cur_sat_prn.erase(sat_id.first);
    		return true;
    	}
    	else return false;
    }

    if (!_isBase && outlier.first != " ")
    {
        if (_ambIF_manager->cur_sats.size() > _minsat)//add for MultiWindow
        {
            pair<string, int> sat_id = outlier;
            //add for MultiWindow
            auto it_IF = _vIF_msg[_rover_count].begin();
            int sat_global_id = outlier.second;
            while (it_IF != _vIF_msg[_rover_count].end())
            {
                if (it_IF->sat_id == sat_global_id)
                {
                    //_vIF_msg[_rover_count].erase(it_IF);
                	it_IF = _vIF_msg[_rover_count].erase(it_IF);
                    continue;
                }
                it_IF++;
            }

            //_ambIF_manager->removeSat_beta(sat_id.second, _rover_count);//delete for MultiWindow
            _ambIF_manager->removeSat(sat_global_id, _rover_count);//add for MultiWindow

            _ambIF_manager->get_last_epoch_sats(_epoch.sow());
            _last_gnss_info->valid = false;
            //update sys info
            if (!_lost_isb_GAL[_rover_count] && outlier.first[0] == 'E')//with GAL sats and delete a GAL sat
            {
                bool is_out = true;
                //for (auto iter = cur_sat_prn.begin(); iter != cur_sat_prn.end(); iter++)//delete for MultiWindow
                for (auto iter = _ambIF_manager->cur_sats.begin(); iter != _ambIF_manager->cur_sats.end(); iter++)//add for MultiWindow
                {
                    if (iter->first[0] == 'E')//still with GAL sats
                    {
                        is_out = false;
                        break;
                    }
                }
                _lost_isb_GAL[_rover_count] = is_out;
            }
            if (!_lost_isb_BDS[_rover_count] && outlier.first[0] == 'C')//with BDS sats and delete a BDS sat
            {
                bool is_out = true;
                //for (auto iter = cur_sat_prn.begin(); iter != cur_sat_prn.end(); iter++)//delete for MultiWindow
                for (auto iter = _ambIF_manager->cur_sats.begin(); iter != _ambIF_manager->cur_sats.end(); iter++)//add for MultiWindow
                {
                    if (iter->first[0] == 'C')//still with BDS sats
                    {
                        is_out = false;
                        break;
                    }
                }
                _lost_isb_BDS[_rover_count] = is_out;
            }
			if (!_lost_isb_GLO[_rover_count] && outlier.first[0] == 'R')
			{
				bool is_out = true;
				for (auto iter = _ambIF_manager->cur_sats.begin(); iter != _ambIF_manager->cur_sats.end(); iter++)
				{
					if (!iter->first.empty() && iter->first[0] == 'R')
					{
						is_out = false;
						break;
					}
				}
				_lost_isb_GLO[_rover_count] = is_out;
			}
            return true;
        }
        else
            return false;
    }

    return true;
}
bool gfgomsf::t_gpvtfgo::_pre_amb_resolution()
{
	t_gallpar construct_para = _all_para_win;
	int nobs_total, npar_number;
	Matrix A_fgo;
	SymmetricMatrix P_fgo;
	ColumnVector l_fgo, dx_fgo;
	SymmetricMatrix Qx0_fgo,Qx_fgo;
	double vtpv_fgo;
	nobs_total = _last_gnss_info->linearized_jacobians.rows();
	npar_number = _last_gnss_info->linearized_jacobians.cols();
	assert(npar_number == construct_para.parNumber());
	//cout << "nobs_total：" << nobs_total << endl;
	//cout << "npar_number：" << npar_number << endl;
	A_fgo.ReSize(nobs_total,npar_number);
	A_fgo = 0.0;
	P_fgo.ReSize(nobs_total);
	P_fgo = 0.0;
	l_fgo.ReSize(nobs_total);
	l_fgo = 0.0;
	dx_fgo.ReSize(npar_number);
	dx_fgo = 0.0;
	Qx0_fgo.ReSize(npar_number);
	Qx0_fgo = 0.0;
	Qx_fgo.ReSize(npar_number);
	Qx_fgo = 0.0;
	_sig_unit= _last_gnss_info->sig_unit;
	vtpv_fgo = _last_gnss_info->vtpv;
	//for Qx
	for (int i = 0; i < npar_number; i++)
	{
		for (int j = 0; j < npar_number; j++)
		{
			//Qx_fgo(i + 1, j + 1) = _last_gnss_info->Qx(i, j);
			Qx0_fgo(i + 1, j + 1) = _last_gnss_info->Qx(i, j);
		}
	}
	Qx_fgo = Qx0_fgo;
	//for A
	for (int i = 0; i < nobs_total; i++)
	{
		for (int j = 0; j < npar_number; j++)
		{
			A_fgo(i + 1, j + 1) = _last_gnss_info->linearized_jacobians(i, j);
		}
	}
	//for P
	for (int i = 0; i < nobs_total; i++)
	{
		for (int j = 0; j < nobs_total; j++)
		{
			P_fgo(i + 1, j + 1) = _last_gnss_info->weight(i,j);
		}
	}
	//for l
	for (int i = 0; i < nobs_total; i++)
	{
		l_fgo(i + 1) = _last_gnss_info->linearized_residuals(i);
	}
	_filter->add_data(construct_para, dx_fgo, Qx_fgo, _sig_unit, Qx0_fgo);
	_filter->add_data(A_fgo, P_fgo, l_fgo);
	_filter->add_data(vtpv_fgo, nobs_total, npar_number);

	return true;;
}

int gfgomsf::t_gpvtfgo::_optimization()
{
	_removed_sats.clear();
	t_tictoc fgo_gnss;
	int count = 0;
	bool iter_flag = false;
	pair<string, int>  outlier = make_pair(" ", -1);
	do
	{
		count++;
		if (!_remove_outlier_sat(outlier))
			break;
		_vector_to_double();
		ceres::Problem problem;
		ceres::LossFunction *loss_function;
		//loss_function = new ceres::HuberLoss(_loss_func_value);
		// change
		loss_function = new ceres::CauchyLoss(2);
		for (int i = 0; i < _rover_count + 1; i++)
		{
			problem.AddParameterBlock(_para_CRD[i], 3);
		}

		_prior_factor(problem);
		//_gnss_obs_index.clear();
		assert(_vDD_msg.size() == _rover_count + 1);
		for (int i = 0; i <= _rover_count; i++)
		{
			t_gallpar params_temp(_para_window[i]);
			vector<DDEquMsg> DD_tmp = _vDD_msg[i];
			vector<t_gsatdata> b_sat_data = _win_base_data[i];
			int dd_equ_count = 0;
			for (auto &dd_iter : DD_tmp)
			{
				if (!_get_DD_data(dd_iter, b_sat_data)) continue;
				pair<string, string> base_rover_site = make_pair(dd_iter.base_site, dd_iter.rover_site);
				pair<FREQ_SEQ, GOBSBAND> freq_band = make_pair(dd_iter.freq, dd_iter.band);
				vector<pair<t_gsatdata, t_gsatdata>> DD_sat_data;
				DD_sat_data.push_back(make_pair(dd_iter.base_ref_sat, dd_iter.rover_ref_sat));  // reference sat
				DD_sat_data.push_back(make_pair(dd_iter.base_nonref_sat, dd_iter.rover_nonref_sat));  // unreference sat
				GOBSTYPE  obstype = dd_iter.obs_type;
				//_gnss_obs_index.push_back(make_pair(make_pair(dd_iter.rover_nonref_sat.sat(), dd_iter.nonref_sat_global_id), make_pair(dd_iter.freq, obstype)));
				if (obstype == GOBSTYPE::TYPE_C)
				{
					PseudorangeDDFactor *pf = new PseudorangeDDFactor(dd_iter.time, base_rover_site, params_temp, DD_sat_data, _gbias_model, freq_band);//DD_Pseudorange
					window_pseudo_factors.push_back(pf);
					problem.AddResidualBlock(pf, NULL, _para_CRD[i]);
				}
				if (obstype == GOBSTYPE::TYPE_L)
				{
					int id1, id2;
					id1 = _amb_manager->getAmbSearchIndex(make_pair(dd_iter.ref_sat_global_id, dd_iter.freq));
					id2 = _amb_manager->getAmbSearchIndex(make_pair(dd_iter.nonref_sat_global_id, dd_iter.freq));

					if (id1 == -1 || id2 == -1)
						continue;

					CarrierphaseDDFactor *lf = new CarrierphaseDDFactor(dd_iter.time, base_rover_site, params_temp, DD_sat_data, _gbias_model, freq_band);//DD_Carrierphase
					window_carrierphase_dd_factors.push_back(lf);
					problem.AddResidualBlock(lf, NULL, _para_CRD[i], _para_amb[id1], _para_amb[id2]);
				}
				dd_equ_count++;
			}
			assert(dd_equ_count == DD_tmp.size());
		}
		//if (_initial_prior)

		if (1)
		{
			for (int i = 0; i < _amb_manager->ambiguity_ids.size(); i++)
			{
				int amb_id = _amb_manager->ambiguity_ids[i];

				double initial_amb = _amb_manager->getInitialAmb(amb_id);
				InitialGnssAMB *amb_prior = new InitialGnssAMB(initial_amb);
				problem.AddResidualBlock(amb_prior, NULL, _para_amb[amb_id]);
			}
		}
		for (int i = 0; i < _rover_count + 1; i++)
		{

			InitialPosFactor *initial_pos = new InitialPosFactor(Eigen::Vector3d(_para_CRD[i]));//initial position factor
			//initial_pos->sqrt_info = 1e-4 * Eigen::Matrix<double, 3, 3>::Identity();
			problem.AddResidualBlock(initial_pos, NULL, _para_CRD[i]);
		}
		//ceres solver
		ceres::Solver::Options options;
		options.linear_solver_type = ceres::DENSE_QR;
		options.max_num_iterations = 10;
		options.trust_region_strategy_type = ceres::DOGLEG;
		ceres::Solver::Summary summary;

		//for (auto* f : window_pseudo_factors) {
		//	f->EnableFreezeJacobian(true);
		//	f->ResetJacobianCache();
		//}
		//for (auto* f : window_carrierphase_dd_factors) {
		//	f->EnableFreezeJacobian(true);
		//	f->ResetJacobianCache();
		//}

		cout << "begin!!!" << endl;
		ceres::Solve(options, &problem, &summary);
		cout << "end!!!" << endl;

		_posteriori_test(problem);
		double idx = _gobs_outlier_detection(outlier);
		if (idx >= 0)
			iter_flag = true;
		else
			iter_flag = false;

	} while (iter_flag);

	std::cout << "Epoch: " << _headers[_rover_count] << " time cost (ms): " << fgo_gnss.toc() <<" iter times: "<<count<< endl;
	if (_last_gnss_info->valid)
	{
		std::cout<<"cur epoch: " << "sigma: " << _last_gnss_info->sig_unit << "  vtpv: " <<  _last_gnss_info->vtpv <<std::endl;
	}
	else
	{
		cout << "soving failed in this epoch!!!" << endl;
	}
	_double_to_vector();
	return 1;
}

int gfgomsf::t_gpvtfgo::_optimization_PPP()
{
    _removed_sats.clear();
    t_tictoc fgo_gnss;
    int count = 0;
    int ceres_iteration = -1;
    double ceres_solving_cost = -1;
    int ceres_state = -1;
    bool iter_flag = false;
    pair<string, int>  outlier = make_pair(" ", -1);
    do
    {
        count++;

        if (!_remove_outlier_sat(outlier))
        {
			// An actual outlier was found, but removing it would leave fewer
			// than min_sat satellites.  Match PPP/INS FGO semantics: reject the
			// entire candidate epoch instead of accepting the bad four-satellite
			// solution into the optimizer window.
			if (outlier.first != " ")
			{
				if (_last_gnss_info)
					_last_gnss_info->valid = false;
				return -1;
			}
            break;
        }
        _vector_to_double();//checked
        ceres::Problem problem;
        ceres::LossFunction* loss_function;
        ceres::LossFunction* loss_function_CP;
        //loss_function = NULL ;
        //loss_function_CP = NULL;
    	loss_function = new ceres::HuberLoss(_loss_func_value);
    	loss_function_CP = new ceres::HuberLoss(_loss_func_value);
        for (int i = 0; i < _rover_count + 1; i++)
        {
            problem.AddParameterBlock(_para_CRD[i], 3);
            problem.AddParameterBlock(_para_CLK[i], 1);
            problem.AddParameterBlock(_para_TRP[i], 1);
            if (!_lost_isb_GAL[i])
                problem.AddParameterBlock(_para_ISB_GAL[i], 1);
            if (!_lost_isb_BDS[i])
                problem.AddParameterBlock(_para_ISB_BDS[i], 1);
			if (!_lost_isb_GLO[i])
				problem.AddParameterBlock(_para_ISB_GLO[i], 1);
        }

        //Posterior Factor
        //_posterior_factor(problem);
        _prior_factor(problem);//checked

        //Troposphere RW Factor
        for (int i = 0; i <= _rover_count; i++) {
	        if (i + 1 > _rover_count)
	        {
	        	continue;
	        }
        	const double graph_dt = std::fabs(_headers[i + 1] - _headers[i]);
        	double sqrt_info = 1.0 / sqrt(graph_interval_random_walk_q(_trpStoModel, graph_dt));
        	RandomWalkFactor* rf = new RandomWalkFactor(sqrt_info);
        	problem.AddResidualBlock(rf, NULL, _para_TRP[i], _para_TRP[i + 1]);
        }

        for (int i = 0; i <= _rover_count; i++)
        {
            double sqrt_info = 1.0 / _sig_init_ztd;
            //InitialFactor* itf = new InitialFactor(_trp[i], sqrt_info);
            InitialFactor* itf = new InitialFactor(_trp_ini, sqrt_info);
            problem.AddResidualBlock(itf, NULL, _para_TRP[i]);
        }
        //ISB RW Factor
        for (int i = 0; i <= _rover_count; i++)
        {
            if (i + 1 > _rover_count)
            {
                continue;
            }
            const double graph_dt = std::fabs(_headers[i + 1] - _headers[i]);
            if (!_lost_isb_GAL[i] && !_lost_isb_GAL[i + 1])
            {
                double sqrt_info_gal = 1.0 / sqrt(graph_interval_random_walk_q(_galStoModel, graph_dt));
                RandomWalkFactor* galf = new RandomWalkFactor(sqrt_info_gal);
                problem.AddResidualBlock(galf, NULL, _para_ISB_GAL[i], _para_ISB_GAL[i + 1]);
            }
            if (!_lost_isb_BDS[i] && !_lost_isb_BDS[i + 1])
            {
                double sqrt_info_bds = 1.0 / sqrt(graph_interval_random_walk_q(_bdsStoModel, graph_dt));
                RandomWalkFactor* bdsf = new RandomWalkFactor(sqrt_info_bds);
                problem.AddResidualBlock(bdsf, NULL, _para_ISB_BDS[i], _para_ISB_BDS[i + 1]);
            }
			if (!_lost_isb_GLO[i] && !_lost_isb_GLO[i + 1])
			{
				double sqrt_info_glo = 1.0 / sqrt(graph_interval_random_walk_q(_gloStoModel, graph_dt));
				RandomWalkFactor* glof = new RandomWalkFactor(sqrt_info_glo);
				problem.AddResidualBlock(glof, NULL, _para_ISB_GLO[i], _para_ISB_GLO[i + 1]);
			}
        }
        //ISB Initial Factor
        for (int i = 0; i <= _rover_count; i++)
        {
            if (!_lost_isb_GAL[i])
            {
                //InitialFactor* igalf = new InitialFactor(_isb_GAL[i], 1.0 / _sig_init_gal);
                InitialFactor* igalf = new InitialFactor(0.0, 1.0 / _sig_init_gal);
                problem.AddResidualBlock(igalf, NULL, _para_ISB_GAL[i]);
            }
            if (!_lost_isb_BDS[i])
            {
                //InitialFactor* ibdsf = new InitialFactor(_isb_BDS[i], 1.0 / _sig_init_bds);
                InitialFactor* ibdsf = new InitialFactor(0.0, 1.0 / _sig_init_bds);
                problem.AddResidualBlock(ibdsf, NULL, _para_ISB_BDS[i]);
            }
			if (!_lost_isb_GLO[i])
			{
				InitialFactor* iglof = new InitialFactor(0.0, 1.0 / _sig_init_glo);
				problem.AddResidualBlock(iglof, NULL, _para_ISB_GLO[i]);
			}
        }
        //Ambiguity Initial Factor
        for (int i = 0; i < _ambIF_manager->ambiguity_ids.size(); i++)
        {
            double initial_amb;
            int amb_id = _ambIF_manager->ambiguity_ids[i];
            initial_amb = _ambIF_manager->getInitialAmb(amb_id);
            InitialGnssAMB* amb_prior = new InitialGnssAMB(initial_amb, 1.0 / _sigAmbig);
            problem.AddResidualBlock(amb_prior, NULL, _para_AMB_IF[amb_id]);
        }
        //Position SPP Factor
        for (int i = 0; i <= _rover_count; i++)
        {
            Eigen::Vector3d xyz = _Pos[i];
            InitialGnssCRD* pos_prior = new InitialGnssCRD(xyz, 1.0 / _sig_init_crd);
            problem.AddResidualBlock(pos_prior, NULL, _para_CRD[i]);
        }
        //Clock Initial Factor
        for (int i = 0; i <= _rover_count; i++)
        {
            double clk = _clk[i];
            InitialFactor* icf = new InitialFactor(clk, 1.0 / _clkStoModel->getQ() );// sqrt(_clkStoModel->getQ()) error
            problem.AddResidualBlock(icf, NULL, _para_CLK[i]);
        }
        //Pseudorange & Carrierphase IONO-FREE Factor
        assert(_vIF_msg.size() == _rover_count + 1);
        for (int i = 0; i <= _rover_count; i++)
        {
            t_gallpar params_temp(_para_window[i]);
            vector<IFEquMsg> IF_tmp = _vIF_msg[i];
            int if_equ_num = 0;
            for (auto& if_iter : IF_tmp)
            {
                GOBSTYPE obstype = if_iter.obs_type;
                string sat_name = if_iter.satdata.sat();
                GSYS gsys = if_iter.satdata.gsys();
                if (obstype == GOBSTYPE::TYPE_C)
                {
                    pair<FREQ_SEQ, GOBSBAND> freq_band1 = make_pair(if_iter.freq1, if_iter.band1);
                    pair<FREQ_SEQ, GOBSBAND> freq_band2 = make_pair(if_iter.freq2, if_iter.band2);
                    if (gsys == GSYS::GPS)
                    {
                        PseudorangeIFFactor* pf = new PseudorangeIFFactor(if_iter.time, if_iter.site, params_temp, if_iter.satdata, _gbias_model, freq_band1, freq_band2);
                        problem.AddResidualBlock(pf, loss_function, _para_CRD[i], _para_CLK[i], _para_TRP[i]);
                    }
                    else if (gsys == GSYS::GAL && !_lost_isb_GAL[i])
                    {
                        MultiPseudorangeIFFactor* pf = new MultiPseudorangeIFFactor(if_iter.time, if_iter.site, params_temp, if_iter.satdata, _gbias_model, freq_band1, freq_band2);
                        problem.AddResidualBlock(pf, loss_function, _para_CRD[i], _para_CLK[i], _para_TRP[i], _para_ISB_GAL[i]);
                    }
                    else if (gsys == GSYS::BDS && !_lost_isb_BDS[i])
                    {
                        MultiPseudorangeIFFactor* pf = new MultiPseudorangeIFFactor(if_iter.time, if_iter.site, params_temp, if_iter.satdata, _gbias_model, freq_band1, freq_band2);
                        problem.AddResidualBlock(pf, loss_function, _para_CRD[i], _para_CLK[i], _para_TRP[i], _para_ISB_BDS[i]);
                    }
					else if (gsys == GSYS::GLO && !_lost_isb_GLO[i])
					{
						MultiPseudorangeIFFactor* pf = new MultiPseudorangeIFFactor(if_iter.time, if_iter.site, params_temp, if_iter.satdata, _gbias_model, freq_band1, freq_band2);
						problem.AddResidualBlock(pf, loss_function, _para_CRD[i], _para_CLK[i], _para_TRP[i], _para_ISB_GLO[i]);
					}
                    else
                    {
                        cerr << "not support the system now : " << gsys << endl;
                    }
                }
                if (obstype == GOBSTYPE::TYPE_L)
                {
                    int id = _ambIF_manager->getAmbSearchIndex(if_iter.sat_id);
                    if (id == -1)
                    {
                        continue;
                    }

                    pair<FREQ_SEQ, GOBSBAND> freq_band1 = make_pair(if_iter.freq1, if_iter.band1);
                    pair<FREQ_SEQ, GOBSBAND> freq_band2 = make_pair(if_iter.freq2, if_iter.band2);
                    if (gsys == GSYS::GPS)
                    {
                        CarrierphaseIFFactor* lf = new CarrierphaseIFFactor(if_iter.time, if_iter.site, params_temp, if_iter.satdata, _gbias_model, freq_band1, freq_band2);
                        problem.AddResidualBlock(lf, loss_function_CP, _para_CRD[i], _para_CLK[i], _para_TRP[i], _para_AMB_IF[id]);
                    }
                    else if (gsys == GSYS::GAL && !_lost_isb_GAL[i])
                    {
                        MultiCarrierphaseIFFactor* lf = new MultiCarrierphaseIFFactor(if_iter.time, if_iter.site, params_temp, if_iter.satdata, _gbias_model, freq_band1, freq_band2);
                        problem.AddResidualBlock(lf, loss_function_CP, _para_CRD[i], _para_CLK[i], _para_TRP[i], _para_ISB_GAL[i], _para_AMB_IF[id]);
                    }
                    else if (gsys == GSYS::BDS && !_lost_isb_BDS[i])
                    {
                        MultiCarrierphaseIFFactor* lf = new MultiCarrierphaseIFFactor(if_iter.time, if_iter.site, params_temp, if_iter.satdata, _gbias_model, freq_band1, freq_band2);
                        problem.AddResidualBlock(lf, loss_function_CP, _para_CRD[i], _para_CLK[i], _para_TRP[i], _para_ISB_BDS[i], _para_AMB_IF[id]);
                    }
					else if (gsys == GSYS::GLO && !_lost_isb_GLO[i])
					{
						MultiCarrierphaseIFFactor* lf = new MultiCarrierphaseIFFactor(if_iter.time, if_iter.site, params_temp, if_iter.satdata, _gbias_model, freq_band1, freq_band2);
						problem.AddResidualBlock(lf, loss_function_CP, _para_CRD[i], _para_CLK[i], _para_TRP[i], _para_ISB_GLO[i], _para_AMB_IF[id]);
					}
                    else
                    {
                        cerr << "not support the system now : " << gsys << endl;
                    }
                    //cout << " sat_name : " << sat_name << " sat_id : " << if_iter.sat_id << " _para_AMB_IF[" << id << "] = " << _para_AMB_IF[id][0] << endl;
                }
                if_equ_num++;
            }
            //assert(if_equ_count == IF_tmp.size());
        }

        // if (false)
        // {
        //     for (int i = 0; i <= _rover_count; i++)
        //     {
        //         double sow = _vIF_msg[i][0].time.sow();
        //         string str_sow = to_string(sow);
        //         str_sow = str_sow.substr(0, 7);
        //         if (_txyz_pvt.count(str_sow) != 0)
        //         {
        //             Eigen::Vector3d pos = _txyz_pvt[str_sow];
        //             Eigen::Vector3d sqrt;
        //             sqrt << 1 / 1, 1 / 1, 1 / 1;
        //             InitialGnssCRD* posf = new InitialGnssCRD(pos, sqrt);
        //             problem.AddResidualBlock(posf, loss_function, _para_CRD[i]);
        //         }
        //     }
        // }

        //ceres solver
        ceres::Solver::Options options;
        options.linear_solver_type = ceres::DENSE_QR;
        options.max_num_iterations = _max_num_iterations;
        //options.parameter_tolerance = 1.0e-10;  //default 1e-8
        options.max_solver_time_in_seconds = 1;
        //options.minimizer_progress_to_stdout = true;
        //options.num_threads = 2; !!!can not be used, caused by newmat
        options.trust_region_strategy_type = ceres::DOGLEG;
        ceres::Solver::Summary summary;
        ceres::Solve(options, &problem, &summary);
        //cout << summary.FullReport() << endl;
        //std::cout << summary.BriefReport() << ", " << "Total time: " << summary.total_time_in_seconds << endl;
        //add for MultiWindow
        if (summary.termination_type == ceres::FAILURE)
        {
            std::cout << "#Epoch: " << _headers[_rover_count] << "ceres FAILURE!!!" << endl;
        }
        ceres_iteration = summary.iterations.size();
        ceres_solving_cost = summary.total_time_in_seconds;
        ceres_state = summary.termination_type;
        //_posteriori_test_PPP();
        //_posteriori_test_PPP(problem);
    	_posteriori_test_PPP(problem);//modify internally for MultiWindow
        if (_gobs_outlier_detection(outlier) >= 0)//checked
            iter_flag = true;
        else
            iter_flag = false;
    } while (iter_flag);

    // //add for MultiWindow
    // _output_estimator_info << setw(15) << fixed << setprecision(1) << _headers[_rover_count] <<
    //     " " << setw(15) << ceres_iteration <<
    //     " " << fixed << setprecision(4) << setw(15) << ceres_solving_cost <<
    //     " " << setw(15) << ceres_state <<
    //     " " << setw(15) << count <<
    //     " " << fixed << setprecision(4) << setw(15) << fgo_gnss.toc() / 1000.0 << endl;
    if (_last_gnss_info->valid)
    {
        _double_to_vector(); //only passed check will be used for update//checked
        std::cout << "sigma: " << _last_gnss_info->sig_unit << "  vtpv: " << _last_gnss_info->vtpv << std::endl;
    }
    else
        std::cout << "Not enough IF observations!!!" << std::endl;

    return 1;
}

void gfgomsf::t_gpvtfgo::_record_ppp_cycle_slips(
	const std::vector<t_gsatdata>& epoch_data)
{
	if (_isBase || !_phase)
		return;

	const double nominal_interval =
		_sampling > 0.0 ? _sampling : 1.0;

	for (const auto& sat_data : epoch_data)
	{
		const std::string sat_name = sat_data.sat();
		const std::vector<GOBSBAND> bands =
			dynamic_cast<t_gsetgnss*>(_set)->band(sat_data.gsys());
		const GSYS system = sat_data.gsys();
		const int frequency_count =
			bands.empty() ? 5 : static_cast<int>(bands.size());
		int required_phase_count = 0;
		int valid_phase_count = 0;
		bool lli_cycle_slip = false;

		for (FREQ_SEQ frequency = FREQ_1;
			 frequency <= frequency_count;
			 frequency = static_cast<FREQ_SEQ>(frequency + 1))
		{
			if (frequency > _frequency)
				continue;

			required_phase_count++;
			const GOBSBAND band =
				bands.size() >= static_cast<size_t>(frequency)
					? bands[frequency - 1]
					: t_gsys::band_priority(system, frequency);
			const t_gobs phase_obs(sat_data.select_phase(band));

			if (!double_eq(sat_data.obs_L(phase_obs), 0.0))
				valid_phase_count++;

			if (sat_data.getlli(phase_obs.gobs()) >= 1)
				lli_cycle_slip = true;
		}

		if (lli_cycle_slip)
		{
			const bool newly_recorded =
				_pending_ppp_slips.insert(sat_name).second;
			if (newly_recorded && _spdlog)
				_spdlog->info(
					"PPP FGO: pending cycle slip recorded for {} at {}",
					sat_name, sat_data.epoch().str_ymdhms());
		}

		const bool has_complete_phase =
			required_phase_count > 0 &&
			valid_phase_count == required_phase_count;
		if (!has_complete_phase)
			continue;

		const t_gtime phase_epoch = sat_data.epoch();
		const auto last_phase = _last_ppp_phase_epoch.find(sat_name);
		if (last_phase != _last_ppp_phase_epoch.end())
		{
			const double phase_gap = phase_epoch.diff(last_phase->second);
			if (phase_gap > 1.5 * nominal_interval)
			{
				const bool newly_recorded =
					_pending_ppp_slips.insert(sat_name).second;
				if (newly_recorded && _spdlog)
					_spdlog->info(
						"PPP FGO: pending ambiguity reset for {} at {} "
						"after {:.3f} s carrier gap",
						sat_name, phase_epoch.str_ymdhms(), phase_gap);
			}
			last_phase->second = phase_epoch;
		}
		else
		{
			_last_ppp_phase_epoch.insert(
				std::make_pair(sat_name, phase_epoch));
		}
	}
}

bool gfgomsf::t_gpvtfgo::_ppp_candidate_has_phase(
	const std::string& sat_name) const
{
	if (_rover_count < 0 ||
		_vIF_msg.size() <= static_cast<size_t>(_rover_count))
		return false;

	for (const auto& if_msg : _vIF_msg[_rover_count])
	{
		if (if_msg.obs_type == GOBSTYPE::TYPE_L &&
			if_msg.satdata.sat() == sat_name)
			return true;
	}
	return false;
}

void gfgomsf::t_gpvtfgo::_commit_ppp_cycle_slips()
{
	for (const auto& sat_name : _candidate_ppp_slips)
	{
		const auto active = std::find_if(
			_ambIF_manager->cur_sats.begin(),
			_ambIF_manager->cur_sats.end(),
			[&sat_name](const pair<string, int>& sat)
			{
				return sat.first == sat_name;
			});

		if (active != _ambIF_manager->cur_sats.end() &&
			_ppp_candidate_has_phase(sat_name))
		{
			_pending_ppp_slips.erase(sat_name);
		}
	}
	_candidate_ppp_slips.clear();
}

void gfgomsf::t_gpvtfgo::_rollback_current_ppp_node()
{
	const int failed_rover = _rover_count;
	if (_isBase || failed_rover < 0)
		return;

	// _set_initial_value() attached this rover index to every satellite in
	// the candidate epoch.  Remove those associations from cur_sats; a
	// satellite already rejected during the outlier loop has already been
	// detached and therefore is intentionally absent here.
	set<int> sat_ids;
	for (const auto& sat : _ambIF_manager->cur_sats)
		sat_ids.insert(sat.second);

	for (const int sat_id : sat_ids)
		_ambIF_manager->removeSat(sat_id, failed_rover);

	_IF_msg.clear();
	if (_vIF_msg.size() > static_cast<size_t>(failed_rover))
		_vIF_msg.pop_back();

	if (_rover_window[failed_rover] != nullptr)
		delete _rover_window[failed_rover];
	_rover_window[failed_rover] = nullptr;
	_para_window[failed_rover].delAllParam();
	_headers[failed_rover] = 0.0;
	_Pos[failed_rover].setZero();
	_Vel[failed_rover].setZero();
	_clk[failed_rover] = 0.0;
	_trp[failed_rover] = 0.0;
	_isb_GAL[failed_rover] = 0.0;
	_isb_BDS[failed_rover] = 0.0;
	_isb_GLO[failed_rover] = 0.0;
	_lost_isb_GAL[failed_rover] = false;
	_lost_isb_BDS[failed_rover] = false;
	_lost_isb_GLO[failed_rover] = false;

	if (failed_rover > 0)
		_ambIF_manager->get_last_epoch_sats(_headers[failed_rover - 1]);
	else
		_ambIF_manager->cur_sats.clear();
	_rover_count--;

	_global_sat_id = _ppp_candidate_global_sat_id;
	_global_amb_id = _ppp_candidate_global_amb_id;
	_candidate_ppp_slips.clear();

	if (_spdlog)
		_spdlog->warn("PPP FGO: rollback rejected GNSS epoch {}", _epoch.str_ymdhms());
}

void gfgomsf::t_gpvtfgo::_slide_window()
{
	if (_rover_count == gwindow_size-1)
	{
		// gwindow_size is the number of active nodes in this port. When the
		// last valid index is gwindow_size - 1, shift only surviving nodes.
		for (int i = 0; i < gwindow_size - 1; i++)
		{
			_headers[i] = _headers[i + 1];
			_Pos[i].swap(_Pos[i + 1]);
			//add
			_Vel[i] = _Vel[i + 1];
			std::swap(_rover_window[i], _rover_window[i + 1]);
			_win_base_data[i].swap(_win_base_data[i + 1]);
			_para_window[i] = _para_window[i + 1];
			//for PPP IF
			if (!_isBase)
			{
				_clk[i] = _clk[i + 1];
				_trp[i] = _trp[i + 1];
				_isb_GAL[i] = _isb_GAL[i + 1];
				_isb_BDS[i] = _isb_BDS[i + 1];
				_isb_GLO[i] = _isb_GLO[i + 1];
				_lost_isb_GAL[i] = _lost_isb_GAL[i + 1];
				_lost_isb_BDS[i] = _lost_isb_BDS[i + 1];
				_lost_isb_GLO[i] = _lost_isb_GLO[i + 1];
			}
		}
		if(_isBase)
		{
			_vDD_msg.erase(_vDD_msg.begin());
		}
		else
		{
			_vIF_msg.erase(_vIF_msg.begin());
		}
		delete _rover_window[gwindow_size - 1];
		_rover_window[gwindow_size - 1] = nullptr;
		// _amb_manager->slidingWindow();
		if (_isBase)
		{
			_amb_manager->slidingWindow();
		}
		else
		{
			_ambIF_manager->slidingWindow();
		}
		_rover_count--;
	}
	if (_isBase)
	{
		_amb_manager->get_last_epoch_sats(_epoch.sow());
	}
	else
	{
		_ambIF_manager->get_last_epoch_sats(_epoch.sow());
	}

}

void gfgomsf::t_gpvtfgo::_marginalization()
{

	// The residual factors added during marginalization include: the marginalization factor retained from the previous iteration, the absolute observation factor of the current dual-difference pseudorange and phase, and the initial ambiguity constraint factor.
	// full == gwindow_size -1
	if (_rover_count == gwindow_size-1)
	{
		_vector_to_double();
		ceres::LossFunction *loss_function;
		loss_function = new ceres::CauchyLoss(2);
		//loss_function = new  ceres::CauchyLoss(_loss_func_value);
		GNSSInfo *gnss_marginalization_info = new GNSSInfo();
		if (_last_gnss_marginalization_info && _last_gnss_marginalization_info->valid)
		{
			vector<int> drop_set;
			vector<int> amb_margin = _amb_manager->getMarginAmb();
			// preevious para_blocks
			for (int i = 0; i < static_cast<int>(_last_gnss_marginalization_para_blocks.size()); i++)
			{
				for (int j = 0; j < amb_margin.size(); j++)
				{
					if (_last_gnss_marginalization_para_blocks[i] == _para_amb[amb_margin[j]])
						drop_set.push_back(i);
				}
			}
			// construct new marginlization_factor
			MarginalizationGNSSFactor *marginalization_gnss_factor = new MarginalizationGNSSFactor(_last_gnss_marginalization_info);
			GNSSResidualBlockInfo *residual_block_info = new GNSSResidualBlockInfo(marginalization_gnss_factor, NULL,
				_last_gnss_marginalization_para_blocks,
				drop_set);
			gnss_marginalization_info->addResidualBlockInfo(residual_block_info);
		}
		/*
		for (int i = 0; i <= _rover_count; i++)
		{
			if (i == 0)
			{
				t_gallpar params_temp = _para_window[i];
				vector<DDEquMsg> dd_msg = _vDD_msg[i];
				vector<t_gsatdata> b_sat_data = _win_base_data[i];
				for (auto &dd_iter : dd_msg)
				{
					if (!_get_DD_data(dd_iter, b_sat_data)) continue;
					pair<string, string> base_rover_site = make_pair(dd_iter.base_site, dd_iter.rover_site);
					pair<FREQ_SEQ, GOBSBAND> freq_band = make_pair(dd_iter.freq, dd_iter.band);
					vector<pair<t_gsatdata, t_gsatdata>> DD_sat_data;
					DD_sat_data.push_back(make_pair(dd_iter.base_ref_sat, dd_iter.rover_ref_sat));
					DD_sat_data.push_back(make_pair(dd_iter.base_nonref_sat, dd_iter.rover_nonref_sat));
					GOBSTYPE  obstype = dd_iter.obs_type;
					if (obstype == GOBSTYPE::TYPE_C)
					{
						PseudorangeDDFactor *pf = new PseudorangeDDFactor(dd_iter.time, base_rover_site, params_temp, DD_sat_data, _gbias_model, freq_band);
						GNSSResidualBlockInfo *residual_block_info = new GNSSResidualBlockInfo(pf, NULL, vector<double *> {_para_CRD[i]}, vector<int>{0});
						gnss_marginalization_info->addResidualBlockInfo(residual_block_info);
					}
					if (obstype == GOBSTYPE::TYPE_L)
					{
						int id1, id2;
						id1 = _amb_manager->getAmbSearchIndex(make_pair(dd_iter.ref_sat_global_id, dd_iter.freq));
						id2 = _amb_manager->getAmbSearchIndex(make_pair(dd_iter.nonref_sat_global_id, dd_iter.freq));

						if (id1 == -1 || id2 == -1)
							continue;

						vector<int> drop_set{ 0 };
						if (_amb_manager->getAmbStartRoverID(id1) == 0 && _amb_manager->getAmbStartRoverID(id2) == 0)
						{
							if (_amb_manager->getAmbEndRoverID(id1) == 0)
								drop_set.push_back(1);
							if (_amb_manager->getAmbEndRoverID(id2) == 0)
								drop_set.push_back(2);

							CarrierphaseDDFactor *lf = new CarrierphaseDDFactor(dd_iter.time, base_rover_site, params_temp, DD_sat_data, _gbias_model, freq_band);
							GNSSResidualBlockInfo *residual_block_info = new GNSSResidualBlockInfo(lf, NULL, vector<double *> {_para_CRD[i], _para_amb[id1], _para_amb[id2]}, drop_set);
							gnss_marginalization_info->addResidualBlockInfo(residual_block_info);
						}

					}
				}
			}
		}	*/

		for (int i = 0; i <= _rover_count; i++)
		{
			if (i == 0)
			{
				t_gallpar params_temp = _para_window[i];
				for (auto &dd_iter : _vDD_msg[i])
				{
					pair<string, string> base_rover_site = make_pair(dd_iter.base_site, dd_iter.rover_site);
					pair<FREQ_SEQ, GOBSBAND> freq_band = make_pair(dd_iter.freq, dd_iter.band);
					vector<pair<t_gsatdata, t_gsatdata>> DD_sat_data;
					DD_sat_data.push_back(make_pair(dd_iter.base_ref_sat, dd_iter.rover_ref_sat));
					DD_sat_data.push_back(make_pair(dd_iter.base_nonref_sat, dd_iter.rover_nonref_sat));
					GOBSTYPE  obstype = dd_iter.obs_type;
					if (obstype == GOBSTYPE::TYPE_C)
					{
						PseudorangeDDFactor *pf = new PseudorangeDDFactor(dd_iter.time, base_rover_site, params_temp, DD_sat_data, _gbias_model, freq_band);
						GNSSResidualBlockInfo *residual_block_info = new GNSSResidualBlockInfo(pf, loss_function, vector<double *> {_para_CRD[i]}, vector<int>{0});
						gnss_marginalization_info->addResidualBlockInfo(residual_block_info);
					}
					if (obstype == GOBSTYPE::TYPE_L)
					{
						int id1, id2;
						id1 = _amb_manager->getAmbSearchIndex(make_pair(dd_iter.ref_sat_global_id, dd_iter.freq));
						id2 = _amb_manager->getAmbSearchIndex(make_pair(dd_iter.nonref_sat_global_id, dd_iter.freq));

						if (id1 == -1 || id2 == -1)
							continue;

						vector<int> drop_set{ 0 };
						if (_amb_manager->getAmbStartRoverID(id1) == 0 && _amb_manager->getAmbStartRoverID(id2) == 0)
						{
							if (_amb_manager->getAmbEndRoverID(id1) == 0)
								drop_set.push_back(1);
							if (_amb_manager->getAmbEndRoverID(id2) == 0)
								drop_set.push_back(2);

							CarrierphaseDDFactor *lf = new CarrierphaseDDFactor(dd_iter.time, base_rover_site, params_temp, DD_sat_data, _gbias_model, freq_band);
							GNSSResidualBlockInfo *residual_block_info = new GNSSResidualBlockInfo(lf, NULL, vector<double *> {_para_CRD[i], _para_amb[id1], _para_amb[id2]}, drop_set);
							gnss_marginalization_info->addResidualBlockInfo(residual_block_info);
						}

					}
				}
			}
		}


		gnss_marginalization_info->preMarginalize();
		gnss_marginalization_info->marginalize();


		std::unordered_map<long, double *> addr_shift;

		vector<int> cur_amb = _amb_manager->getCurWinAmb();
		for (int i = 0; i < cur_amb.size(); i++)
		{
			addr_shift[reinterpret_cast<long>(_para_amb[cur_amb[i]])] = _para_amb[cur_amb[i]];
		}

		vector<double *> parameter_blocks = gnss_marginalization_info->getParameterBlocks(addr_shift);
		if (_last_gnss_marginalization_info) delete _last_gnss_marginalization_info;
		// save marginalization info
		_last_gnss_marginalization_info = gnss_marginalization_info;
		_last_gnss_marginalization_para_blocks = parameter_blocks;
	}
}

void gfgomsf::t_gpvtfgo::_marginalization_PPP()
{
    //if (_rover_count == GWIN_SIZE)//delete for MultiWindow
    if (_rover_count == gwindow_size-1)//add for MultiWindow
    {
        cout << "###oldest_pos " << _headers[0] << std::fixed << std::setprecision(4) << " " << _Pos[0].coeff(0) << " " << _Pos[0].coeff(1) << " " << _Pos[0].coeff(2) << endl;
        _vector_to_double();//fgo2est
        ceres::LossFunction* loss_function;
        ceres::LossFunction* loss_function_CP;
        //loss_function =NULL;// new ceres::HuberLoss(_loss_func_value);
        //loss_function_CP = NULL;//new ceres::HuberLoss(_loss_func_value);
    	loss_function = new ceres::HuberLoss(_loss_func_value);
    	loss_function_CP = new ceres::HuberLoss(_loss_func_value);
        GNSSInfo* gnss_marginalization_info = new GNSSInfo();
        if (_last_gnss_marginalization_info && _last_gnss_marginalization_info->valid)
        {
            vector<int> drop_set;
            vector<int> amb_margin = _ambIF_manager->getMarginAmb();
            for (int i = 0; i < static_cast<int>(_last_gnss_marginalization_para_blocks.size()); i++)
            {
                //?????
                for (int j = 0; j < amb_margin.size(); j++)
                {
                    if (_last_gnss_marginalization_para_blocks[i] == _para_AMB_IF[amb_margin[j]])
                        drop_set.push_back(i);
                }
                if (_last_gnss_marginalization_para_blocks[i] == _para_TRP[0])
                {
                    drop_set.push_back(i);
                }
                if (_last_gnss_marginalization_para_blocks[i] == _para_ISB_GAL[0] && !_lost_isb_GAL[0])
                {
                    drop_set.push_back(i);
                }
                if (_last_gnss_marginalization_para_blocks[i] == _para_ISB_BDS[0] && !_lost_isb_BDS[0])
                {
                    drop_set.push_back(i);
                }
				if (_last_gnss_marginalization_para_blocks[i] == _para_ISB_GLO[0] && !_lost_isb_GLO[0])
				{
					drop_set.push_back(i);
				}
            }
            MarginalizationGNSSFactor* marginalization_gnss_factor = new MarginalizationGNSSFactor(_last_gnss_marginalization_info);
            GNSSResidualBlockInfo* residual_block_info = new GNSSResidualBlockInfo(marginalization_gnss_factor, NULL,
                _last_gnss_marginalization_para_blocks,
                drop_set);
            gnss_marginalization_info->addResidualBlockInfo(residual_block_info);
        }
        double s = 1;
        vector<int> drop_set{ 0 };
        const double graph_dt = std::fabs(_headers[1] - _headers[0]);
        double sqrt_info = 1.0 / sqrt(graph_interval_random_walk_q(_trpStoModel, graph_dt));
        RandomWalkFactor* rf = new RandomWalkFactor(s * sqrt_info);
        GNSSResidualBlockInfo* residual_block_info = new GNSSResidualBlockInfo(rf, NULL, vector<double*> {_para_TRP[0], _para_TRP[1]}, drop_set);
        gnss_marginalization_info->addResidualBlockInfo(residual_block_info);

        if (!_lost_isb_GAL[0] && !_lost_isb_GAL[1])
        {
            double sqrt_info_gal = 1.0 / sqrt(graph_interval_random_walk_q(_galStoModel, graph_dt));
            RandomWalkFactor* galf = new RandomWalkFactor(s * sqrt_info_gal);
            GNSSResidualBlockInfo* residual_block_info = new GNSSResidualBlockInfo(galf, NULL, vector<double*> {_para_ISB_GAL[0], _para_ISB_GAL[1]}, drop_set);
            gnss_marginalization_info->addResidualBlockInfo(residual_block_info);
        }
        if (!_lost_isb_BDS[0] && !_lost_isb_BDS[1])
        {
            double sqrt_info_bds = 1.0 / sqrt(graph_interval_random_walk_q(_bdsStoModel, graph_dt));
            RandomWalkFactor* bdsf = new RandomWalkFactor(s * sqrt_info_bds);
            GNSSResidualBlockInfo* residual_block_info = new GNSSResidualBlockInfo(bdsf, NULL, vector<double*> {_para_ISB_BDS[0], _para_ISB_BDS[1]}, drop_set);
            gnss_marginalization_info->addResidualBlockInfo(residual_block_info);
        }
		if (!_lost_isb_GLO[0] && !_lost_isb_GLO[1])
		{
			double sqrt_info_glo = 1.0 / sqrt(graph_interval_random_walk_q(_gloStoModel, graph_dt));
			RandomWalkFactor* glof = new RandomWalkFactor(s * sqrt_info_glo);
			GNSSResidualBlockInfo* residual_block_info = new GNSSResidualBlockInfo(glof, NULL, vector<double*> {_para_ISB_GLO[0], _para_ISB_GLO[1]}, drop_set);
			gnss_marginalization_info->addResidualBlockInfo(residual_block_info);
		}

        for (int i = 0; i <= _rover_count; i++)
        {
            if (i == 0)
            {
                t_gallpar params_temp = _para_window[i];
                vector<IFEquMsg> IF_tmp = _vIF_msg[i];
                for (auto& if_iter : IF_tmp)
                {
                    GOBSTYPE  obstype = if_iter.obs_type;
                    string sat_name = if_iter.satdata.sat();
                    GSYS gsys = if_iter.satdata.gsys();
                    if (obstype == GOBSTYPE::TYPE_C)
                    {
                        // 第0历元码因子的CRD/CLK/TRP都将离窗，三者全部消去。
                        vector<int> drop_set{ 0,1,2 };
                        pair<FREQ_SEQ, GOBSBAND> freq_band1 = make_pair(if_iter.freq1, if_iter.band1);
                        pair<FREQ_SEQ, GOBSBAND> freq_band2 = make_pair(if_iter.freq2, if_iter.band2);
                        if (gsys == GSYS::GPS)
                        {
                            PseudorangeIFFactor* pf = new PseudorangeIFFactor(if_iter.time, if_iter.site, params_temp, if_iter.satdata, _gbias_model, freq_band1, freq_band2);
                            GNSSResidualBlockInfo* residual_block_info = new GNSSResidualBlockInfo(pf, loss_function, vector<double*> {_para_CRD[i], _para_CLK[i], _para_TRP[i]}, drop_set);
                            gnss_marginalization_info->addResidualBlockInfo(residual_block_info);
                        }
                        else if (gsys == GSYS::GAL && !_lost_isb_GAL[i])
                        {
                            drop_set.push_back(3);
                            MultiPseudorangeIFFactor* pf = new MultiPseudorangeIFFactor(if_iter.time, if_iter.site, params_temp, if_iter.satdata, _gbias_model, freq_band1, freq_band2);
                            GNSSResidualBlockInfo* residual_block_info = new GNSSResidualBlockInfo(pf, loss_function, vector<double*> {_para_CRD[i], _para_CLK[i], _para_TRP[i], _para_ISB_GAL[i]}, drop_set);
                            gnss_marginalization_info->addResidualBlockInfo(residual_block_info);
                        }
                        else if (gsys == GSYS::BDS && !_lost_isb_BDS[i])
                        {
                            drop_set.push_back(3);
                            MultiPseudorangeIFFactor* pf = new MultiPseudorangeIFFactor(if_iter.time, if_iter.site, params_temp, if_iter.satdata, _gbias_model, freq_band1, freq_band2);
                            GNSSResidualBlockInfo* residual_block_info = new GNSSResidualBlockInfo(pf, loss_function, vector<double*> {_para_CRD[i], _para_CLK[i], _para_TRP[i], _para_ISB_BDS[i]}, drop_set);
                            gnss_marginalization_info->addResidualBlockInfo(residual_block_info);
                        }
						else if (gsys == GSYS::GLO && !_lost_isb_GLO[i])
						{
							drop_set.push_back(3);
							MultiPseudorangeIFFactor* pf = new MultiPseudorangeIFFactor(if_iter.time, if_iter.site, params_temp, if_iter.satdata, _gbias_model, freq_band1, freq_band2);
							GNSSResidualBlockInfo* residual_block_info = new GNSSResidualBlockInfo(pf, loss_function, vector<double*> {_para_CRD[i], _para_CLK[i], _para_TRP[i], _para_ISB_GLO[i]}, drop_set);
							gnss_marginalization_info->addResidualBlockInfo(residual_block_info);
						}
                        else
                        {
                            cerr << "not support the system now : " << gsys << endl;
                        }
                    }//TYPE_C
                    if (obstype == GOBSTYPE::TYPE_L)
                    {
                        vector<int> drop_set{ 0,1,2 };
                        int amb_id = _ambIF_manager->getAmbSearchIndex(if_iter.sat_id);
                        if (amb_id == -1)
                        {
                            continue;
                        }

                        if (_ambIF_manager->getAmbStartRoverID(amb_id) == 0)
                        {
                            bool drop_amb = false;
                            if (_ambIF_manager->getAmbEndRoverID(amb_id) == 0)
                            {
                                drop_amb = true;
                            }
                            pair<FREQ_SEQ, GOBSBAND> freq_band1 = make_pair(if_iter.freq1, if_iter.band1);
                            pair<FREQ_SEQ, GOBSBAND> freq_band2 = make_pair(if_iter.freq2, if_iter.band2);
                            if (gsys == GSYS::GPS)
                            {
                                if (drop_amb)
                                {
                                    drop_set.push_back(3);
                                }
                                CarrierphaseIFFactor* lf = new CarrierphaseIFFactor(if_iter.time, if_iter.site, params_temp, if_iter.satdata, _gbias_model, freq_band1, freq_band2);
                                GNSSResidualBlockInfo* residual_block_info = new GNSSResidualBlockInfo(lf, loss_function_CP, vector<double*> {_para_CRD[i], _para_CLK[i], _para_TRP[i], _para_AMB_IF[amb_id]}, drop_set);
                                gnss_marginalization_info->addResidualBlockInfo(residual_block_info);
                            }
                            else if (gsys == GSYS::GAL && !_lost_isb_GAL[i])
                            {
                                drop_set.push_back(3);//for ISB
                                if (drop_amb)
                                {
                                    drop_set.push_back(4);
                                }
                                MultiCarrierphaseIFFactor* lf = new MultiCarrierphaseIFFactor(if_iter.time, if_iter.site, params_temp, if_iter.satdata, _gbias_model, freq_band1, freq_band2);
                                GNSSResidualBlockInfo* residual_block_info = new GNSSResidualBlockInfo(lf, loss_function_CP, vector<double*> {_para_CRD[i], _para_CLK[i], _para_TRP[i], _para_ISB_GAL[i], _para_AMB_IF[amb_id]}, drop_set);
                                gnss_marginalization_info->addResidualBlockInfo(residual_block_info);
                            }
                            else if (gsys == GSYS::BDS && !_lost_isb_BDS[i])
                            {
                                drop_set.push_back(3);//for ISB
                                if (drop_amb)
                                {
                                    drop_set.push_back(4);
                                }
                                MultiCarrierphaseIFFactor* lf = new MultiCarrierphaseIFFactor(if_iter.time, if_iter.site, params_temp, if_iter.satdata, _gbias_model, freq_band1, freq_band2);
                                GNSSResidualBlockInfo* residual_block_info = new GNSSResidualBlockInfo(lf, loss_function_CP, vector<double*> {_para_CRD[i], _para_CLK[i], _para_TRP[i], _para_ISB_BDS[i], _para_AMB_IF[amb_id]}, drop_set);
                                gnss_marginalization_info->addResidualBlockInfo(residual_block_info);
                            }
							else if (gsys == GSYS::GLO && !_lost_isb_GLO[i])
							{
								drop_set.push_back(3);
								if (drop_amb)
									drop_set.push_back(4);
								MultiCarrierphaseIFFactor* lf = new MultiCarrierphaseIFFactor(if_iter.time, if_iter.site, params_temp, if_iter.satdata, _gbias_model, freq_band1, freq_band2);
								GNSSResidualBlockInfo* residual_block_info = new GNSSResidualBlockInfo(lf, loss_function_CP, vector<double*> {_para_CRD[i], _para_CLK[i], _para_TRP[i], _para_ISB_GLO[i], _para_AMB_IF[amb_id]}, drop_set);
								gnss_marginalization_info->addResidualBlockInfo(residual_block_info);
							}
                            else
                            {
                                cerr << "not support the system now : " << gsys << endl;
                            }
                        }//getAmbStartRoverID
                    }//TYPE_L
                }//if_iter
            }//i = 0
        }//_rover_count
        gnss_marginalization_info->preMarginalize();
        gnss_marginalization_info->marginalize();
        std::unordered_map<long, double*> addr_shift;

        vector<int> cur_amb = _ambIF_manager->getCurWinAmb();
        for (int i = 0; i < cur_amb.size(); i++)
        {
            addr_shift[reinterpret_cast<long>(_para_AMB_IF[cur_amb[i]])] = _para_AMB_IF[cur_amb[i]];
        }
        addr_shift[reinterpret_cast<long>(_para_TRP[1])] = _para_TRP[0];
        if (!_lost_isb_GAL[0] && !_lost_isb_GAL[1])
        {
            addr_shift[reinterpret_cast<long>(_para_ISB_GAL[1])] = _para_ISB_GAL[0];
        }
        if (!_lost_isb_BDS[0] && !_lost_isb_BDS[1])
        {
            addr_shift[reinterpret_cast<long>(_para_ISB_BDS[1])] = _para_ISB_BDS[0];
        }
		if (!_lost_isb_GLO[0] && !_lost_isb_GLO[1])
		{
			addr_shift[reinterpret_cast<long>(_para_ISB_GLO[1])] = _para_ISB_GLO[0];
		}

        vector<double*> parameter_blocks = gnss_marginalization_info->getParameterBlocks(addr_shift);
        if (_last_gnss_marginalization_info) delete _last_gnss_marginalization_info;
        _last_gnss_marginalization_info = gnss_marginalization_info;
        _last_gnss_marginalization_para_blocks = parameter_blocks;
    }//GWIN_SIZE
}

int gfgomsf::t_gpvtfgo::_gobs_outlier_detection(pair<string, int> & outlier)
{
	pair<string, int> sat_id;
	int idx = -1;
	if (_last_gnss_info->valid)
	{
		Eigen::VectorXd v = _last_gnss_info->v_norm;
		//std::cout << v << endl;
		int nobs = v.rows();
		double max = 0.0;

		for (int i = 0; i < nobs; i++)
		{
			if (fabs(v(i)) > max && fabs(v(i)) > _max_res_norm)
			{
				max = fabs(v(i));
				idx = i;
			}
		}
		if (idx >= 0)
		{
			sat_id = _gnss_obs_index[idx].first;
			outlier = sat_id;
			int id = sat_id.second;


			string obsType = gobstype2str(_gnss_obs_index[idx].second.second);
			int freq = static_cast<int>(_gnss_obs_index[idx].second.first);

			int cur_sats = -1;
			int cur_obs = -1;

			if (!_isBase)
			{
				cur_sats = static_cast<int>(_ambIF_manager->cur_sats.size());
				cur_obs = static_cast<int>(_vIF_msg[_rover_count].size());
			}
			else
			{
				cur_sats = static_cast<int>(_amb_manager->cur_sats.size());
				cur_obs = static_cast<int>(_vDD_msg[_rover_count].size());
			}

			for (int j = 0; j < nobs; ++j)
			{
				const auto& obs = _gnss_obs_index[j];
				const std::string& sat = obs.first.first;
				int sat_id = obs.first.second;
				int freq = static_cast<int>(obs.second.first);
				std::string type = gobstype2str(obs.second.second);

				// out << _epoch.str_ymdhms() << ","
				// 	<< _epoch.sow() + _epoch.dsec() << ","
				// 	<< _cntrep << ","
				// 	<< sat << ","
				// 	<< sat_id << ","
				// 	<< type << ","
				// 	<< freq << ","
				// 	<< _last_gnss_info->linearized_residuals(j) << ","
				// 	<< _last_gnss_info->postfit_qv_diag(j) << ","
				// 	<< _last_gnss_info->v_norm(j) << ","
				// 	<< _last_gnss_info->postfit_nobs << ","
				// 	<< _last_gnss_info->postfit_npar << ","
				// 	<< _last_gnss_info->postfit_dof_raw << ","
				// 	<< cur_sats << ","
				// 	<< (j == idx ? 1 : 0)
				// 	<< "\n";
			}

			auto it_find = find_if(_removed_sats.begin(), _removed_sats.end(), [id](pair<string, int> & sat_id)
			{
				return sat_id.second == id;

			});
			if (it_find == _removed_sats.end())
				_removed_sats.push_back(sat_id);

			//string obsType = gobstype2str(_gnss_obs_index[idx].second.second);
			ostringstream os;
			/*os << _site << " outlier (" << obsType << _obs_index[idx].second.first << ") " << sat
				<< " v: " << fixed << setw(16) << right << setprecision(3) << max;
			if (_log)
				_log->comment(1, "gpvtfgo", _epoch.str_ymdhms(" epoch ") + os.str());*/


			std::cout << _site << " outlier (" << obsType << _gnss_obs_index[idx].second.first << ") " << sat_id.first
				<< " v: " << fixed << setw(16) << right << setprecision(3) << max << endl;

			if (_spdlog)
				SPDLOG_LOGGER_ERROR(_spdlog, string("gpvtfgo "), _epoch.str_ymdhms(" epoch ") + os.str());
		}

		// if (_vDD_msg[_rover_count].size() - _removed_sats.size() <= 2)
		// 	idx = -1;
		if (_isBase)
		{
			if (_vDD_msg[_rover_count].size() - _removed_sats.size() <= 2)
				idx = -1;
		}
		else
		{
			 // (_vIF_msg[_rover_count].size() - _removed_sats.size() <= _minsat )
				// idx = -1;
		}
	}

	return idx;
}



bool gfgomsf::t_gpvtfgo::_check_outlier(const string & sat)
{
	auto it_find = find_if(outlier_sats.begin(), outlier_sats.end(), [sat](const string & sat_name)
	{
		return sat_name == sat;

	});
	if (it_find != outlier_sats.end())
		return true;
	else
		return false;
}

void gfgomsf::t_gpvtfgo::_posteriori_test(ceres::Problem& problem)
{
	cout << "begin to gnss posteriori_test" << endl;
	//for current gnss epoch
	_all_para_win.delAllParam();
	ceres::LossFunction* loss_function;
	//loss_function = new ceres::HuberLoss(_loss_func_value);
	loss_function = new ceres::CauchyLoss(_loss_func_value);
	GNSSInfo* gnss_info = new GNSSInfo();
	if (_vDD_msg[_rover_count].size() > 2)
	{
		//constrcut window para_index
		vector<vector<int>> crd_para_col_index;
		vector<vector<int>> amb_para_col_index;
		vector<double*> _parameter_blocks;
		int total_para_size = 0;
		vector<int> crdi;
		vector<par_type> partype{ par_type::CRD_X ,par_type::CRD_Y ,par_type::CRD_Z };
		_parameter_blocks.push_back(_para_CRD[_rover_count]);
		for (int j = 0; j < 3; j++)
		{
			t_gpar par_crd;
			par_crd.site = _site;
			par_crd.parType = partype[j];
			par_crd.value(_para_CRD[_rover_count][j]);
			par_crd.beg = _rover_window[_rover_count]->cur_time;
			par_crd.end = _rover_window[_rover_count]->cur_time;
			par_crd.index = j + 1;
			_all_para_win.addParam(par_crd);
			int id = j;
			crdi.push_back(id);
			total_para_size = total_para_size + 1;
		}
		crd_para_col_index.push_back(crdi);

		int amb_index_start = 3;
		map<int, int> amb_col_id;
		int amb_size = -1;
		_gnss_obs_index.clear();
		t_gallpar params_temp = _para_window[_rover_count];
		vector<DDEquMsg> dd_msg = _vDD_msg[_rover_count];
		vector<t_gsatdata> b_sat_data = _win_base_data[_rover_count];
		for (auto& dd_iter : dd_msg)
		{
			if (!_get_DD_data(dd_iter, b_sat_data)) continue;
			pair<string, string> base_rover_site = make_pair(dd_iter.base_site, dd_iter.rover_site);
			pair<FREQ_SEQ, GOBSBAND> freq_band = make_pair(dd_iter.freq, dd_iter.band);
			vector<pair<t_gsatdata, t_gsatdata>> DD_sat_data;
			DD_sat_data.push_back(make_pair(dd_iter.base_ref_sat, dd_iter.rover_ref_sat));
			DD_sat_data.push_back(make_pair(dd_iter.base_nonref_sat, dd_iter.rover_nonref_sat));
			GOBSTYPE  obstype = dd_iter.obs_type;
			_gnss_obs_index.push_back(make_pair(make_pair(dd_iter.rover_nonref_sat.sat(), dd_iter.nonref_sat_global_id), make_pair(dd_iter.freq, obstype)));
			if (obstype == GOBSTYPE::TYPE_C)
			{
				PseudorangeDDFactor* pf = new PseudorangeDDFactor(dd_iter.time, base_rover_site, params_temp, DD_sat_data, _gbias_model, freq_band);
				GNSSResidualBlockInfo* residual_block_info = new GNSSResidualBlockInfo(pf, NULL, vector<double*> {_para_CRD[_rover_count]});
				map<long, vector<int>> para_index;
				para_index[reinterpret_cast<long>(_para_CRD[_rover_count])] = crd_para_col_index[0];
				gnss_info->addResidualBlockInfo(residual_block_info, para_index);
			}
			if (obstype == GOBSTYPE::TYPE_L)
			{
				vector<int> amb_id12(2);
				amb_id12[0] = _amb_manager->getAmbSearchIndex(make_pair(dd_iter.ref_sat_global_id, dd_iter.freq));
				amb_id12[1] = _amb_manager->getAmbSearchIndex(make_pair(dd_iter.nonref_sat_global_id, dd_iter.freq));

				if (amb_id12[0] == -1 || amb_id12[1] == -1)
					continue;

				for (int i = 0; i < amb_id12.size(); i++)
				{
					int amb_id = amb_id12[i];
					if (amb_col_id.find(amb_id) == amb_col_id.end())
					{
						amb_size++;
						int amb_para_id = amb_index_start + amb_size;
						amb_col_id[amb_id] = amb_para_id;
						amb_para_col_index.push_back(vector<int> {amb_para_id});
						_amb_manager->addGpara(_all_para_win, amb_id, _para_amb[amb_id][0]);
						_parameter_blocks.push_back(_para_amb[amb_id]);
					}
				}
				map<long, vector<int>> para_index;
				//cout << "para_col_index: " << para_col_index.size() << " id1: " << _rover_count + 1+ id1 << " id2: " << _rover_count + 1+id2 << endl;
				double* addr1 = _para_amb[amb_id12[0]];
				double* addr2 = _para_amb[amb_id12[1]];
				para_index[reinterpret_cast<long>(_para_CRD[_rover_count])] = crd_para_col_index[0];
				para_index[reinterpret_cast<long>(addr1)] = vector<int>{ amb_col_id[amb_id12[0]] };
				para_index[reinterpret_cast<long>(addr2)] = vector<int>{ amb_col_id[amb_id12[1]] };
				CarrierphaseDDFactor* lf = new CarrierphaseDDFactor(dd_iter.time, base_rover_site, params_temp, DD_sat_data, _gbias_model, freq_band);
				GNSSResidualBlockInfo* residual_block_info = new GNSSResidualBlockInfo(lf, NULL, vector<double*> {_para_CRD[_rover_count], _para_amb[amb_id12[0]], _para_amb[amb_id12[1]]});
				gnss_info->addResidualBlockInfo(residual_block_info, para_index);
			}
		}
		total_para_size = total_para_size + amb_size + 1;
		//construct variances by ceres solver
		ceres::Covariance::Options options_co;
		//options_co.algorithm_type = ceres::SPARSE_QR;
		//options_co.algorithm_type = ceres::DENSE_SVD;
		//options_co.sparse_linear_algebra_library_type = ceres::SparseLinearAlgebraLibraryType::SUITE_SPARSE;
		//options_co.apply_loss_function = false; //optional, true or false is depended on the reliability of covariance
		ceres::Covariance covariance(options_co);
		std::vector<const double*> covariance_blocks; //all parameter_blocks
		for (int i = 0; i < _parameter_blocks.size(); i++)
		{
			covariance_blocks.push_back(_parameter_blocks[i]);
		}
		Eigen::MatrixXd Qx = Eigen::MatrixXd::Zero(total_para_size, total_para_size);
		try
		{
			covariance.Compute(covariance_blocks, &problem);
			covariance.GetCovarianceMatrix(covariance_blocks, Qx.data());
		}
		catch (...)
		{
			cout << "Covariance Compute Failed" << endl;
		}
		gnss_info->constructEqu_fromCeres(Qx);
	}
	if (_last_gnss_info) delete _last_gnss_info;
	_last_gnss_info = gnss_info;

}





void gfgomsf::t_gpvtfgo::_posteriori_test_PPP(ceres::Problem& problem)
{
    //for current gnss epoch
    _all_para_win.delAllParam();
    ceres::LossFunction* loss_function;
    ceres::LossFunction* loss_function_CP;
    loss_function = new ceres::HuberLoss(_loss_func_value);
    loss_function = NULL;
    loss_function_CP = NULL;
    GNSSInfo* gnss_info = new GNSSInfo();

    //if (cur_sat_prn.size() < _minsat)//delete for MultiWindow
    if (_ambIF_manager->cur_sats.size() < _minsat)//add for MultiWindow
    {
        cout << "sat num : " << _ambIF_manager->cur_sats.size() << " < min sat num : " << _minsat << endl;
    }
    //constrcut window para_index
    vector<vector<int>> crd_para_col_index;
    vector<vector<int>> clk_para_col_index;
    vector<vector<int>> trp_para_col_index;
    vector<vector<int>> isb_gal_para_col_index;
    vector<vector<int>> isb_bds_para_col_index;
	vector<vector<int>> isb_glo_para_col_index;
    vector<vector<int>> amb_para_col_index;
    int total_para_size = 0;

    vector<int> crdi;
    vector<par_type> partype{ par_type::CRD_X ,par_type::CRD_Y ,par_type::CRD_Z };
    _parameter_blocks.clear();
    for (int j = 0; j < 3; j++)
    {
        t_gpar par_crd;
        par_crd.site = _site;
        par_crd.parType = partype[j];
        par_crd.value(_para_CRD[_rover_count][j]);
        //par_crd.beg = _rover_window[_rover_count]->cur_time;//delete for MultiWindow
        //par_crd.end = _rover_window[_rover_count]->cur_time;//delete for MultiWindow
        par_crd.beg = _epoch;//add for MultiWindow
        par_crd.end = _epoch;//add for MultiWindow
        par_crd.index = j + 1;
        _all_para_win.addParam(par_crd);
        int id = j;
        crdi.push_back(id);
        total_para_size = total_para_size + 1;
    }
    crd_para_col_index.push_back(crdi);
    _parameter_blocks.push_back(_para_CRD[_rover_count]);

    vector<int> clki;
    t_gpar par_clk;
    par_clk.site = _site;
    par_clk.parType = par_type::CLK;
    par_clk.value(_para_CLK[_rover_count][0]);
    //par_clk.beg = _rover_window[_rover_count]->cur_time;//delete for MultiWindow
    //par_clk.end = _rover_window[_rover_count]->cur_time;//delete for MultiWindow
    par_clk.beg = _epoch;//add for MultiWindow
    par_clk.end = _epoch;//add for MultiWindow
    par_clk.index = 4;
    _all_para_win.addParam(par_clk);
    clki.push_back(3);
    total_para_size = total_para_size + 1;
    clk_para_col_index.push_back(clki);
    _parameter_blocks.push_back(_para_CLK[_rover_count]);

    vector<int> trpi;
    t_gpar par_trp;
    par_trp.site = _site;
    par_trp.parType = par_type::TRP;
    par_trp.value(_para_TRP[_rover_count][0]);
    //par_trp.beg = _rover_window[_rover_count]->cur_time;//delete for MultiWindow
    //par_trp.end = _rover_window[_rover_count]->cur_time;//delete for MultiWindow
    par_trp.beg = _epoch;//add for MultiWindow
    par_trp.end = _epoch;//add for MultiWindow
    par_trp.index = 5;
    _all_para_win.addParam(par_trp);
    trpi.push_back(4);
    total_para_size = total_para_size + 1;
    trp_para_col_index.push_back(trpi);
    _parameter_blocks.push_back(_para_TRP[_rover_count]);

    int isb_num = 0;
    if (!_lost_isb_GAL[_rover_count])
    {
        isb_num++;
        vector<int> isb_gal_i;
        t_gpar par_isb_gal;
        par_isb_gal.site = _site;
        par_isb_gal.parType = par_type::GAL_ISB;
        par_isb_gal.value(_para_ISB_GAL[_rover_count][0]);
        //par_isb_gal.beg = _rover_window[_rover_count]->cur_time;//delete for MultiWindow
        //par_isb_gal.end = _rover_window[_rover_count]->cur_time;//delete for MultiWindow
        par_isb_gal.beg = _epoch;//add for MultiWindow
        par_isb_gal.end = _epoch;//add for MultiWindow
        par_isb_gal.index = 5 + isb_num;
        _all_para_win.addParam(par_isb_gal);
        isb_gal_i.push_back(4 + isb_num);
        total_para_size = total_para_size + 1;
        isb_gal_para_col_index.push_back(isb_gal_i);
        _parameter_blocks.push_back(_para_ISB_GAL[_rover_count]);
    }
    if (!_lost_isb_BDS[_rover_count])
    {
        isb_num++;
        vector<int> isb_bds_i;
        t_gpar par_isb_bds;
        par_isb_bds.site = _site;
        par_isb_bds.parType = par_type::BDS_ISB;
        par_isb_bds.value(_para_ISB_BDS[_rover_count][0]);
        //par_isb_bds.beg = _rover_window[_rover_count]->cur_time;//delete for MultiWindow
        //par_isb_bds.end = _rover_window[_rover_count]->cur_time;//delete for MultiWindow
        par_isb_bds.beg = _epoch;//add for MultiWindow
        par_isb_bds.end = _epoch;//add for MultiWindow
        par_isb_bds.index = 5 + isb_num;
        _all_para_win.addParam(par_isb_bds);
        isb_bds_i.push_back(4 + isb_num);
        total_para_size = total_para_size + 1;
        isb_bds_para_col_index.push_back(isb_bds_i);
        _parameter_blocks.push_back(_para_ISB_BDS[_rover_count]);
    }
	if (!_lost_isb_GLO[_rover_count])
	{
		isb_num++;
		vector<int> isb_glo_i;
		t_gpar par_isb_glo;
		par_isb_glo.site = _site;
		par_isb_glo.parType = par_type::GLO_ISB;
		par_isb_glo.value(_para_ISB_GLO[_rover_count][0]);
		par_isb_glo.beg = _epoch;
		par_isb_glo.end = _epoch;
		par_isb_glo.index = 5 + isb_num;
		_all_para_win.addParam(par_isb_glo);
		isb_glo_i.push_back(4 + isb_num);
		total_para_size++;
		isb_glo_para_col_index.push_back(isb_glo_i);
		_parameter_blocks.push_back(_para_ISB_GLO[_rover_count]);
	}

    int amb_index_start = 5 + isb_num;
    map<int, int> amb_col_id;
    int amb_size = -1;
    _gnss_obs_index.clear();
    t_gallpar params_temp = _para_window[_rover_count];
    vector<IFEquMsg> IF_tmp = _vIF_msg[_rover_count];
    for (auto& if_iter : IF_tmp)
    {
        GOBSTYPE  obstype = if_iter.obs_type;
        string sat_name = if_iter.satdata.sat();
        GSYS gsys = if_iter.satdata.gsys();
        _gnss_obs_index.push_back(make_pair(make_pair(if_iter.satdata.sat(), if_iter.sat_id), make_pair(if_iter.freq1, obstype)));//freq2 is not uesd here
        if (obstype == GOBSTYPE::TYPE_C)
        {
            map<long, vector<int>> para_index;
            para_index[reinterpret_cast<long>(_para_CRD[_rover_count])] = crd_para_col_index[0];
            para_index[reinterpret_cast<long>(_para_CLK[_rover_count])] = clk_para_col_index[0];
            para_index[reinterpret_cast<long>(_para_TRP[_rover_count])] = trp_para_col_index[0];
            pair<FREQ_SEQ, GOBSBAND> freq_band1 = make_pair(if_iter.freq1, if_iter.band1);
            pair<FREQ_SEQ, GOBSBAND> freq_band2 = make_pair(if_iter.freq2, if_iter.band2);
            if (gsys == GSYS::GPS)
            {
                PseudorangeIFFactor* pf = new PseudorangeIFFactor(if_iter.time, if_iter.site, params_temp, if_iter.satdata, _gbias_model, freq_band1, freq_band2);
                GNSSResidualBlockInfo* residual_block_info = new GNSSResidualBlockInfo(pf, loss_function, vector<double*> {_para_CRD[_rover_count], _para_CLK[_rover_count], _para_TRP[_rover_count]});
                gnss_info->addResidualBlockInfo(residual_block_info, para_index);
            }
            else if (gsys == GSYS::GAL && !_lost_isb_GAL[_rover_count])
            {
                MultiPseudorangeIFFactor* pf = new MultiPseudorangeIFFactor(if_iter.time, if_iter.site, params_temp, if_iter.satdata, _gbias_model, freq_band1, freq_band2);
                GNSSResidualBlockInfo* residual_block_info = new GNSSResidualBlockInfo(pf, loss_function, vector<double*> {_para_CRD[_rover_count], _para_CLK[_rover_count], _para_TRP[_rover_count], _para_ISB_GAL[_rover_count]});
                para_index[reinterpret_cast<long>(_para_ISB_GAL[_rover_count])] = isb_gal_para_col_index[0];
                gnss_info->addResidualBlockInfo(residual_block_info, para_index);
            }
            else if (gsys == GSYS::BDS && !_lost_isb_BDS[_rover_count])
            {
                MultiPseudorangeIFFactor* pf = new MultiPseudorangeIFFactor(if_iter.time, if_iter.site, params_temp, if_iter.satdata, _gbias_model, freq_band1, freq_band2);
                GNSSResidualBlockInfo* residual_block_info = new GNSSResidualBlockInfo(pf, loss_function, vector<double*> {_para_CRD[_rover_count], _para_CLK[_rover_count], _para_TRP[_rover_count], _para_ISB_BDS[_rover_count]});
                para_index[reinterpret_cast<long>(_para_ISB_BDS[_rover_count])] = isb_bds_para_col_index[0];
                gnss_info->addResidualBlockInfo(residual_block_info, para_index);
            }
			else if (gsys == GSYS::GLO && !_lost_isb_GLO[_rover_count])
			{
				MultiPseudorangeIFFactor* pf = new MultiPseudorangeIFFactor(if_iter.time, if_iter.site, params_temp, if_iter.satdata, _gbias_model, freq_band1, freq_band2);
				GNSSResidualBlockInfo* residual_block_info = new GNSSResidualBlockInfo(pf, loss_function, vector<double*> {_para_CRD[_rover_count], _para_CLK[_rover_count], _para_TRP[_rover_count], _para_ISB_GLO[_rover_count]});
				para_index[reinterpret_cast<long>(_para_ISB_GLO[_rover_count])] = isb_glo_para_col_index[0];
				gnss_info->addResidualBlockInfo(residual_block_info, para_index);
			}
            else
            {
                cerr << "not support the system now : " << gsys << endl;
            }
        }
        if (obstype == GOBSTYPE::TYPE_L)
        {
            int amb_id = _ambIF_manager->getAmbSearchIndex(if_iter.sat_id);
            if (amb_id == -1)
            {
                continue;
            }

            if (amb_col_id.find(amb_id) == amb_col_id.end())
            {
                amb_size++;
                int amb_para_id = amb_index_start + amb_size;
                amb_col_id[amb_id] = amb_para_id;
                amb_para_col_index.push_back(vector<int> {amb_para_id});
                _ambIF_manager->addGpara(_all_para_win, amb_id, _para_AMB_IF[amb_id][0]);
                _parameter_blocks.push_back(_para_AMB_IF[amb_id]);
            }
            map<long, vector<int>> para_index;
            double* addr = _para_AMB_IF[amb_id];
            para_index[reinterpret_cast<long>(_para_CRD[_rover_count])] = crd_para_col_index[0];
            para_index[reinterpret_cast<long>(_para_CLK[_rover_count])] = clk_para_col_index[0];
            para_index[reinterpret_cast<long>(_para_TRP[_rover_count])] = trp_para_col_index[0];
            para_index[reinterpret_cast<long>(addr)] = vector<int>{ amb_col_id[amb_id] };
            pair<FREQ_SEQ, GOBSBAND> freq_band1 = make_pair(if_iter.freq1, if_iter.band1);
            pair<FREQ_SEQ, GOBSBAND> freq_band2 = make_pair(if_iter.freq2, if_iter.band2);
            if (gsys == GSYS::GPS)
            {
                CarrierphaseIFFactor* lf = new CarrierphaseIFFactor(if_iter.time, if_iter.site, params_temp, if_iter.satdata, _gbias_model, freq_band1, freq_band2);
                GNSSResidualBlockInfo* residual_block_info = new GNSSResidualBlockInfo(lf, loss_function_CP, vector<double*> {_para_CRD[_rover_count], _para_CLK[_rover_count], _para_TRP[_rover_count], _para_AMB_IF[amb_id]});
                gnss_info->addResidualBlockInfo(residual_block_info, para_index);
            }
            else if (gsys == GSYS::GAL && !_lost_isb_GAL[_rover_count])
            {
                MultiCarrierphaseIFFactor* lf = new MultiCarrierphaseIFFactor(if_iter.time, if_iter.site, params_temp, if_iter.satdata, _gbias_model, freq_band1, freq_band2);
                GNSSResidualBlockInfo* residual_block_info = new GNSSResidualBlockInfo(lf, loss_function_CP, vector<double*> {_para_CRD[_rover_count], _para_CLK[_rover_count], _para_TRP[_rover_count], _para_ISB_GAL[_rover_count], _para_AMB_IF[amb_id]});
                para_index[reinterpret_cast<long>(_para_ISB_GAL[_rover_count])] = isb_gal_para_col_index[0];
                gnss_info->addResidualBlockInfo(residual_block_info, para_index);
            }
            else if (gsys == GSYS::BDS && !_lost_isb_BDS[_rover_count])
            {
                MultiCarrierphaseIFFactor* lf = new MultiCarrierphaseIFFactor(if_iter.time, if_iter.site, params_temp, if_iter.satdata, _gbias_model, freq_band1, freq_band2);
                GNSSResidualBlockInfo* residual_block_info = new GNSSResidualBlockInfo(lf, loss_function_CP, vector<double*> {_para_CRD[_rover_count], _para_CLK[_rover_count], _para_TRP[_rover_count], _para_ISB_BDS[_rover_count], _para_AMB_IF[amb_id]});
                para_index[reinterpret_cast<long>(_para_ISB_BDS[_rover_count])] = isb_bds_para_col_index[0];
                gnss_info->addResidualBlockInfo(residual_block_info, para_index);
            }
			else if (gsys == GSYS::GLO && !_lost_isb_GLO[_rover_count])
			{
				MultiCarrierphaseIFFactor* lf = new MultiCarrierphaseIFFactor(if_iter.time, if_iter.site, params_temp, if_iter.satdata, _gbias_model, freq_band1, freq_band2);
				GNSSResidualBlockInfo* residual_block_info = new GNSSResidualBlockInfo(lf, loss_function_CP, vector<double*> {_para_CRD[_rover_count], _para_CLK[_rover_count], _para_TRP[_rover_count], _para_ISB_GLO[_rover_count], _para_AMB_IF[amb_id]});
				para_index[reinterpret_cast<long>(_para_ISB_GLO[_rover_count])] = isb_glo_para_col_index[0];
				gnss_info->addResidualBlockInfo(residual_block_info, para_index);
			}
            else
            {
                cerr << "not support the system now : " << gsys << endl;
            }
        }
    }
    //construct covariance
    total_para_size = total_para_size + amb_size + 1;
    //construct variances by ceres solver
    ceres::Covariance::Options options_co;
    options_co.algorithm_type = ceres::SPARSE_QR;
    //options_co.algorithm_type = ceres::DENSE_SVD;
    options_co.apply_loss_function = false; //optional, true or false is depended on the reliability of covariance
    ceres::Covariance covariance(options_co);
    std::vector<const double*> covariance_blocks; //all parameter_blocks
    for (int i = 0; i < _parameter_blocks.size(); i++)
    {
        covariance_blocks.push_back(_parameter_blocks[i]);
    }
    //CHECK(covariance.Compute(covariance_blocks, &problem));
    bool cmp = covariance.Compute(covariance_blocks, &problem);
    Eigen::MatrixXd Qx(total_para_size, total_para_size);
    Qx.setZero();
    if (cmp == false)
    {
        for (int i = 0; i < crd_para_col_index.size(); i++)
        {
            for (int j = 0; j < 3; j++)
            {
                int id = crd_para_col_index[i][j];
                Qx(id, id) = _sig_init_crd * _sig_init_crd;
            }
        }
        for (int i = 0; i < clk_para_col_index.size(); i++)
        {
            int id = clk_para_col_index[i][0];
            Qx(id, id) = _clkStoModel->getQ() * _clkStoModel->getQ();
        }
        for (int i = 0; i < trp_para_col_index.size(); i++)
        {
            int id = trp_para_col_index[i][0];
            Qx(id, id) = _sig_init_ztd * _sig_init_ztd;
        }
        for (int i = 0; i < isb_gal_para_col_index.size(); i++)
        {
            int id = isb_gal_para_col_index[i][0];
            Qx(id, id) = _sig_init_gal * _sig_init_gal;
        }
        for (int i = 0; i < isb_bds_para_col_index.size(); i++)
        {
            int id = isb_bds_para_col_index[i][0];
            Qx(id, id) = _sig_init_bds * _sig_init_bds;
        }
		for (int i = 0; i < isb_glo_para_col_index.size(); i++)
		{
			int id = isb_glo_para_col_index[i][0];
			Qx(id, id) = _sig_init_glo * _sig_init_glo;
		}
        for (int i = 0; i < amb_para_col_index.size(); i++)
        {
            int id = amb_para_col_index[i][0];
            Qx(id, id) = _sigAmbig * _sigAmbig;
        }
        gnss_info->constructEqu(Qx, 0);
    }

    else
    {
        covariance.GetCovarianceMatrix(covariance_blocks, Qx.data());
        gnss_info->constructEqu(Qx, 1);
    }

    if (_last_gnss_info) delete _last_gnss_info;
    _last_gnss_info = gnss_info;
}

gfgomsf::t_gpvtfgo::DDEquMsg::DDEquMsg(const t_gsatdata& _ref_sat, const t_gsatdata& _nonref_sat, const GOBSTYPE& _obs_type, const FREQ_SEQ& _freq)
	: rover_ref_sat(_ref_sat),
	rover_nonref_sat(_nonref_sat),
	obs_type(_obs_type),
	freq(_freq),
	ref_sat(_ref_sat.sat()),
	nonref_sat(_nonref_sat.sat())
{
	this->time = _nonref_sat.epoch();
}

