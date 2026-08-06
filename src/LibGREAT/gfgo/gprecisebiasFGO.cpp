#include"gprecisebiasFGO.h"
//#include "gproc/glsqmatrix.h"

namespace gfgo
{

	t_gprecisebiasFGO::t_gprecisebiasFGO(t_gallproc *data, t_spdlog spdlog, t_gsetbase *setting) :t_gprecisebiasGPP(data, spdlog, setting)
	{
	}
	t_gprecisebiasFGO::~t_gprecisebiasFGO()
	{
	}

	bool t_gprecisebiasFGO::cmb_equ(bool isFGO, bool calculate_equ, t_gtime &epoch, t_gallpar &params, t_gsatdata &obsdata, t_gobs &gobs, t_gbaseEquation &result)
	{
		// Serialize residual evaluation on the shared scratch state when the
		// Ceres solve uses num_threads > 1.  The FGO factors call cmb_equ from
		// parallel worker threads and _prepare_obs_GPP_FGO/_omc_obs_ALL/...
		// mutate member buffers (_crt_obs, _crs_sat_crd, _crs_rec_crd, ...).
		std::lock_guard<std::mutex> lock(_cmb_equ_mutex);

		//cout << "call : cmb_equ" << endl;		
		// check obs_type valid
		double Obs_value = obsdata.getobs(gobs.gobs());
		// double snr_value = obsdata.getobs(pl2snr(gobs.gobs()));
		if (double_eq(Obs_value, 0.0))
		{
			if (_spdlog)
				SPDLOG_LOGGER_ERROR(_spdlog, "Obs_value is 0.0");
			return false;
		}

		tuple<string, string, t_gtime>* flag = &_rec_sat_before;
		if (make_tuple(obsdata.site(), obsdata.sat(), epoch) != *flag)
		{
			bool update_valid = t_gprecisebiasGPP::_update_obs_info_GPP(epoch, _gall_nav, _gallobj, obsdata, params);
			if (!update_valid)
			{
				if (_spdlog)
					SPDLOG_LOGGER_ERROR(_spdlog, "update obs information failed" + epoch.str_ymdhms("", false));
				return false;
			}

			/*cout << fixed << setprecision(8)
				<< ">> " << _crt_obs.site() << " " << _crt_obs.sat() << " " << _crt_obs.epoch().sow();
			cout << " " << gobs2str(gobs.gobs());

			if (_crt_obs.is_process == true)
			{
				cout << "\n TRUE " << endl;
			}
			else
			{
				cout << "\n FALSE " << endl;
			}*/

			// prepare Caculate some common bias[sat_pos,rec_pos,relative,rho]
			bool pre_valid;
			if (isFGO)
			{
				pre_valid = t_gprecisebiasFGO::_prepare_obs_GPP_FGO(epoch, _gall_nav, _gallobj, params);
			}
			else
			{
				pre_valid = t_gprecisebiasGPP::_prepare_obs_GPP(epoch, _gall_nav, _gallobj, params);
			}
			//bool pre_valid = t_gprecisebiasGPP::_prepare_obs_GPP(epoch, _gall_nav, _gallobj, params, obsdata);
			if (!pre_valid)
			{
				if (_spdlog)
					SPDLOG_LOGGER_ERROR(_spdlog, "prepare obs information failed" + epoch.str_ymdhms("", false));
				return false;
			}
			*flag = make_tuple(obsdata.site(), obsdata.sat(), epoch);
		}
		//cout << obsdata.site() << " " << obsdata.sat() << endl;

		if (calculate_equ)
		{
			// combine equ
			//t_glsqEquationMatrix equ;
			double omc = 0.0, wgt = 0.0;
			vector<pair<int, double>> coef;

			if (!_omc_obs_ALL(epoch, _crt_obs, params, gobs, omc))
			{
				if (_spdlog)
					SPDLOG_LOGGER_ERROR(_spdlog, "omc obs failed");
				return false;
			};

			if (!_wgt_obs_ALL(t_gdata::REC, gobs, _crt_obs, 1.0, wgt))
			{
				if (_spdlog)
					SPDLOG_LOGGER_ERROR(_spdlog, "weight obs failed");
				return false;
			};

			if (!_prt_obs_ALL(epoch, _crt_obs, params, gobs, coef))
			{
				if (_spdlog)
					SPDLOG_LOGGER_ERROR(_spdlog, "partialrange obs failed");
				return false;
			}

			static std::set<std::string> logged_equations;
			std::string diag_key = std::to_string(epoch.sow()) + "," + obsdata.sat() + "," + gobs2str(gobs.gobs());
			if (epoch.sow() <= 31741.5 && logged_equations.insert(diag_key).second)
			{
				std::ofstream diag("ppp_first_epoch_equations.csv", std::ios::app);
				diag << std::setprecision(15) << "OPEN," << epoch.sow() << "," << obsdata.sat()
				     << "," << gobs2str(gobs.gobs()) << "," << Obs_value << "," << omc << "," << wgt
				     << "," << _crt_obs.rho() << "," << _crt_obs.clk() << "," << _crt_obs.ele() << "\n";
			}

			/*cout << " omc= " << omc << " wgt= " << wgt << endl;
			for (int i = 0; i < coef.size(); i++)
			{
				cout << " coef= " << "( " << coef[i].first << " , " << coef[i].second << " )" << endl;
			}*/

			result.B.push_back(coef);
			result.P.push_back(wgt);
			result.l.push_back(omc);
		}
		_update_obs_info(obsdata);
		//obsdata.is_process = true;

		return true;
	}

	bool t_gprecisebiasFGO::_prepare_obs_GPP_FGO(const t_gtime &crt_epo, t_gallnav *gallnav, t_gallobj *gallobj, t_gallpar &pars)
	{
		//cout << "call : _prepare_obs_GPP" << endl;
		if (!gallnav || !gallobj)
		{
			if (_spdlog)
				SPDLOG_LOGGER_ERROR(_spdlog, "no navgation data or atx data for epoch " + crt_epo.str_ymdhms());
			return false;
		}

		// compute reciver time
		_crt_rec_epo = crt_epo - _crt_rec_clk;
		_crt_obs.addrecTime(_crt_rec_epo);

		// get Rec crd
		bool calculate_tides = false;// Correct tides or not
		bool apply_obj_valid = t_gprecisebias::_apply_rec_RTK(crt_epo, _crt_rec_epo, pars, calculate_tides);
		if (!apply_obj_valid)
		{
			if (_spdlog)
				SPDLOG_LOGGER_ERROR(_spdlog, "can not apply site in " + crt_epo.str_ymdhms());
			return false;
		}

		// get sat crd
		bool apply_sat_valid = _crt_obs.is_process;
		if (!_crt_obs.is_process)
		{
			apply_sat_valid = _apply_sat(_crt_rec_epo, _crt_sat_epo, gallnav);
			_crt_obs.is_process = true;
		}
		else
		{
			_crs_sat_crd = _crt_obs.satcrdcrs();
			_crs_sat_vel = _crt_obs.satvel_crs();
			_crt_sat_epo = _crt_obs.satTime();
			// addrho
			double tmp = (_crs_sat_crd - _crs_rec_crd).norm();
			_crt_obs.addrho(tmp);
			return true;
		}


		if (!apply_sat_valid)
		{
			if (_spdlog)
				SPDLOG_LOGGER_ERROR(_spdlog, "can not apply sat in " + crt_epo.str_ymdhms());
			return false;
		}
		_crt_obs.addsatTime(_crt_sat_epo);

#ifdef DEBUG
		cout << "site crs" << scientific << setw(25) << setprecision(15) << _crdSiteCrs[0] << "  " << _crdSiteCrs[1] << "  " << _crdSiteCrs[2] << "  " << endl;
		cout << "sat crs" << scientific << setw(25) << setprecision(15) << _crdSatCrs[0] << "  " << _crdSatCrs[1] << "  " << _crdSatCrs[2] << "  " << endl;
#endif // DEBUG21

		// get f01 PCO
		this->_crs_rec_pco = this->_crs_rec_crd;
		this->_crs_sat_pco = this->_crs_sat_crd;

		shared_ptr<t_gobj> sat_obj = this->_gallobj->obj(_crt_sat);
		shared_ptr<t_gobj> rec_obj = this->_gallobj->obj(_crt_rec);

		shared_ptr<t_gpcv> sat_pcv = (sat_obj != 0) ? sat_obj->pcv(crt_epo) : nullptr;
		shared_ptr<t_gpcv> rec_pcv = (rec_obj != 0) ? rec_obj->pcv(crt_epo) : nullptr;

		//zzwu changes according to lvhb
		if (!_isCalSatPCO)
		{
			sat_pcv = nullptr;
		}

		if (sat_pcv ) //add _realtime by xiongyun
		{
			// Satellite phase center offset
			t_gtriple pco(0, 0, 0);
			if (sat_pcv->pcoS(_crt_obs, pco, _observ, _band_index[_crt_sys][FREQ_1], _band_index[_crt_sys][FREQ_2]) > 0)
			{
				//Matrix _rot_matrix = _RotMatrix_Ant(_crt_obs, _crt_sat_epo, sat_obj, true);
				Matrix _rot_matrix = _RotMatrix_Ant(_crt_obs, _crt_epo, _crt_sat_epo, sat_obj, true);
				this->_crs_sat_pco += t_gtriple(_rot_matrix * pco.crd_cvect());
			}
		}

		if (rec_pcv)
		{
			// Receiver phase center offset
			t_gtriple pco(0.0, 0.0, 0.0);
			if (rec_pcv->pcoR(_crt_obs, pco, _observ, _band_index[_crt_sys][FREQ_1], _band_index[_crt_sys][FREQ_2]) > 0)
			{
				//Matrix _rot_matrix = _RotMatrix_Ant(_crt_obs, _crt_rec_epo, rec_obj, true);
				Matrix _rot_matrix = _RotMatrix_Ant(_crt_obs, _crt_epo, _crt_rec_epo, rec_obj, true);
				this->_crs_rec_pco += t_gtriple(_rot_matrix * pco.crd_cvect());
			}
		}

		// compute satclk[s]
		double sat_clk = _crt_sat_clk;

		// compute reldelay[m]
		double reldelay = relDelay(this->_crs_rec_pco, this->_crs_rec_vel, this->_crs_sat_pco, this->_crs_sat_vel);
		_crt_obs.addclk(sat_clk * CLIGHT - reldelay); // include reldelay
		_crt_obs.addreldelay(reldelay);                  // for debug

		// addrho
		double tmp = (_crs_sat_crd - _crs_rec_crd).norm();
		_crt_obs.addrho(tmp);

		// add drate
		_crt_obs.adddrate((DotProduct((_crs_sat_vel - _crs_rec_vel).crd_cvect(), (_crs_sat_pco - _crs_rec_pco).crd_cvect())) / (CLIGHT * tmp)); //yjqin

		// add azim && elev
		t_gtriple xyz_rho = _crs_sat_pco - _crs_rec_pco;
		t_gtriple ell_r, neu_s;

		Matrix _rot_ant2crs = _RotMatrix_Ant(_crt_obs, _crt_rec_epo, _crt_rec_epo, _crt_obj, true);
		neu_s = t_gtriple(_rot_ant2crs.t() * xyz_rho.crd_cvect());
		double NE2 = neu_s[0] * neu_s[0] + neu_s[1] * neu_s[1];
		double ele = acos(sqrt(NE2) / _crt_obs.rho());
		if (sqrt(NE2) / _crt_obs.rho() > 1.0)
		{
			_crt_obs.addele(0.0);
		}
		else
		{
			_crt_obs.addele(ele);
		}
		// for off-nadir angle, yqyuan
		double offnadir = dotproduct(xyz_rho.crd_cvect(), _crs_sat_pco.crd_cvect()) / xyz_rho.norm() / _crs_sat_pco.norm();
		offnadir = acos(offnadir);
		_crt_obs.addnadir(offnadir);
		_crt_obs.addzen_sat(offnadir);

		double azi = atan2(neu_s[1], neu_s[0]);
		if (azi < 0)
		{
			azi += 2 * G_PI;
		}
		_crt_obs.addazi_rec(azi);
		_crt_obs.addzen_rec((G_PI / 2.0 - ele));

		/// for satellite-side azimuth, added by yqyuan for satellite-side aizmuth-dependent PCV
		//Matrix _rot_matrix_scf2crs = _RotMatrix_Ant(_crt_obs, _crt_sat_epo, sat_obj, true);
		Matrix _rot_matrix_scf2crs = _RotMatrix_Ant(_crt_obs, _crt_epo, _crt_sat_epo, sat_obj, true);
		t_gtriple xyz_s2r = t_gtriple(_rot_matrix_scf2crs.t() * ((-1) * xyz_rho.crd_cvect())); // from sat. to rec. in SCF XYZ
		double azi_sat = atan2(xyz_s2r[0], xyz_s2r[1]);
		if (azi_sat < 0)
			azi_sat += 2 * G_PI;
		_crt_obs.addazi_sat(azi_sat);

		///add for another elev and azi used in calculating weight matric
		t_gtriple xyz_rh = _trs_sat_crd - _trs_rec_crd;
		t_gtriple ell_(0, 0, 0), neu_sa(0, 0, 0), xRec(0, 0, 0), xyz_s(0, 0, 0);
		xyz2ell(_trs_rec_crd, ell_, false);
		xyz2neu(ell_, xyz_rh, neu_sa);
		double rho0 = sqrt(pow(_trs_rec_crd[0] - xyz_s[0], 2) + pow(_trs_rec_crd[1] - xyz_s[1], 2) + pow(_trs_rec_crd[2] - xyz_s[2], 2));
		double dPhi = OMEGA * rho0 / CLIGHT;
		xRec[0] = _trs_rec_crd[0] * cos(dPhi) - _trs_rec_crd[1] * sin(dPhi);
		xRec[1] = _trs_rec_crd[1] * cos(dPhi) + _trs_rec_crd[0] * sin(dPhi);
		xRec[2] = _trs_rec_crd[2];
		double NE2_ = neu_sa[0] * neu_sa[0] + neu_sa[1] * neu_sa[1];
		double ele_ = acos(sqrt(NE2_) / _crt_obs.rho());
		_crt_obs.addele_leo(ele_);

		// check elevation cut-off
		if (_crt_obj->id_type() == t_gdata::REC && _crt_obs.ele_deg() < _minElev)
		{
			if (_spdlog)
				SPDLOG_LOGGER_WARN(_spdlog, "Prepare fail! the elevation is too small");
			return false;
		}

		return true;
	}

	bool t_gprecisebiasFGO::_prt_obs_ALL(const t_gtime &crt_epo, t_gsatdata &obsdata, t_gallpar &pars, t_gobs &gobs, vector<pair<int, double>> &coeff)
	{
		auto par_list = pars.getPartialIndex(_crt_rec, _crt_sat);
		for (const int& ipar : par_list)
		{
			t_gpar par = pars.getPar(ipar);
			double coeff_value = 0.0;
			t_gbiasmodel::_Partial_basic(crt_epo, _crt_obs, gobs, par, coeff_value);
			if (coeff_value != 0.0)
			{
				//note: This is bad from 1 in coeff vector
				coeff.push_back(make_pair(ipar + 1, coeff_value));
			}
		}
		return true;
	}
}
