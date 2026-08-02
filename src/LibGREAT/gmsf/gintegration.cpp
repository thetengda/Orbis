/**
 * @file         gintegration.cpp
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        loosely/tightly coupled integration
 * @version      1.0
 * @date         2025-01-01
 *
 * @copyright Copyright (c) 2025, Wuhan University. All rights reserved.
 *
 */


#include "gintegration.h"
#include "gproc/gpreproc.h"
#include "gutils/gtimesync.h"
#include "gins/gins.h"
#include <iostream>
using namespace great;
using namespace std;

/******************************************* t_gintegration ***********************************************/
great::t_gintegration::t_gintegration(string site, string site_base, t_gsetbase* gset, std::shared_ptr<spdlog::logger> spdlog,t_gallproc* allproc) :
    t_gsinskf(gset,spdlog, (/*site_base+*/site)),
    t_gspp(site, gset, spdlog),
    t_gpvtflt(site, site_base, gset, spdlog, allproc)
{
    _spdlogkf = _spdlog;
    kftk = 0.0; measflag = 0;
    Ft = Pk = Eigen::MatrixXd::Zero(nq, nq);
    Hk = Eigen::MatrixXd::Zero(nr, nq);
    Qt = Xk = Eigen::VectorXd::Zero(nq); 
    Zk = Eigen::VectorXd::Zero(nr); 
    Rk = Eigen::MatrixXd::Zero(nr, nr);
    Phik = Eigen::MatrixXd::Identity(nq, nq);
    _global_variance = Eigen::MatrixXd::Identity(nq, nq);
    _gnss_beg = dynamic_cast<t_gsetgen*>(gset)->beg();
    _gnss_end = dynamic_cast<t_gsetgen*>(gset)->end();
    double start = dynamic_cast<t_gsetins*>(gset)->start();
    double end = dynamic_cast<t_gsetins*>(gset)->end();
    _ins_beg = t_gtime(_gnss_beg.gwk(), start);
    _ins_end = t_gtime(_gnss_end.gwk(), end);
    _ign_type = dynamic_cast<t_gsetign*>(gset)->ign_type();
    _imudata = nullptr;
    _initial_merge = false;
    _publisher.Initialize();
}

int great::t_gintegration::processBatchFB(const t_gtime& beg, const t_gtime& end, bool beg_end)
{
#ifdef BMUTEX   
    boost::mutex::scoped_lock lock(_mutex);
#endif 


    if (_grec == nullptr) {
        ostringstream os;
        os << "ERROR: No object found (" << _site << "). Processing terminated!!! " << beg.str_ymdhms() << " -> " << end.str_ymdhms() << endl;
        if (_spdlog) SPDLOG_LOGGER_ERROR(_spdlog, string("gintegration:  ") +  os.str());
        //_gmutex.unlock(); 
        return -1;
    }

    t_gpvtflt::InitProc(beg, end);

    if (!this->_init()) {
        if (_spdlog) SPDLOG_LOGGER_ERROR(_spdlog, string("gintegration:  ") +  "init failed");
        return -1;
    }
    _gobs->setepoches(_site);

	// time align
	if (_ins_beg < _gnss_beg) {
		_ins_crt = _imudata->erase_bef(_gnss_beg);
	}
	else
	{
		t_gpvtflt::processBatch(_gnss_crt, _ins_crt, false);
		if (_gnss_crt < _ins_crt) {
			int nEpo = round(_ins_crt.diff(_gnss_crt) / (_sampling) + 0.5);
			if (_sampling > 1) _gnss_crt.add_secs(int(_sampling * nEpo));  //  < 1Hz data
			else               _gnss_crt.add_dsec(_sampling * nEpo);       //  >=1Hz data
		}
	}

    std::cerr << endl << _site << ": Start MSF Processing:" << _ins_crt.str_ymdhms() << "  " << _end_time.str_ymdhms() << endl;
    cout << "EPOCH " << _ins_crt.str_ymdhms() << endl;
    t_gposdata::data_pos posdata;
    bool time_loop = true;

    // ins time loop
    while (time_loop)
    {



		


        if (_ins_crt > _ins_end || _gnss_crt > _gnss_end)
        {
            cerr << "\n" << "Processing Finished!" << endl;
            if (_spdlog) SPDLOG_LOGGER_INFO(_spdlog, string("gintegration:  ") +  "Processing Finished!");
            time_loop = false; break;
        }

        if (!_imudata->load(_wm, _vm, _shm.t, _shm.ts, _shm.nSamples))
        {
            cerr << "\n" << "Processing Finished!" << endl;
            break;
        }
        _ins_crt = t_gtime(_gnss_crt.gwk(), _shm.t);
        if (!_aligned)
        {
            if (_time_valid(_gnss_crt, _ins_crt))
            {
                Flag = _getPOS(posdata);
            }
            _aligned = cascaded_align(posdata.pos, posdata.vn);
        }
        else
        {
            if (_ign_type == IGN_TYPE::TCI && !_initial_merge) _merge_init();
            sins.Update(_wm, _vm, _shm);
            time_update(sins.nts);
            kftk = sins.t;
            int imeas = _getMeas(), irc = 0;

            string outinfo = "";
            set<MEAS_TYPE>::const_iterator it = _Meas_Type.begin();


			
            while (it != _Meas_Type.end())
            {
                switch (*it)
                {
                case GNSS_MEAS:
					irc = _GNSS_Update();
					break;
                default:break;
                }
                if (irc > 0) {
                    outinfo += meas2str(*it) + " ";
                    t_gintegration::_feedback(*it, true);
                    ++it;
                }
                else {
                    t_gintegration::_feedback(*it, false);
                    it = _Meas_Type.erase(it);
                }


				



            }
			


            Pk.conservativeResize(nq, nq); Xk.conservativeResize(nq);

            if (_Meas_Type.size()) {
                _imu_state.orientation = Eigen::Quaterniond(t_gbase::q2mat(sins.qnb));
                _imu_state.position = Geod2Cart(sins.pos, false);
                _publisher.UpdateNewState(_imu_state);
                sins.pure_ins_time = 0.0;
                
                double percent = _ins_crt.diff(_ins_beg) / _ins_end.diff(_ins_beg) * 100.0;
                cerr << "\r" << _ins_crt.str_ymdhms("Processing Epoch: ") << " Meas = " << meas2str(*_Meas_Type.begin()) << fixed << setprecision(1) << setw(6) << percent << "%";
                if (_spdlog) SPDLOG_LOGGER_INFO(_spdlog, string("gintegration:  ") + _site + ": ( " + outinfo + ")" +
                    _ins_crt.str_ymdhms(" integrated navigation processing epoch: "));
            }

            if ((_ign_type == IGN_TYPE::LCI || (_ign_type == IGN_TYPE::TCI && !_amb_state))
                && fabs(kftk - int(kftk)) < _shm.delay) {
                t_gintegration::_write();
            }






        }
        while (_ins_crt.diff(_gnss_crt) > _sampling)
        {
            if (_sampling > 1) _gnss_crt.add_secs(int(_sampling));  // =<1Hz data
            else               _gnss_crt.add_dsec(_sampling);       //  >1Hz data
        }



    }

    return 0;
}

int great::t_gintegration::_gnss_feedback(const ColumnVector& dx)
{
    try {
        // parameter dx
        int ign_size = _param.parNumber();
        int ins_size = nq;

        // trophere
        for (unsigned int iPar = 0; iPar < _param.parNumber(); iPar++) {
            if (_param[iPar].parType == par_type::TRP)
            {
                string site = _param[iPar].site;
                t_gtriple Ell, XYZ;
                if (_param.getCrdParam(site, XYZ) > 0) {
                }
                else { XYZ = _gallobj->obj(site)->crd_arp(_epoch); }
                xyz2ell(XYZ, Ell, false);
                if (site == _site && _gModel->tropoModel() != 0)
                    _param[iPar].apriori(_gModel->tropoModel()->getZHD(Ell, _epoch));
                else if (site == _site_base && _gModel_base->tropoModel() != 0)
                    _param[iPar].apriori(_gModel_base->tropoModel()->getZHD(Ell, _epoch));
            }
        }

        // parameter transmit
        for (unsigned int iPar = 0; iPar < _param.parNumber(); iPar++) {
            if (_param[iPar].parType == par_type::CRD_X || _param[iPar].parType == par_type::CRD_Y || _param[iPar].parType == par_type::CRD_Z)
            {
                _param[iPar].value(_param[iPar].value() - dx(_param[iPar].index));
            }
            else
            {
                _param[iPar].value(_param[iPar].value() + dx(_param[iPar].index));
            }
        }
        return 1;
    }
    catch (...)
    {
        if (_spdlog) SPDLOG_LOGGER_INFO(_spdlog, string("gintegration:  ") +  "GNSS feedback failed!");
        //_reset();
        Flag = NO_MEAS;
        return -1;
    }
}

// Separate the correction of each subsystem from Xk and sync var-matrix
void great::t_gintegration::_feedback(MEAS_TYPE type, bool successed)
{
	if (type == MEAS_TYPE::GNSS_MEAS && _ign_type == TCI)
	{
        t_gsins sins_out = sins; Eigen::VectorXd Xk_out = Xk;
        _global_variance = BaseMatrix2Eigen(_Qx);
		if (_amb_state) Xk = Columns2VectorXd(_filter->dx()).block(0, 0, nq, 1);
		this->feedback();
		if (_amb_state)this->_write();
		sins = sins_out; Xk = Xk_out;
	}
    else
    {
        if (successed) _global_variance = Pk;
        else  _global_variance = _gv_sav;
    }

    int global_size = _global_variance.size();
    if(!successed) Xk = Eigen::VectorXd::Zero(global_size);
    Pk = _global_variance; 
    if (_ign_type == TCI) _Qx = Eigen2BaseMatrix(_global_variance);

    // GNSS feedback
    if (_ign_type == IGN_TYPE::TCI)
    {
        ColumnVector dx; dx.ReSize(global_size);
        dx = VectorXd2Columns(Xk);
        _gnss_feedback(dx);
    }

    Xk.conservativeResize(nq);
    t_gsinskf::feedback();
    Xk = Eigen::VectorXd::Zero(global_size); 
    _gv_sav = _global_variance; 
}
void great::t_gintegration::_write(string info)
{
    set<string> ambs = _param.amb_prns();
    int nsat = ambs.size();
    // get amb status
    string amb = "Float";
    if (_amb_state)amb = "Fixed";
    string meas = "INS";
	
    if (_Meas_Type.size())    meas = meas2str(*_Meas_Type.begin()); //  ( GNSS -> ZUPT -> ODO -> NHC )
    _prt_ins_kml();
    ostringstream os;
    sins.prt_sins(os);
    os << fixed << setprecision(0)
        << " " << setw(10) << meas            // meas
        << " " << setw(5) << nsat            // nsat
        << fixed << setprecision(2)
        << " " << setw(7) << _dop.pdop()     // pdop
        << fixed << setprecision(2)
        << " " << setw(8) << amb
        << setw(10) << (_amb_state ? _ambfix->get_ratio() : 0.0);
    os << endl;
    _fins->write(os.str().c_str(), os.str().size());
    os.str("");
}

int great::t_gintegration::_prt_ins_kml()
{
    Eigen::Vector3d Geo_pos = Eigen::Vector3d(sins.pos(0) / t_gglv::deg, sins.pos(1) / t_gglv::deg, sins.pos(2));
    bool ins = dynamic_cast<t_gsetinp*>(_set)->input_size("imu") > 0 ? true : false;

    double crt = _ins_crt.sow() + _ins_crt.dsec();
    Eigen::Vector3d Qpos = Pk.block(6, 6, 3, 3).diagonal(), Qvel = Pk.block(3, 3, 3, 3).diagonal();
    if (_amb_state && _ign_type == IGN_TYPE::TCI)
    {
        Qpos = _global_variance.block(6, 6, 3, 3).diagonal();
        Qvel = _global_variance.block(3, 3, 3, 3).diagonal();
    }
    Eigen::Vector3d position = sins.pos_ecef, velocity = sins.ve;
    t_gposdata::data_pos posdata = t_gposdata::data_pos{ crt, position, velocity, Qpos, Qvel, 1.3, int(_data.size()), _amb_state };

    // write kml
    if (_kml && ins) {
        ostringstream out;
        out << fixed << setprecision(11) << " " << setw(0) << Geo_pos[1] << ',' << Geo_pos[0];
        string val = out.str();

        xml_node root = _doc;
        xml_node node = this->_default_node(root, _root.c_str());
        xml_node document = node.child("Document");
        xml_node last_child = document.last_child();
        xml_node placemark = document.insert_child_after("Placemark", last_child);
        string q = "#P" + _quality_grade(posdata);
        this->_default_node(placemark, "styleUrl", q.c_str());
        this->_default_node(placemark, "time", int2str(_ins_crt.sow()).c_str());
        xml_node point = this->_default_node(placemark, "Point");
        this->_default_node(point, "coordinates", val.c_str()); // for point
        xml_node description = placemark.append_child("description");
        description.append_child(pugi::node_cdata).set_value(_gen_kml_description(_ins_crt, posdata).c_str());
        xml_node TimeStamp = placemark.append_child("TimeStamp");
        string time = trim(_ins_crt.str_ymd()) + "T" + trim(_ins_crt.str_hms()) + "Z";
        this->_default_node(TimeStamp, "when", time.c_str());

        xml_node Placemark = document.child("Placemark");
        xml_node LineString = Placemark.child("LineString");
        this->_default_node(LineString, "coordinates", val.c_str(), false);  // for line
    }

    return 1;
}

int great::t_gintegration::_init()
{

    if (!_ins_init() || !_gnss_init())
    {
        if (_spdlog) SPDLOG_LOGGER_ERROR(_spdlog, string("gintegration:  ") +  "initialize integration filter failed! ");
        return -1;
    }

    _ins_crt = _gnss_crt = t_gtime::current_time(t_gtime::GPS);
    _gnss_beg = _beg_time;
    _gnss_end = _end_time;

    if (_gobs)
    {
        t_gtime start = _gobs->beg_obs(_site);
        t_gtime end = _gobs->end_obs(_site);
        t_gtime start_set = dynamic_cast<t_gsetgen*>(_setkf)->beg();
        t_gtime end_set = dynamic_cast<t_gsetgen*>(_setkf)->end();
        start = (start < start_set) ? start_set : start;
        end = (end < end_set) ? end : end_set;
        _gnss_beg = start;
        _gnss_end = end;
        _gnss_crt = _gnss_beg;
    }

    if (_imudata)
    {
        double start = _imudata->beg_obs();
        double end = _imudata->end_obs();
        double start_set = dynamic_cast<t_gsetins*>(_set)->start();
        double end_set = dynamic_cast<t_gsetins*>(_set)->end();
        start = (start < start_set) ? start_set : start;
        end = (end < end_set) ? end : end_set;

		_ins_beg = t_gtime(_gnss_beg.gwk(), start);
		_ins_end = t_gtime(_gnss_end.gwk(), end);
		_ins_crt = _ins_beg;
    }
    else
    {
        if (_spdlog) SPDLOG_LOGGER_ERROR(_spdlog, string("gintegration:  ") +  "IMU observation is not existing!!! ");
        return -1;
    }
    
    _ins_end = _gnss_end < _ins_end ? _gnss_end : _ins_end;
    _gnss_end = _ins_end;

    return 1;
}


bool great::t_gintegration::_gnss_init()
{
    if ((_ign_type == IGN_TYPE::TCI || _ign_type == IGN_TYPE::VIS_TCI || _ign_type == IGN_TYPE::LIDAR_TCI || _ign_type == IGN_TYPE::MULTIGN_TCI || _ign_type == IGN_TYPE::UWB_TCI))
    {
        try
        {
            SymmetricMatrix Qx_extended;

            int last_ins_idx = t_gsinskf::get_last_idx();

            int icrdz_gnss = _param.getParam(_site, par_type::CRD_Z, "");
            for (int i = icrdz_gnss + 1; i < _param.parNumber(); i++)
            {
                param_of_sins.addParam(_param[i]);
            }

            Qx_extended.ReSize(param_of_sins.parNumber()); Qx_extended = 0.0;

            if (last_ins_idx + 2 < param_of_sins.parNumber())
            {
                Qx_extended.SymSubMatrix(last_ins_idx + 2, param_of_sins.parNumber()) =
                    _Qx.SymSubMatrix(icrdz_gnss + 2, _param.parNumber());
            }

            for (int i = 0; i < last_ins_idx + 1; i++)
            {
                Qx_extended(i + 1, i + 1) = 1.0;
            }

            _param = param_of_sins;
            _Qx = Qx_extended;

            _param.reIndex();
            param_of_sins.reIndex(); 

        }
        catch (...)
        {
            return false;
        }
    }
    return true;
}



int great::t_gintegration::_processEpoch(const t_gtime& runEpoch)
{

    if (_grec == nullptr) {
        if (_spdlog) SPDLOG_LOGGER_INFO(_spdlog, string("gintegration:  ") +  "No receiver settings available!!!");
        return NO_MEAS;
    }

    _epoch = runEpoch;
    _amb_state = false;

    if (!_crd_xml_valid())
        _sig_init_crd = 100.0;

    Matrix A;
    SymmetricMatrix P;
    ColumnVector l, dx;
    ColumnVector v_orig, v_norm;
    SymmetricMatrix Qsav, QsavBP;
    t_gallpar XsavBP;
    double vtpv;
    int nobs_total, npar_number;
    string outlier = "";

    _cntrep = 0;    // number of iterations caused by outliers


    do
    {
        // erase sat satellite because of outlier
        _remove_sat(outlier);

		if (_valid_ins_constraint()) {
			int crd_idx = _param.getParam(_site, par_type::CRD_X, "");
			t_gtriple xyz_gnss(sins.pos_ecef + sins.Ceb * lever);
			t_gtriple Qxyz(sqrt(_Qx(crd_idx + 1, crd_idx + 1)), sqrt(_Qx(crd_idx + 2, crd_idx + 2)), sqrt(_Qx(crd_idx + 3, crd_idx + 3)));
			_external_pos(xyz_gnss, Qxyz);
		}

        if (_prepareData() < 0) {
            _predict(runEpoch);
            return NO_MEAS;
        }
        QsavBP = _Qx;
        XsavBP = _param;

        _predict(runEpoch);



        _initialized = true;
        _pos_constrain = false;

        if (_data.size() < _minsat) {
            if (_spdlog) SPDLOG_LOGGER_INFO(_spdlog, string("gintegration:  ") +  "Not enough visible satellites!");
            t_gsppflt::_restore(QsavBP, XsavBP);
            return NO_MEAS;
        }

        // define a number of measurements 
        unsigned int nObs = _data.size();
        unsigned int mult = 1;
        if (_observ == OBSCOMBIN::RAW_ALL )
        {
            mult = 2;
            nObs *= 5;
        } 
        if (_observ == OBSCOMBIN::RAW_MIX)
        {
            mult = 1;
            nObs *= 5;
        }
        if (_phase)
        {
            mult *= 2;
            nObs *= 2;
        } 

        unsigned int nPar = _param.parNumber();
        unsigned int iobs = 1;
        _frqNum.clear(); _obs_index.clear();

        if (_isBase)
        {
            dynamic_cast<t_gcombDD*>(&(*_base_model))->set_base_data(&_data_base);
            dynamic_cast<t_gcombDD*>(&(*_base_model))->set_rec_info(_gallobj->obj(_site_base)->crd(_epoch), _vBanc(4), _vBanc_base(4));
        }




        // use combmodel
        t_gfltEquationMatrix equ;
        iobs = _cmp_equ(equ);
        equ.chageNewMat(A, P, l, nPar); 

        dx.ReSize(nPar);
        dx = 0.0;
        // generate obs_index
        _obs_index.clear();
        _generateObsIndex(equ);


        if (iobs < _minsat * mult) {
            if (_spdlog) SPDLOG_LOGGER_INFO(_spdlog, string("gintegration:  ") + "Not enough processable observations!");
            t_gsppflt::_restore(QsavBP, XsavBP);
            return NO_MEAS;
        }

        if (_isBase) {
            if (_combineDD(A, P, l) < 0) {
                if (_spdlog) SPDLOG_LOGGER_INFO(_spdlog, string("gintegration:  ") + _site + runEpoch.str_ymdhms(" combining DD observation failed at epoch: "));
                return NO_MEAS;
            }
        }
        Qsav = _Qx;

        if (_merge_pose(A) < 0) {
            if (_spdlog) SPDLOG_LOGGER_INFO(_spdlog, string("gintegration:  ") + _site + runEpoch.str_ymdhms(" merge failed at epoch: "));
            dx = 0.0; _Qx = Qsav;
            return NO_MEAS;
        }
        try
        {
            _filter->update(A, P, l, dx, _Qx);
        }
        catch (...)
        {
            if (_spdlog) SPDLOG_LOGGER_INFO(_spdlog, string("gintegration:  ") +  _site + runEpoch.str_ymdhms(" meas_update failed at epoch: "));
            dx = 0.0; _Qx = Qsav;
            return NO_MEAS;
        }

		

        // increasing variance after update in case of introducing new ambiguity
        //for (size_t iPar = 0; iPar < _param.parNumber(); iPar++) {
        //    if (_param[iPar].parType == par_type::AMB_IF ||
        //        _param[iPar].parType == par_type::AMB_L1 ||
        //        _param[iPar].parType == par_type::AMB_L2 ||
        //        _param[iPar].parType == par_type::AMB_L3 ||
        //        _param[iPar].parType == par_type::AMB_L4 ||
        //        _param[iPar].parType == par_type::AMB_L5) {
        //        string sat = _param[iPar].prn;
        //        if (_newAMB.find(sat) != _newAMB.end() && _cntrep == 1) {
        //            if (_newAMB[sat] == 1) _Qx(iPar + 1, iPar + 1) += 10;
        //            if (_newAMB[sat] == 2 && _Qx(iPar + 1, iPar + 1) > 0.01) _Qx(iPar + 1, iPar + 1) += 1;
        //            _newAMB[sat]++;
        //        }
        //    }
        //}

        // post-fit residuals
        v_orig = l - A * dx;
        ColumnVector vtPv = v_orig.t() * P * v_orig;
        vtpv = vtPv(1);
        nobs_total = A.Nrows();
        npar_number = A.ncols();
        int freedom = A.Nrows() - A.Ncols();
        if (freedom < 0)freedom = 1;
        _sig_unit = vtPv(1) / freedom;

        // normalized post-fit residuals
        Matrix Qv = A * _Qx * A.t() + P.i();
        ////////
        v_norm.ReSize(v_orig.Nrows());
        for (int i = 1; i <= v_norm.Nrows(); i++) {
            v_norm(i) = sqrt(1 / Qv(i, i)) * v_orig(i);
        }
    } while (_outlierDetect(v_norm, Qsav, outlier) != 0);


    _filter->add_data(_param, dx, _Qx, _sig_unit, Qsav);
    _filter->add_data(A, P, l);
    _filter->add_data(vtpv, nobs_total, npar_number);

    Xk = Columns2VectorXd(dx);


    _amb_resolution();
	
    return GNSS_MEAS;
}

int great::t_gintegration::_merge_pose(Matrix& A)
{
    try
    {
        int nobs = A.Nrows();
        // for position coff
        int icrdx = _param[_param.getParam(_site, par_type::CRD_X, "")].index;
        A.SubMatrix(1, nobs, icrdx, icrdx + 2) = -A.SubMatrix(1, nobs, icrdx, icrdx + 2); 
        // for attitude coff
        int iattx = _param[_param.getParam(_site, par_type::ATT_X, "")].index;
        A.SubMatrix(1, nobs, iattx, iattx + 2) = A.SubMatrix(1, nobs, icrdx, icrdx + 2) * Eigen2newMat(t_gbase::askew(sins.Ceb * lever));
        tmeas = _gnss_crt.sow() + _gnss_crt.dsec();
        if (fabs(sins.t - tmeas) > 0.5) throw("merge failed");
    }
    catch (...)
    {
        if (_spdlog) SPDLOG_LOGGER_INFO(_spdlog, string("gintegration:  ") +  "GNSS and INS infomation fusion failed!");
        //_reset();
        Flag = NO_MEAS;
        return -1;
    }
    return 1;
}

int great::t_gintegration::_merge_init()
{
    t_gtriple gnss_pos;
    _param.getCrdParam(_site, gnss_pos);   
    if (_shm.align)
    {
        sins.eth.Update(Cart2Geod(Eigen::Vector3d(gnss_pos[0], gnss_pos[1], gnss_pos[2]), false), Eigen::Vector3d::Zero());
        sins.Cnb = t_gbase::q2mat(sins.qnb);
        sins.pos_ecef = Eigen::Vector3d(gnss_pos[0], gnss_pos[1], gnss_pos[2]) - sins.eth.Cen * sins.Cnb * lever;
        sins.pos = Cart2Geod(sins.pos_ecef, false);
    }
    int icrdx = _param[_param.getParam(_site, par_type::CRD_X, "")].index;
    Pk.block(icrdx - 1, icrdx - 1, 3, 3) = NewMat2Eigen(_Qx).block(icrdx - 1, icrdx - 1, 3, 3);
    for (int i = icrdx + 3; i <= _param.parNumber(); i++)
        for (int j = icrdx; j <= icrdx + 2; j++)
            _Qx(i, j) = -_Qx(i, j);
    _global_variance = BaseMatrix2Eigen(_Qx);
    _initial_merge = true;
    return 1;
}

int great::t_gintegration::_GNSS_Update()
{

    t_gposdata::data_pos posdata;

    if (_ign_type == IGN_TYPE::LCI)
    {
        _global_variance = Pk;
        Flag = _getPOS(posdata);
        if (Flag != NO_MEAS) {
            _meas_update();
        }
    }
    else if (_ign_type == IGN_TYPE::TCI)
    {
        _timeUpdate(_gnss_crt);
        Flag = MEAS_TYPE(t_gintegration::_processEpoch(_gnss_crt));
        if (Flag == NO_MEAS)
        {
            if (_spdlog) SPDLOG_LOGGER_INFO(_spdlog, string("gintegration:  ") + _gnss_crt.str_ymdhms("TC GNSS/INS filter updating Failed at epoch: "));
        }
        else
        {
            if (_spdlog) SPDLOG_LOGGER_INFO(_spdlog, string("gintegration:  ") + _gnss_crt.str_ymdhms("TC GNSS/INS filter updating Succeed at epoch: "));
        }
    }

    if (Flag == NO_MEAS)return -1;

    return 1;
}

MEAS_TYPE great::t_gintegration::_getPOS(t_gposdata::data_pos& pos)
{
    MEAS_TYPE res_type;
    double crt = _ins_crt.sow() + _ins_crt.dsec();
    t_gtime runEpoch = _gobs->load(_site, crt);
    int irc = ProcessOneEpoch(runEpoch);
    if (irc < 0) {
        return MEAS_TYPE::NO_MEAS;
    }
    _get_result(runEpoch, pos);

    MeasVel = pos.vn; MeasPos = pos.pos; tmeas = pos.t;
    _Cov_MeasVn = pos.Rvn; _Cov_MeasPos = pos.Rpos;

    if (!_aligned)
    {
        sins.pos = Cart2Geod(MeasPos, false);
        sins.vn = t_gbase::Cen(sins.pos).transpose() * MeasVel;
    }
    res_type = MEAS_TYPE::POS_MEAS;

    if (pos.PDOP > _shm.max_pdop)res_type = NO_MEAS;
    if (pos.nSat < _shm.min_sat)res_type = NO_MEAS;
    if (_isBase && !pos.amb_state)res_type = NO_MEAS;

    return res_type;

}

int great::t_gintegration::_getMeas()
{
    Flag = NO_MEAS; 
    _Meas_Type.clear(); _amb_state = false;

    if (t_gintegration::_time_valid(_gnss_crt, _ins_crt))
    {
        _Meas_Type.insert(MEAS_TYPE::GNSS_MEAS);
    }

    if (_Meas_Type.size())
    {
        int size = _global_variance.rows();
        _global_variance.block(0, 0, nq, nq) = Pk;
        _global_variance.block(nq, 0, size - nq, nq) = _global_variance.block(nq, 0, size - nq, nq) * Phik.transpose();
        _global_variance.block(0, nq, nq, size - nq) = Phik * _global_variance.block(0, nq, nq, size - nq);
        if (_ign_type == TCI) _Qx = Eigen2BaseMatrix(_global_variance);
        Pk = _global_variance; Xk.conservativeResize(Pk.rows()); Xk.block(nq, 0, size - nq, 1) = Eigen::VectorXd::Zero(size - nq);
        Phik = Eigen::MatrixXd::Identity(nq, nq);
        _gv_sav = _global_variance; // restore
    }
    return _Meas_Type.size();
}

bool great::t_gintegration::_time_valid(t_gtime gt, t_gtime inst)
{
    bool res_valid = false;
    double crt = inst.sow() + inst.dsec();
    _data.erase(_data.begin(), _data.end());
    t_gtime obsEpo = _gobs->load(_site, crt);

    if ((abs(inst.diff(obsEpo)) < _shm.delay && inst >= obsEpo)
        || abs(inst.diff(obsEpo)) < 1e-6)
    {
		this->_slip_detect(obsEpo);
		
        _data = _gobs->obs(_site, obsEpo);
		if (_data.size() > 0) {
			res_valid = true;
			if (_gallbias) {
				for (auto& itdata : _data) {
					itdata.apply_bias(_gallbias);
				}
			}
		}
        else {
            if (_spdlog) {
                SPDLOG_LOGGER_INFO(_spdlog, string("gintegration:  ") +  _site + gt.str_ymdhms(" no observation found at epoch: "));
            }
            res_valid = false;
        }
        if (_isBase)
        {
            if (true) {
                _data_base.erase(_data_base.begin(), _data_base.end());
                _data_base = _gobs->obs(_site_base, obsEpo);
            }
            if (_data_base.size() > 0) {
                res_valid = true;
                // apply dcb
                if (_gallbias) {
                    for (auto& itdata : _data_base) {
                        itdata.apply_bias(_gallbias);
                    }
                }
            }
            else {
                if (_spdlog) {
                    SPDLOG_LOGGER_INFO(_spdlog, string("gintegration:  ") +  _site_base + gt.str_ymdhms(" no base observation found at epoch: "));
                }
                res_valid = false;
            }
        }
    }
    // judge ultimately(for rtk)
    if (_data.size() == 0)res_valid = false;
    if (res_valid)_gnss_crt = obsEpo;
    return res_valid;
}
