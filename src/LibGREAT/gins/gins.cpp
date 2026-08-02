/**
 * @file         gins.cpp
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        INS core algorithm
 * @version      1.0
 * @date         2025-01-01
 *
 * @copyright Copyright (c) 2025, Wuhan University. All rights reserved.
 *
 */

#include "gins/gins.h"
#include <iostream>
#include <iomanip>
#include "gset/gsetout.h"
#include "gset/gsetins.h"
#include "gset/gsetign.h"
#include <Eigen/src/Core/util/DisableStupidWarnings.h>
#include <gmsf/gintegration.h>
#include <map>
#include <mutex>
using namespace Eigen;
using namespace std;
using namespace great;

namespace
{
    struct ProcessCovarianceDiagnosticState
    {
        Eigen::MatrixXd covariance;
        double dt;
        int steps;

        ProcessCovarianceDiagnosticState() : dt(0.0), steps(0) {}
    };

    std::map<const great::t_gsinskf*, ProcessCovarianceDiagnosticState>
        process_covariance_diagnostics;
    std::mutex process_covariance_diagnostic_mutex;

    ProcessCovarianceDiagnosticState& process_covariance_state_unlocked(
        const great::t_gsinskf* filter)
    {
        return process_covariance_diagnostics[filter];
    }
}


/******************************************* t_gsins ***********************************************/
t_gsins::t_gsins(const t_gquat& qnb0, const Vector3d& vn0, const Vector3d& pos0, double t0)
{
    qnb = qnb0; vn = vn0; pos = pos0; t = t0;
    eth.Update(pos0, vn0);
    Cnb = t_gbase::q2mat(qnb); att = t_gbase::q2att(qnb); Cbn = Cnb.transpose(); vb = Cbn*vn;
    Kg = Ka = Vector3d::Ones();
    eb = db = _tauG = _tauA = _tauGScale = _tauAScale = Vector3d::Zero(); _tauOdo = 0;
    wib = web = wnb = fb = fn = an = Vector3d::Zero();
    xyz_out[0] = 0.0;
}

great::t_gsins::t_gsins(gnut::t_gsetbase * set)
{
    Vector3d Cart= dynamic_cast<t_gsetins*>(set)->pos();
    pos = Cart2Geod(Cart, false);
    vn = t_gbase::Cen(pos).transpose()*(dynamic_cast<t_gsetins*>(set)->vel());
    att = dynamic_cast<t_gsetins*>(set)->att();
    eth.Update(pos, vn);
    qnb = t_gbase::a2qua(att);
    Cnb = t_gbase::q2mat(qnb); Cbn = Cnb.transpose(); vb = Cbn*vn;
    Kg = Ka = Vector3d::Ones();
    eb = db = _tauG = _tauA = _tauGScale = _tauAScale = Vector3d::Zero(); _tauOdo = 0;
    eb = dynamic_cast<t_gsetins*>(set)->gyro_bias() * t_gglv::dph;
    db = dynamic_cast<t_gsetins*>(set)->acce_bias()* t_gglv::mg;
    wib = web = wnb = fb = fn = an = Vector3d::Zero();
    Eigen::Vector3d installation_att = dynamic_cast<t_gsetins*>(set)->imu_installation_rotation();
    qvb = t_gbase::a2qua(installation_att); Cvb = t_gbase::q2mat(qvb); Cbv = Cvb.transpose();
    pure_ins_time = 1.0; 
}

void great::t_gsins::set_posvel(const Eigen::Vector3d & pos0, const Eigen::Vector3d & vn0)
{
    pos = pos0; vn = vn0;
    eth.Update(pos, vn);
}

Vector3d t_gsins::align_coarse(const Vector3d& wmm,const Vector3d& vmm)
{
    double latitude = pos(0);
    double T11, T12, T13, T21, T22, T23, T31, T32, T33;
    double cb = cos(latitude), tb = tan(latitude), nn;
    Vector3d wbib = wmm / wmm.norm(), fb = vmm / vmm.norm();
    T31 = fb(0), T32 = fb(1), T33 = fb(2);
    T21 = wbib(0) / cb - T31*tb, T22 = wbib(1) / cb - T32*tb, T23 = wbib(2) / cb - T33*tb;       
    nn = sqrt(T21*T21 + T22*T22 + T23*T23);  T21 /= nn, T22 /= nn, T23 /= nn;
    T11 = T22*T33 - T23*T32, T12 = T23*T31 - T21*T33, T13 = T21*T32 - T22*T31;
    nn = sqrt(T11*T11 + T12*T12 + T13*T13);  T11 /= nn, T12 /= nn, T13 /= nn;
    Matrix3d Cnb;
    Cnb << T11, T12, T13, T21, T22, T23, T31, T32, T33;
    return t_gbase::m2att(Cnb);
}

void t_gsins::Update(const vector<Vector3d>& wm, const vector<Vector3d>& vm, const t_scheme& scm)
{
    nts = scm.ts;
    t = scm.t;
    nts = abs(nts);
    double nts_2 = nts / 2.0;
    imu.Update(wm, vm, scm);
    imu.phim = Kg.asDiagonal()*imu.phim - eb*nts; imu.dvbm = Ka.asDiagonal() *imu.dvbm - db*nts;
    Vector3d vn1_2 = vn + an*nts_2, pos1_2 = pos + eth.v2dp(vn1_2, nts_2);
    eth.Update(pos1_2, vn1_2);
    wib = imu.phim / nts; fb = imu.dvbm / nts;
    web = wib - Cbn*eth.wnie;
    wnb = wib - t_gquat::conj(qnb * t_gbase::rv2q(imu.phim / 2)) * eth.wnin;
    fn = qnb*fb;
    an = t_gbase::rv2q(-eth.wnin*nts_2)*fn + eth.gcc;
    Vector3d vn1 = vn + an*nts;
    pos = pos + eth.v2dp(vn + vn1, nts_2);    vn = vn1;
    qnb = t_gbase::rv2q(-eth.wnin*nts)*qnb*t_gbase::rv2q(imu.phim);
    Cnb = t_gbase::q2mat(qnb); att = t_gbase::m2att(Cnb); Cbn = Cnb.transpose(); vb = Cbn*vn;
    eth.Update(pos, vn); pos_ecef = Geod2Cart(pos, false);
    Ceb = eth.Cen*Cnb; Cbe = Ceb.transpose(); qeb = t_gbase::m2qua(Ceb); ve = eth.Cen*vn; ae = eth.Cen*an;
    pure_ins_time += nts;
}

void great::t_gsins::prt_header(ostringstream& os, bool imu_scale, bool use_odo)
{
	// the first line
	os << "# ";
	os << setw(15) << "Seconds of Week";

	os << setw(18) << "X-ECEF" <<
		setw(18) << "Y-ECEF" <<
		setw(18) << "Z-ECEF";
	os << setw(10) << "VX" <<
		setw(10) << "VY" <<
		setw(10) << "VZ";
	os << setw(10) << "Pitch" <<
		setw(10) << "Roll" <<
		setw(10) << "Yaw";
	os << setw(12) << "GyroBiasX" <<
		setw(12) << "GyroBiasY" <<
		setw(12) << "GyroBiasZ" <<
		setw(12) << "AcceBiasX" <<
		setw(12) << "AcceBiasY" <<
		setw(12) << "AcceBiasZ";
	if (imu_scale)
		os << setw(12) << "GyroScaleX" <<
		setw(12) << "GyroScaleY" <<
		setw(12) << "GyroScaleZ" <<
		setw(12) << "AcceScaleX" <<
		setw(12) << "AcceScaleY" <<
		setw(12) << "AcceScaleZ";
	if (use_odo) os << setw(12) << "OdoScale";
    os << setw(12) << "MeasType" <<
        setw(7) << "Nsat" <<
        setw(7) << "PDOP" <<
        setw(12) << "AmbStatus" <<
        setw(7) << "Ratio";
	os << "         " << endl;


	// the second line
	os << "# ";
	os << setw(15) << "(s)";
	os << setw(18) << "(m)" <<
		setw(18) << "(m)" <<
		setw(18) << "(m)";
	os << setw(10) << "(m/s)" <<
		setw(10) << "(m/s)" <<
		setw(10) << "(m/s)";
	os << setw(10) << "(deg)" <<
		setw(10) << "(deg)" <<
		setw(10) << "(deg)";
	os << setw(12) << "(deg/h)";
	os << setw(12) << "(deg/h)";
	os << setw(12) << "(deg/h)";
	os << setw(12) << "(mg)";
	os << setw(12) << "(mg)";
	os << setw(12) << "(mg)";

	if (imu_scale)
		os << setw(12) << "#" <<
		setw(12) << "#" <<
		setw(12) << "#" <<
		setw(12) << "#" <<
		setw(12) << "#" <<
		setw(12) << "#";
	if (use_odo)
		os << setw(12) << "#";
	os << setw(12) << " " <<
		setw(7) << "#" <<
		setw(7) << "#" <<
		setw(12) << " ";
	os << endl;
}

void great::t_gsins::prt_sins(ostringstream& os)
{
    // format
    Vector3d pos_out;
    std::string content = os.str();
    if (abs(xyz_out[0]) < 0.01||content[0] == '0')
    {
        pos_out = pos_ecef;
    }
    else pos_out = xyz_out;
    string os_content = os.str();

    if (!os_content.empty() && (os_content[0] == '0' || os_content[0] == '1')) {
        os.str(os_content.substr(1));
    }

    Vector3d vel_out = ve;
    Vector3d att_out = att / t_gglv::deg;
    Vector3d eb_out = eb / t_gglv::deg * t_gglv::hur;
    Vector3d db_out = db / t_gglv::mg;

    // for other format

    // output
	os << fixed << setprecision(6) << setw(18) << t;
	os << fixed << setprecision(3) <<
		setw(18) << pos_out(0) <<
		setw(18) << pos_out(1) <<
		setw(18) << pos_out(2);
	os << fixed << setprecision(3) <<
		setw(10) << vel_out(0) <<
		setw(10) << vel_out(1) <<
		setw(10) << vel_out(2);
	os << fixed << setprecision(4) <<
		setw(10) << att_out(0) <<
		setw(10) << att_out(1) <<
		setw(10) << att_out(2);
	os << fixed << setprecision(4) <<
		setw(12) << eb_out(0) <<
		setw(12) << eb_out(1) <<
		setw(12) << eb_out(2);
	os << fixed << setprecision(4) <<
		setw(12) << db_out(0) <<
		setw(12) << db_out(1) <<
		setw(12) << db_out(2);
}



void great::t_gsins::debug_ins_info()
{
    cout << "ins info : \n"
        << "     time: " << setw(16) << t << "\n"
        << "     att: " << setw(10) << att.transpose() / t_gglv::deg << "\n"
        << "     vel: " << setw(10) << ve.transpose() << "\n"
        << "     pos: " << setw(10) << pos_ecef.transpose() << "\n"
        << "     an: " << setw(10) << an.transpose() << "\n"
        << "     ae: " << setw(10) << ae.transpose() << "\n"
        << "     web: " << setw(10) << web.transpose() / t_gglv::dps << "\n"
        << "     wnb: " << setw(10) << wnb.transpose() / t_gglv::dps << "\n"
        << endl;
}


/******************************************* t_gsinskf ***********************************************/


great::t_gsinskf::t_gsinskf(gnut::t_gsetbase* gset, t_spdlog spdlog, string site)
{
    kftk = 0.0; measflag = 0;
    nq = dynamic_cast<t_gsetign*>(gset)->nq();
    nr = dynamic_cast<t_gsetign*>(gset)->nr();
    Ft = Pk = MatrixXd::Zero(nq, nq);
    Hk = MatrixXd::Zero(nr, nq);
    Qt = Xk = VectorXd::Zero(nq);
    Zk = VectorXd::Zero(nr);
    Rk = MatrixXd::Zero(nr, nr);
    reset_diagnostic_process_covariance();

    this->sins = t_gsins(gset);
    this->gset(gset);
    _shm = t_scheme(gset);
    _spdlogkf = spdlog;
    if (!site.empty()) _name = site;
    this->set_out(); this->init_par(_name);
    wmm = vmm = Eigen::Vector3d::Zero();
    if (_shm.align)_aligned = false;
    else _aligned = true;
    _align_count = 0;
    _first_align = true;
    if (_shm._imu_scale)nq += 6;  // gyro and acce scale factors
    if (_shm._imu_inst_rot)nq += 3;  // gyro and acce scale factors
    //if (_shm._odo) nq += 1;// default(in xml): 15, if(_odo)nq++;
    _Cov_MeasOdo = pow(dynamic_cast<t_gsetign*>(_setkf)->odo_std(), 2);
    _Cov_MeasNHC = dynamic_cast<t_gsetign*>(_setkf)->NHC_std().array().abs2();
    _Cov_MeasZUPT = dynamic_cast<t_gsetign*>(_setkf)->ZUPT_std().array().abs2();
    _Cov_MeasZIHR = pow(dynamic_cast<t_gsetign*>(_setkf)->ZIHR_std(), 2);
    _Cov_MeasAtt = dynamic_cast<t_gsetign*>(_setkf)->ATT_std().array().abs2();
    _map_maxnorm = dynamic_cast<t_gsetign*>(_setkf)->max_norm();
    _resample_intv = 1.0 / dynamic_cast<t_gsetins*>(_setkf)->resampled_freq();
}

great::t_gsinskf::~t_gsinskf()
{
    if (_fins) {delete _fins; _fins = nullptr;}
    std::lock_guard<std::mutex> lock(
        process_covariance_diagnostic_mutex);
    process_covariance_diagnostics.erase(this);
}

Eigen::MatrixXd
great::t_gsinskf::diagnostic_process_covariance() const
{
    std::lock_guard<std::mutex> lock(
        process_covariance_diagnostic_mutex);
    return process_covariance_state_unlocked(this).covariance;
}

double great::t_gsinskf::diagnostic_process_dt() const
{
    std::lock_guard<std::mutex> lock(
        process_covariance_diagnostic_mutex);
    return process_covariance_state_unlocked(this).dt;
}

int great::t_gsinskf::diagnostic_process_steps() const
{
    std::lock_guard<std::mutex> lock(
        process_covariance_diagnostic_mutex);
    return process_covariance_state_unlocked(this).steps;
}

void great::t_gsinskf::reset_diagnostic_process_covariance()
{
    std::lock_guard<std::mutex> lock(
        process_covariance_diagnostic_mutex);
    ProcessCovarianceDiagnosticState& state =
        process_covariance_state_unlocked(this);
    state.covariance = Eigen::MatrixXd::Zero(nq, nq);
    state.dt = 0.0;
    state.steps = 0;
}

void great::t_gsinskf::Add_IMU(t_gimudata* imu)
{
    _imudata = imu;
    return;
}

void great::t_gsinskf::init_par(string site)
{
    _name = site;
    // attitude
    for (int ipar = (int)par_type::ATT_X; ipar <= (int)par_type::ATT_Z; ipar++)
        param_of_sins.addParam(t_gpar(site, par_type(ipar), ipar, ""));
    // velocity
    for (int ipar = (int)par_type::VEL_X; ipar <= (int)par_type::VEL_Z; ipar++)
        param_of_sins.addParam(t_gpar(site, par_type(ipar), ipar, ""));
    // position 
    for (int ipar = (int)par_type::CRD_X; ipar <= (int)par_type::CRD_Z; ipar++)
        param_of_sins.addParam(t_gpar(site, par_type(ipar), ipar, ""));
    // bias
    for (int ipar = (int)par_type::eb_X; ipar <= (int)par_type::db_Z; ipar++)
        param_of_sins.addParam(t_gpar(site, par_type(ipar), ipar, ""));

    param_of_sins.reIndex();
}

// The last parameter should be consistent with the initialization order of the parameter list !!!
// see also great::t_gsinskf::init_par(string site)
int great::t_gsinskf::get_last_idx()
{
    int idx = param_of_sins.getParam(_name, par_type::db_Z, "");

    return idx;
}

bool great::t_gsinskf::align_coarse(const vector<Eigen::Vector3d>& wm, const vector<Eigen::Vector3d>& vm, const t_scheme & scm)
{
    sins.t = scm.t;
    for (int j = 0; j < scm.nSamples; j++) {
        wmm = wmm + wm[j];
        vmm = vmm + vm[j];
        _align_count++;
    }
    if (abs(_align_count - scm.align_time * scm.freq) < 1e-5)
    {
        cerr << "static alignment" << endl;
        wmm = wmm / _align_count; vmm = vmm / _align_count;
        sins.qnb = t_gbase::a2qua(sins.align_coarse(wmm, vmm));
        return true;
    }
    return false;
}

bool great::t_gsinskf::align_pva(const Eigen::Vector3d& pos, const t_scheme & scm)
{
    if (_first_align)
    {
        _first_pos = _pre_pos = pos;
        _first_align = false;
        return false;
    }
    Vector3d endpos = pos;
    Vector3d vel = endpos - _pre_pos; _pre_pos = endpos;
    Eigen::Vector3d blh = Cart2Geod(pos, false);
    Eigen::Vector3d vn = t_gbase::Cen(blh).transpose()*vel;
    sins.set_posvel(blh, vn);

    Vector3d baseline = XYZ2ENU(endpos, _first_pos);
    double dist = SQRT(SQR(baseline(0)) + SQR(baseline(1)));

    // double yaw;
    double dyaw = 10*t_gglv::deg;
    double pos_dist = dynamic_cast<t_gsetins*>(_setkf)->pos_dist();
    if (dist > pos_dist)
    {
        double yaw = acos(fabs(baseline(1)) / dist);
        if (baseline(1) > 0 && baseline(0) > 0)yaw = -yaw;
        if (baseline(1) > 0 && baseline(0) < 0)yaw = yaw;
        if (baseline(1) < 0 && baseline(0) > 0)yaw = -(t_gglv::PI - yaw);
        if (baseline(1) < 0 && baseline(0) < 0)yaw = t_gglv::PI - yaw;
        sins.qnb = t_gbase::a2qua(Vector3d(0, 0, yaw));
        if (_yaw0 == 0.0)
        {
            _yaw0 = yaw;
            _first_pos = endpos;
            return false;
        }
        if(fabs(yaw - _yaw0) < dyaw) return true;
        _yaw0 = 0.0;
        _first_align = true;
        cerr << "\r" << _ins_crt.str_ymdhms("Position vector alignment failed. Please keep straight: ") << endl;
    }

    return false;
}

bool great::t_gsinskf::align_vva(const Eigen::Vector3d& vel, const t_scheme & scm)
{
    double vel_norm = dynamic_cast<t_gsetins*>(_setkf)->vel_norm();
    if (vel.norm() > vel_norm)
    {
        double yaw = atan2(fabs(vel(0)), fabs(vel(1)));
        if (vel(1) > 0 && vel(0) > 0)yaw = -yaw;
        if (vel(1) > 0 && vel(0) < 0)yaw = yaw;
        if (vel(1) < 0 && vel(0) > 0)yaw = -(t_gglv::PI - yaw);
        if (vel(1) < 0 && vel(0) < 0)yaw = t_gglv::PI - yaw;
        sins.qnb = t_gbase::a2qua(Vector3d(0, 0, yaw));
        return true;
    }
    cerr << "\r" << _ins_crt.str_ymdhms("Velocity vector alignment failed. Please speed up: ") << endl;
    return false;
}

bool great::t_gsinskf::cascaded_align(const Eigen::Vector3d& pos, const Eigen::Vector3d& vel)
{
    sins.imu.Update(_wm, _vm, _shm);

    Eigen::Vector3d blh = Cart2Geod(pos, false);
    Eigen::Vector3d vn = t_gbase::Cen(blh).transpose()*vel;
    sins.set_posvel(blh, vn);

    ALIGN_TYPE align_type = dynamic_cast<t_gsetins*>(_setkf)->align_type();
    bool ok = false;
    if (align_type == STC_AGN)
    {
        ok = align_coarse(_wm, _vm, _shm);
    }
    else if (align_type == VEL_AGN)
    {
        if (SQRT(SQR(vn(0)) + SQR(vn(1))) > 2)ok = align_vva(vn, _shm);
    }
    else if (align_type == POS_AGN)
    {
        ok = align_pva(pos, _shm);
    }
    if (ok)
        cerr << "\r" << _ins_crt.str_ymdhms("Alignment finished successfully: ") << endl;

    return ok;
}

void great::t_gsinskf::set_Ft()
{
    int iodo = param_of_sins.getParam(_name, par_type::ODO_k, "");
    int iscale = param_of_sins.getParam(_name, par_type::gyro_scale_X, "");

    Ft.block(0, 0, 3, 3) = -t_gbase::askew(sins.eth.weie);   // att-att
    Ft.block(0, 6, 3, 3) = Eigen::Matrix3d::Zero();          // att-pos
    Ft.block(0, 9, 3, 3) = -sins.Ceb;                        // att-eb

    Ft.block(3, 0, 3, 3) = t_gbase::askew(sins.eth.Cen*sins.fn);  // vel-att
    Ft.block(3, 3, 3, 3) = -2 * t_gbase::askew(sins.eth.weie);      // vel-vel
    Ft.block(3, 12, 3, 3) = sins.Ceb;                              // vel-db

    Ft.block(6, 3, 3, 3) = Eigen::Matrix3d::Identity();              // pos-vel

    Ft.block(9, 9, 3, 3) = sins._tauG.array().matrix().asDiagonal();    // eb-eb

    Ft.block(12, 12, 3, 3) = sins._tauA.array().matrix().asDiagonal();    // db-db
}
 
void great::t_gsinskf::set_Hk()
{
    Eigen::Vector3d vartheta = sins.Ceb * t_gbase::askew(lever) * sins.wib + t_gbase::askew(sins.eth.weie) * sins.Ceb * lever;
    switch (Flag)
    {
    case POS_VEL_MEAS:
        Hk.block(0, 3, 3, 3) = Eigen::Matrix3d::Identity();      // vel-vel
        Hk.block(0, 0, 3, 3) = -t_gbase::askew(vartheta);         // vel-att
        Hk.block(0, 9, 3, 3) = -sins.Ceb * t_gbase::askew(lever);  // vel-eb
        Hk.block(3, 6, 3, 3) = Eigen::Matrix3d::Identity();      // pos-pos
        Hk.block(3, 0, 3, 3) = t_gbase::askew(sins.Ceb * lever);   // pos-att
        break;
    case POS_MEAS:
        Hk.block(0, 6, 3, 3) = Eigen::Matrix3d::Identity();
        Hk.block(0, 0, 3, 3) = t_gbase::askew(sins.Ceb * lever);     // pos-att
        break;
    case VEL_MEAS:
        Hk.block(0, 3, 3, 3) = Eigen::Matrix3d::Identity();      // vel-vel
        Hk.block(0, 0, 3, 3) = -t_gbase::askew(vartheta);         // vel-att
        Hk.block(0, 9, 3, 3) = -sins.Ceb * t_gbase::askew(lever);  // vel-eb
        break;
    default:
        break;
    }
}

MOTION_TYPE great::t_gsinskf::motion_state()
{
	if (sins.wnb.norm() < 1.5 * t_gglv::dps && abs(sins.vb(1)) > 5)
		return m_straight;
	else if (sins.an.norm() < 1 && sins.wnb.norm() < 1 * t_gglv::dps && sins.vn.norm() < 0.1)
		return m_static;
	else
		return m_default;
}

MEAS_TYPE great::t_gsinskf::meas_state()
{
    MOTION_TYPE motion = motion_state();
    switch (motion)
    {
    case m_static:
        return ZUPT_MEAS; break;
    case m_straight:
        return NHC_MEAS; break;
    default:
        return NO_MEAS;break;
    }
}

void great::t_gsinskf::feedback()
{
    sins.qeb = sins.qeb - Xk.block(0, 0, 3, 1);
    sins.ve = sins.ve - Xk.block(3, 0, 3, 1);
    sins.pos_ecef = sins.pos_ecef - Xk.block(6, 0, 3, 1);
    sins.eb = sins.eb + Xk.block(9, 0, 3, 1);
    sins.db = sins.db + Xk.block(12, 0, 3, 1);

    int iscale = param_of_sins.getParam(_name, par_type::gyro_scale_X, "");
    int iodo = param_of_sins.getParam(_name, par_type::ODO_k, "");
    int irot = param_of_sins.getParam(_name, par_type::IMU_INST_ATT_X, "");
    if (_shm._imu_scale)
    {
        sins.Kg = sins.Kg - Xk.block(iscale, 0, 3, 1);
        sins.Ka = sins.Ka - Xk.block(iscale + 3, 0, 3, 1);
    }

    sins.pos = Cart2Geod(sins.pos_ecef, false);
    sins.eth.Update(sins.pos, sins.vn);
    sins.vn = sins.eth.Cne * sins.ve;
    sins.qnb = t_gbase::m2qua(sins.eth.Cne) * sins.qeb;
    sins.att = t_gbase::q2att(sins.qnb);
    sins.Ceb = t_gbase::q2mat(sins.qeb);
    sins.Cnb = t_gbase::q2mat(sins.qnb);
}

void great::t_gsinskf::write()
{
    ostringstream os;
    _prt_ins_kml();
    sins.prt_sins(os);

    _fins->write(os.str().c_str(), os.str().size());
    os.str("");

}


void great::t_gsinskf::set_out()
{
    string tmp;
    tmp = dynamic_cast<t_gsetout*>(_setkf)->outputs("ins");
    if (tmp.empty())
    {
        tmp = _name +"result.ins";
    }
    _fins = new t_giof;
    _fins->tsys(t_gtime::GPS);
    if (_name != "") {
        substitute(tmp, "$(rec)", _name, false);
    }
    _fins->mask(tmp);
    _fins->append(dynamic_cast<t_gsetout*>(_setkf)->append());

    ostringstream os;
    sins.prt_header(os);
    _fins->write(os.str().c_str(), os.str().size());
}

bool great::t_gsinskf::_ins_init()
{
    _imudata->interpolate(_resample_intv);
    kftk = sins.t; MeasVel = MeasPos = Eigen::Vector3d::Zero(), Flag = NO_MEAS;
    this->lever = dynamic_cast<t_gsetign*>(_setkf)->lever();
    this->odo_lever = dynamic_cast<t_gsetign*>(_setkf)->odo_lever();
    this->uwb_lever = dynamic_cast<t_gsetign*>(_setkf)->uwb_lever();
    Eigen::VectorXd v(nq), vmax(nq), vmin(nq);

    int iodo = param_of_sins.getParam(_name, par_type::ODO_k, "");
    int iscale = param_of_sins.getParam(_name, par_type::gyro_scale_X, "");
    int irot = param_of_sins.getParam(_name, par_type::IMU_INST_ATT_X, "");

    Eigen::Vector3d angle_tmp, v_tmp, pos_tmp, eb_tmp, db_tmp, gscale_tmp, ascale_tmp, inst_att_tmp;
    double dk_tmp;
    angle_tmp = dynamic_cast<t_gsetign*>(_setkf)->initial_misalignment_std();
    v_tmp = dynamic_cast<t_gsetign*>(_setkf)->initial_vel_std();
    pos_tmp = dynamic_cast<t_gsetign*>(_setkf)->initial_pos_std();
    eb_tmp = dynamic_cast<t_gsetign*>(_setkf)->initial_gyro_std();
    db_tmp = dynamic_cast<t_gsetign*>(_setkf)->initial_acce_std();
    dk_tmp = dynamic_cast<t_gsetign*>(_setkf)->initial_odoscale_std();
    gscale_tmp = dynamic_cast<t_gsetign*>(_setkf)->initial_gyro_scale_std();
    ascale_tmp = dynamic_cast<t_gsetign*>(_setkf)->initial_acce_scale_std();
    inst_att_tmp = dynamic_cast<t_gsetign*>(_setkf)->initial_imu_installation_att_std();
    v << angle_tmp * t_gglv::deg, v_tmp, pos_tmp, eb_tmp* t_gglv::dph, db_tmp* t_gglv::mg;

    v = v.array().abs2();
    Pk = v.array().matrix().asDiagonal();
    Eigen::Matrix3d Patt = Pk.block<3, 3>(0, 0);
    Pk.block<3, 3>(0, 0) = sins.eth.Cen * Patt * sins.eth.Cen.transpose();
   
    if (_spdlogkf)
    {
        SPDLOG_LOGGER_INFO(_spdlogkf, "initial misalignment std(deg): " + dbl2str(angle_tmp(0)) + dbl2str(angle_tmp(1)) + dbl2str(angle_tmp(2)));
        SPDLOG_LOGGER_INFO(_spdlogkf, "initial velocity std(m/s): " + dbl2str(v_tmp(0)) + dbl2str(v_tmp(1)) + dbl2str(v_tmp(2)));
        SPDLOG_LOGGER_INFO(_spdlogkf, "initial posotion std(m): " + dbl2str(pos_tmp(0)) + dbl2str(pos_tmp(1)) + dbl2str(pos_tmp(2)));
        SPDLOG_LOGGER_INFO(_spdlogkf, "initial gyro bias std(deg/h): " + dbl2str(eb_tmp(0)) + dbl2str(eb_tmp(1)) + dbl2str(eb_tmp(2)));
        SPDLOG_LOGGER_INFO(_spdlogkf, "initial acce bias std(mg): "+ dbl2str(db_tmp(0)) + dbl2str(db_tmp(1)) + dbl2str(db_tmp(2)));
    }

    angle_tmp = dynamic_cast<t_gsetign*>(_setkf)->misalignment_psd();
    v_tmp = dynamic_cast<t_gsetign*>(_setkf)->vel_psd();
    pos_tmp = dynamic_cast<t_gsetign*>(_setkf)->pos_psd();
    eb_tmp = dynamic_cast<t_gsetign*>(_setkf)->gyro_psd();
    db_tmp = dynamic_cast<t_gsetign*>(_setkf)->acce_psd();
    gscale_tmp = dynamic_cast<t_gsetign*>(_setkf)->gyro_scale_psd();
    ascale_tmp = dynamic_cast<t_gsetign*>(_setkf)->acce_scale_psd();
    inst_att_tmp = dynamic_cast<t_gsetign*>(_setkf)->imu_inst_rot_psd();
    dk_tmp = dynamic_cast<t_gsetign*>(_setkf)->odo_psd();
    Qt << angle_tmp * t_gglv::dpsh, v_tmp* t_gglv::mgpsHz, pos_tmp* t_gglv::mpsh, eb_tmp* t_gglv::dphpsh, db_tmp* t_gglv::mgpsh;

    Qt = Qt.array().abs2();
    if (_spdlogkf)
    {
        SPDLOG_LOGGER_INFO(_spdlogkf, "misalignment psd(deg/sqrt(h)): " + dbl2str(angle_tmp(0)) + dbl2str(angle_tmp(1)) + dbl2str(angle_tmp(2)));
        SPDLOG_LOGGER_INFO(_spdlogkf, "velocity psd(mg/sqrt(Hz)): " + dbl2str(v_tmp(0)) + dbl2str(v_tmp(1)) + dbl2str(v_tmp(2)));
        SPDLOG_LOGGER_INFO(_spdlogkf, "posotion psd(m/sqrt(h)): " + dbl2str(pos_tmp(0)) + dbl2str(pos_tmp(1)) + dbl2str(pos_tmp(2)));
        SPDLOG_LOGGER_INFO(_spdlogkf, "gyro bias psd(deg/h/sqrt(h)): " + dbl2str(eb_tmp(0)) + dbl2str(eb_tmp(1)) + dbl2str(eb_tmp(2)));
        SPDLOG_LOGGER_INFO(_spdlogkf, "acce bias psd(mg/sqrt(h)): " + dbl2str(db_tmp(0)) + dbl2str(db_tmp(1)) + dbl2str(db_tmp(2)));
    }

    return true;
}


bool great::t_gsinskf::_valid_ins_constraint()
{
    // todo: customizable

    return true;
}

void great::t_gsinskf::time_update(double kfts, double inflation)
{
    set_Ft();
    Eigen::MatrixXd Fk = Eigen::MatrixXd::Identity(nq, nq) + (Ft * kfts);
    Phik = Fk * Phik;
    Xk = Fk * Xk;
    Eigen::MatrixXd Qk = (Qt * kfts * inflation).array().matrix().asDiagonal();
    {
        std::lock_guard<std::mutex> lock(
            process_covariance_diagnostic_mutex);
        ProcessCovarianceDiagnosticState& diagnostic =
            process_covariance_state_unlocked(this);
        if (diagnostic.covariance.rows() != nq ||
            diagnostic.covariance.cols() != nq)
        {
            diagnostic.covariance = Eigen::MatrixXd::Zero(nq, nq);
            diagnostic.dt = 0.0;
            diagnostic.steps = 0;
        }
        diagnostic.covariance =
            Fk * diagnostic.covariance * Fk.transpose() + Qk;
        diagnostic.dt += kfts;
        ++diagnostic.steps;
    }
    Pk = Fk * Pk * (Fk.transpose());  Pk += Qk;
}

int great::t_gsinskf::_meas_update()
{
    if (Flag == NO_MEAS)  return -1;
    Eigen::VectorXd Pxz = Eigen::VectorXd::Zero(Pk.rows()),
        Kk = Eigen::VectorXd::Zero(Pk.rows()), Hi = Xk = Eigen::VectorXd::Zero(Pk.rows());
    this->_set_meas();

    for (int i = 0; i < nr; i++)
    {
        Hi = Hk.block(i, 0, 1, Pk.rows()).transpose();
        Pxz = Pk * Hi;
        double Pz0 = (Hi.transpose() * Pxz), r = Zk(i) - (Hi.transpose() * Xk);
        double Pzz = Pz0 + Rk(i, i);
        Kk = Pxz * (1.0 / Pzz);
        Xk += Kk * r;
        Pk = Pk - Kk * Pxz.transpose();
    }

    if (outlier_detect()!=0)
    {
        Flag = NO_MEAS; return -1;
    }

    t_gbase::symmetry(Pk);
    //set_meas_flag(0);
    if (_spdlogkf)
        SPDLOG_LOGGER_INFO(_spdlogkf, "gintegration: " + _ins_crt.str_ymdhms("LC " + meas2str(Flag) + "/INS filter updating Succeed at epoch: "));

    return 1;
}

bool great::t_gsinskf::_set_meas()
{
    double dt = sins.t - tmeas;
    if (fabs(dt) > _shm.delay) return false;
    _resize();
    Eigen::Vector3d vartheta = sins.Ceb * t_gbase::askew(lever) * sins.wib + t_gbase::askew(sins.eth.weie) * sins.Ceb * lever;
    switch (Flag)
    {
    case POS_VEL_MEAS:
        Zk.block(0, 0, 3, 1) = sins.ve - MeasVel - sins.ae * dt - vartheta;
        Rk.block(0, 0, 3, 3) = _Cov_MeasVn.array().matrix().asDiagonal();
        Zk.block(3, 0, 3, 1) = sins.pos_ecef - MeasPos - sins.ve * dt + sins.Ceb * lever;
        Rk.block(3, 3, 3, 3) = _Cov_MeasPos.array().matrix().asDiagonal();
        measflag = 1;
        MeasVel = MeasPos = Eigen::Vector3d::Zero();
        break;
    case POS_MEAS:
        Zk.block(0, 0, 3, 1) = sins.pos_ecef - MeasPos - sins.ve * dt + sins.Ceb * lever;
        Rk.block(0, 0, 3, 3) = _Cov_MeasPos.array().matrix().asDiagonal();
        measflag = 1;
        MeasVel = MeasPos = Eigen::Vector3d::Zero();
        break;
    case VEL_MEAS:
        Zk.block(0, 0, 3, 1) = sins.ve - MeasVel - sins.ae * dt + sins.Ceb * t_gbase::askew(sins.web) * lever;
        Rk.block(0, 0, 3, 3) = _Cov_MeasVn.array().matrix().asDiagonal();
        measflag = 1;
        MeasVel = MeasPos = Eigen::Vector3d::Zero();
        break;
    case NO_MEAS:
        measflag = 0;
        break;
    default:
        break;
    }
    if (this->measflag)set_Hk();
    return true;
}

bool great::t_gsinskf::_resize()
{
    switch (Flag)
    {
    case POS_VEL_MEAS:
        nr = 6; break;
    case POS_MEAS:
        nr = 3; break;
    case VEL_MEAS:
        nr = 3; break;
    case ZUPT_MEAS:
        nr = 3; break;
    case ZIHR_MEAS:
        nr = 1; break;
    case NHC_MEAS:
        nr = 2; break;
    case ODO_MEAS:
        nr = 1; break;
    case Hgt_MEAS:
        nr = 1; break;
    case YAW_MEAS:
        nr = 1; break;
    case ATT_MEAS:
        nr = 1; break;
    case ZUPT_POS_MEAS:
        nr = 3; break;
    case NO_MEAS:
        nr = 0; break;
    default:
        break;
    }
    Hk = Eigen::MatrixXd::Zero(nr, Pk.rows()); Zk = Eigen::VectorXd::Zero(nr);
    Rk = Eigen::MatrixXd::Zero(nr, nr);
    return true;
}

int great::t_gsinskf::outlier_detect()
{
    double max_set = _map_maxnorm[Flag];

    if (double_eq(max_set, 0))max_set = 3;
    _v_post = Zk - Hk * Xk;
    Eigen::MatrixXd Qz = Hk * Pk * Hk.transpose() + Rk;
    Eigen::VectorXd v_norm = Eigen::VectorXd::Zero(nr);
    for (int i = 0; i < nr; i++)
    {
        v_norm(i) = sqrt(1.0 / Qz(i, i)) * abs(_v_post(i));
    }
    for (int i = 0; i < nr; i++)
    {
        if (v_norm[i] > max_set)
        {
            cout <<  v_norm[i] <<" " << endl;
            SPDLOG_LOGGER_INFO(_spdlogkf, "gintegration: " +
                _ins_crt.str_ymdhms("LC " + meas2str(Flag) + "/INS filter updating Failed (v: " + to_string(v_norm[i]) + ") at epoch : "));
            return 1;
        }
    }
    return 0;
}
