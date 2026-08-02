/**
 * @file         gsetins.cpp
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        control set from XML for INS
 * @version      1.0
 * @date         2025-01-01
 *
 * @copyright Copyright (c) 2025, Wuhan University. All rights reserved.
 *
 */

#include "gset/gsetins.h"
#include "gutils/gsysconv.h"
UNIT_TYPE great::str2Unit(const string& s)
{
    string tmp = s;
    transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
    UNIT_TYPE res = UNDF;
    if (tmp == "RAD")
        res = RAD;
    else if (tmp == "DEG")
        res = DEG;
    else if (tmp == "RPS")
        res = RPS;
    else if (tmp == "DPS")
        res = DPS;
    else if (tmp == "RPH")
        res = RPH;
    else if (tmp == "DPH")
        res = DPH;
    else if (tmp == "MPS")
        res = MPS;
    else if (tmp == "MPS2")
        res = MPS2;
    else
        res = UNDF;
    return res;
}

ALIGN_TYPE great::str2align(const string& s)
{
    string tmp = s;
    transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
    if (tmp == "AUTO")
        return ALIGN_TYPE::AUTO;
    if (tmp == "STATIC")
        return ALIGN_TYPE::STC_AGN;
    if (tmp == "MA")
        return ALIGN_TYPE::MA;
    if (tmp == "POS")
        return ALIGN_TYPE::POS_AGN;
    if (tmp == "VEL")
        return ALIGN_TYPE::VEL_AGN;
    if (tmp == "MULTI_ANT")
        return ALIGN_TYPE::Multi_Ant;
    return STC_AGN;
}

CPS_TYPE great::str2cps(const string& s)
{
    string tmp = s;
    transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
    if (tmp == "POLY")
        return POLY;
    if (tmp == "CONE")
        return CONE;
    if (tmp == "ONE_PLUS_PRE")
        return ONE_PLUS_PRE;
    return POLY;
}

//const t_gglv great::glv = t_gglv();
//t_gglv::t_gglv(double Re, double f, double wie0, double g0)
//{
//}
/*
const t_gglv great::glv = t_gglv();
t_gglv::t_gglv(double Re, double f, double wie0, double g0)
{
    this->Re = Re;
    this->f = f;
    this->wie = wie0;
    this->g0 = g0;
    mg = 1.0e-3 * g0;
    ug = 1.0e-6 * glv.g0;
    mgpsHz = mg / sqrt(1.0);
    ugpsHz = ug / sqrt(1.0);
    mgpsh = mg / sqrt(hur);
    ugpsh = ug / sqrt(hur);
}
*/

great::t_gsetins::t_gsetins() : t_gsetbase()
{
    _set.insert(XMLKEY_INS);
    _freq = 200;
    _ts = 1.0 / _freq;
    _pos = _vel = _att = Eigen::Vector3d::Zero();
    _align_time = 5 * 60;
    _GyroUnit = "DPS";
    _AcceUnit = "MPS2";
    _cps = "POLY";
    _subsample = 1;
    _order = "garfu";
    _out_order = "XYZ_XYZ_PRY_NW";
    _out_freq = 1;
    _int_sec = true;
    _align_type = "SA";
    _start = 0.0;
    _end = _start + 7 * 24 * 60 * 60;
    _acce_bias = Eigen::Vector3d::Zero();
    _gyro_bias = Eigen::Vector3d::Zero();
}

great::t_gsetins::~t_gsetins()
{
}

void great::t_gsetins::check()
{
    _gmutex.lock();
    xml_node parent = _doc.child(XMLKEY_ROOT);
    xml_node node = _default_node(parent, XMLKEY_INS);

    string str = to_string(_pos(0)) + ' ' + to_string(_pos(1)) + ' ' + to_string(_pos(2));
    _default_attr(node, "pos", str);

    str.clear();
    str = "";
    str = to_string(_vel(0)) + ' ' + to_string(_vel(1)) + ' ' + to_string(_vel(2));
    _default_attr(node, "vel", str);

    str.clear();
    str = "";
    str = to_string(_att(0)) + ' ' + to_string(_att(1)) + ' ' + to_string(_att(2));
    _default_attr(node, "att", str);

    _default_attr(node, "freq", 1 / _ts);
    _default_attr(node, "gyro_unit", _GyroUnit);
    _default_attr(node, "acce_unit", _AcceUnit);
    _default_attr(node, "align_type", _align_type);
    _default_attr(node, "align_time", _align_time);
    _default_attr(node, "compensation", _cps);
    _default_attr(node, "subsample", _subsample);
    _default_attr(node, "order", _order);
    _default_attr(node, "out_freq", _out_freq);
    _default_attr(node, "int_sec", _int_sec);
    _default_attr(node, "out_order", _out_order);
    _default_attr(node, "start", _start);
    _default_attr(node, "end", _end);
    _default_attr(node, "imu_scale", false);

    _default_node(node, "motion", "");

    str.clear();
    str = "";
    str = to_string(_gyro_bias(0)) + ' ' + to_string(_gyro_bias(1)) + ' ' + to_string(_gyro_bias(2));
    _default_attr(node, "gyro_bias", str);

    str.clear();
    str = "";
    str = to_string(_acce_bias(0)) + ' ' + to_string(_acce_bias(1)) + ' ' + to_string(_acce_bias(2));
    _default_attr(node, "acce_bias", str);

    _gmutex.unlock();
}

void great::t_gsetins::help()
{
    _gmutex.lock();
    cerr << " <ins \n"
        << "   ts=\"" << _ts << "\" \n"
        << "   gyro_unit=\"" << _GyroUnit << "\" \n"
        << "   acce_unit=\"" << _AcceUnit << "\" \n"

        << "  />\n";

    _gmutex.unlock();
}

string great::t_gsetins::order()
{
    _gmutex.lock();

    string res = _doc.child(XMLKEY_ROOT).child(XMLKEY_INS).child("DataFormat").child("AxisOrder").attribute("Type").value();
    transform(res.begin(), res.end(), res.begin(), ::tolower);

    _gmutex.unlock();

    return res;
}

int great::t_gsetins::resampled_freq()
{
    _gmutex.lock();

    int res = _doc.child(XMLKEY_ROOT).child(XMLKEY_INS).child("DataFormat").child("Resample").attribute("Value").as_int();

    _gmutex.unlock();
    return res;
}

great::UNIT_TYPE great::t_gsetins::GyroUnit()
{
    _gmutex.lock();

    string gyrounit = _doc.child(XMLKEY_ROOT).child(XMLKEY_INS).child("DataFormat").child("GyroUnit").attribute("Type").value();

    _gmutex.unlock();

    return str2Unit(gyrounit);
}

great::UNIT_TYPE great::t_gsetins::AcceUnit()
{
    _gmutex.lock();

    string acceunit = _doc.child(XMLKEY_ROOT).child(XMLKEY_INS).child("DataFormat").child("AcceUnit").attribute("Type").value();

    _gmutex.unlock();

    return str2Unit(acceunit);
}

great::UNIT_TYPE great::t_gsetins::AttUnit()
{
    _gmutex.lock();

    string attunit = _doc.child(XMLKEY_ROOT).child(XMLKEY_INS).child("DataFormat").child("AttUnit").attribute("Type").value();

    _gmutex.unlock();

    return str2Unit(attunit);
}

int great::t_gsetins::freq()
{
    _gmutex.lock();

    int res = _doc.child(XMLKEY_ROOT).child(XMLKEY_INS).child("DataFormat").child("Frequency").attribute("Value").as_int();

    _gmutex.unlock();
    return res;
}

double great::t_gsetins::ts()
{
    _gmutex.lock();

    _freq = _doc.child(XMLKEY_ROOT).child(XMLKEY_INS).child("DataFormat").child("Frequency").attribute("Value").as_int();
    double res = 1.0 / _freq;

    _gmutex.unlock();
    return res;
}


great::UNIT_TYPE great::t_gsetins::MagUnit()
{
    _gmutex.lock();

    string acceunit = _doc.child(XMLKEY_ROOT).child(XMLKEY_INS).attribute("mag_unit").value();

    _gmutex.unlock();

    return str2Unit(acceunit);
}

Eigen::Vector3d great::t_gsetins::pos()
{
    _gmutex.lock();
    xml_node PosNode = _doc.child(XMLKEY_ROOT).child(XMLKEY_INS).child("InitialStates").child("Position");
    string type = PosNode.attribute("Type").value();
    double X = PosNode.attribute("X").as_double();
    double Y = PosNode.attribute("Y").as_double();
    double Z = PosNode.attribute("Z").as_double();
    if (type == "Cartesian")
    {
        // X; Y; Z;
    }
    else if (type == "Geodetic")
    {
        double ell[3] = { X, Y, Z }, XYZ[3];
        gnut::ell2xyz(ell, XYZ, true);
        X = XYZ[0]; Y = XYZ[1]; Z = XYZ[2];
    }
    else
    {
        X = Y = Z = 0.0;
    }
    _gmutex.unlock();
    return Eigen::Vector3d(X, Y, Z);
}

Eigen::Vector3d great::t_gsetins::vel()
{
    _gmutex.lock();
    xml_node VelNode = _doc.child(XMLKEY_ROOT).child(XMLKEY_INS).child("InitialStates").child("Velocity");
    string type = VelNode.attribute("Type").value();
    double X = VelNode.attribute("X").as_double();
    double Y = VelNode.attribute("Y").as_double();
    double Z = VelNode.attribute("Z").as_double();
    if (type == "Cartesian")
    {
        // X; Y; Z;
    }
    else
    {
        X = Y = Z = 0.0;
    }
    _gmutex.unlock();
    return Eigen::Vector3d(X, Y, Z);
}

Eigen::Vector3d great::t_gsetins::att()
{
    _gmutex.lock();
    xml_node AttNode = _doc.child(XMLKEY_ROOT).child(XMLKEY_INS).child("InitialStates").child("Attitude");
    string type = AttNode.attribute("Type").value();
    double X = AttNode.attribute("Pitch").as_double();
    double Y = AttNode.attribute("Roll").as_double();
    double Z = AttNode.attribute("Yaw").as_double();
    if (type == "ON")
    {
        // X; Y; Z;
    }
    else
    {
        X = Y = Z = 0.0;
    }
    _gmutex.unlock();
    return Eigen::Vector3d(X, Y, Z) * t_gglv::deg;
}


Eigen::Vector3d great::t_gsetins::gyro_bias()
{
    _gmutex.lock();
    xml_node GyroBiasNode = _doc.child(XMLKEY_ROOT).child(XMLKEY_INS).child("InitialStates").child("GyroBias");
    string type = GyroBiasNode.attribute("Type").value();
    double X = GyroBiasNode.attribute("X").as_double();
    double Y = GyroBiasNode.attribute("Y").as_double();
    double Z = GyroBiasNode.attribute("Z").as_double();
    if (type == "ON")
    {
        // X; Y; Z;
    }
    else if (type == "OFF")
    {
        X = Y = Z = 0.0;
    }
    else
    {
        //todo: MEAN
    }
    _gmutex.unlock();
    return Eigen::Vector3d(X, Y, Z);
}

Eigen::Vector3d great::t_gsetins::imu_installation_rotation()
{
    _gmutex.lock();
    xml_node InstallationNode = _doc.child(XMLKEY_ROOT).child(XMLKEY_INS).child("Installation");
    string type = InstallationNode.attribute("Type").value();
    if (type == "OFF") { _gmutex.unlock(); return Eigen::Vector3d::Zero(); }
    double X = InstallationNode.child("Rotation").attribute("Pitch").as_double();
    double Y = InstallationNode.child("Rotation").attribute("Roll").as_double();
    double Z = InstallationNode.child("Rotation").attribute("Yaw").as_double();
    _gmutex.unlock();
    return Eigen::Vector3d(X, Y, Z) * t_gglv::deg;
}

Eigen::Vector3d great::t_gsetins::imu_installation_translation()
{
    _gmutex.lock();
    xml_node InstallationNode = _doc.child(XMLKEY_ROOT).child(XMLKEY_INS).child("Installation");
    string type = InstallationNode.attribute("Type").value();
    if (type == "OFF") { _gmutex.unlock(); return Eigen::Vector3d::Zero(); }
    double X = InstallationNode.child("Lever").attribute("X").as_double();
    double Y = InstallationNode.child("Lever").attribute("Y").as_double();
    double Z = InstallationNode.child("Lever").attribute("Z").as_double();
    _gmutex.unlock();
    return Eigen::Vector3d(X, Y, Z);
}

Eigen::Vector3d great::t_gsetins::acce_bias()
{
    _gmutex.lock();
    xml_node AcceBiasNode = _doc.child(XMLKEY_ROOT).child(XMLKEY_INS).child("InitialStates").child("AcceBias");
    string type = AcceBiasNode.attribute("Type").value();
    double X = AcceBiasNode.attribute("X").as_double();
    double Y = AcceBiasNode.attribute("Y").as_double();
    double Z = AcceBiasNode.attribute("Z").as_double();
    if (type == "ON")
    {
        // X; Y; Z;
    }
    else if (type == "OFF")
    {
        X = Y = Z = 0.0;
    }
    else
    {
        //todo: MEAN
    }

    _gmutex.unlock();
    return Eigen::Vector3d(X, Y, Z);
}


int great::t_gsetins::cps()
{
    _gmutex.lock();
    string x = _doc.child(XMLKEY_ROOT).child(XMLKEY_INS).child("Compensation").attribute("Type").as_string();
    _gmutex.unlock();
    return str2cps(x);
}

int great::t_gsetins::subsample()
{
    _gmutex.lock();
    int x = _doc.child(XMLKEY_ROOT).child(XMLKEY_INS).child("Compensation").attribute("Subsample").as_int();
    _gmutex.unlock();
    if (!x)x = 1;
    return x;
}

double great::t_gsetins::start()
{
    _gmutex.lock();

    double res = _doc.child(XMLKEY_ROOT).child(XMLKEY_INS).child("ProcTime").attribute("Start").as_double();

    _gmutex.unlock();
    return res;
}

double great::t_gsetins::end()
{
    _gmutex.lock();

    double res = _doc.child(XMLKEY_ROOT).child(XMLKEY_INS).child("ProcTime").attribute("End").as_double();

    _gmutex.unlock();
    return res;
}


bool great::t_gsetins::align()
{
    _gmutex.lock();
    string res = _doc.child(XMLKEY_ROOT).child(XMLKEY_INS).child("Alignment").attribute("Type").value();
    bool x = false;
    if (res != "OFF")x = true;
    _gmutex.unlock();
    return x;
}

ALIGN_TYPE great::t_gsetins::align_type()
{
    _gmutex.lock();

    string res = _doc.child(XMLKEY_ROOT).child(XMLKEY_INS).child("Alignment").attribute("Type").value();

    _gmutex.unlock();

    return str2align(res);
}

double great::t_gsetins::align_time()
{
    _gmutex.lock();

    double res = _doc.child(XMLKEY_ROOT).child(XMLKEY_INS).child("Alignment").child("CoarseAlignTime").attribute("Value").as_double();

    _gmutex.unlock();
    return res;
}

double great::t_gsetins::pos_dist()
{
    _gmutex.lock();

    double res = _doc.child(XMLKEY_ROOT).child(XMLKEY_INS).child("Alignment").child("PositionVector").attribute("Value").as_double();

    _gmutex.unlock();
    return res;
}

double great::t_gsetins::vel_norm()
{
    _gmutex.lock();

    double res = _doc.child(XMLKEY_ROOT).child(XMLKEY_INS).child("Alignment").child("VelocityVector").attribute("Value").as_double();

    _gmutex.unlock();
    return res;
}

bool great::t_gsetins::out_intsec()
{
    _gmutex.lock();
    bool x = _doc.child(XMLKEY_ROOT).child(XMLKEY_INS).attribute("int_sec").as_bool();
    _gmutex.unlock();
    return x;
}

int great::t_gsetins::out_freq()
{
    _gmutex.lock();

    int res = _doc.child(XMLKEY_ROOT).child(XMLKEY_INS).attribute("out_freq").as_double();

    _gmutex.unlock();
    return res;
}

vector<string> great::t_gsetins::motion()
{
    _gmutex.lock();

    vector<string> motion_strvec = t_gsetbase::_vecval(XMLKEY_INS, "Motion");

    _gmutex.unlock();
    return motion_strvec;
}

/********************************** IMUErrorModel *****************************************/
map<IMU_TYPE, ErrorModel> great::IMUErrorModels()
{
    IMU_TYPE imu_type; map<IMU_TYPE, ErrorModel> map_imu_error_models;


    imu_type = IMU_TYPE::NovAtel_SPAN_FSAS;
    map_imu_error_models.insert(make_pair(imu_type, ErrorModel()));
    map_imu_error_models[imu_type].AttInitialSTD = Eigen::Vector3d(0.5, 0.5, 5);
    map_imu_error_models[imu_type].VelInitialSTD = Eigen::Vector3d(0.5, 0.5, 0.5);
    map_imu_error_models[imu_type].PosInitialSTD = Eigen::Vector3d(3, 3, 3);

    map_imu_error_models[imu_type].GyroBiasInitialSTD = Eigen::Vector3d(0.324, 0.324, 0.324);
    map_imu_error_models[imu_type].AcceBiasInitialSTD = Eigen::Vector3d(2.000e-002, 2.000e-002, 5.000e-002) / t_gglv::mg;

    map_imu_error_models[imu_type].AttProcNoisePSD = Eigen::Vector3d(6.350, 6.350, 6.350) * t_gglv::mpsh;
    map_imu_error_models[imu_type].VelProcNoisePSD = Eigen::Vector3d(2.500e-007, 2.500e-007, 8.100e-007) / t_gglv::mgpsHz;
    map_imu_error_models[imu_type].PosProcNoisePSD = Eigen::Vector3d(1.000e-006, 1.000e-006, 1.000e-006) / t_gglv::mpsh;
    map_imu_error_models[imu_type].GyroBiasProcNoisePSD = Eigen::Vector3d(3.240e-012, 3.240e-012, 3.240e-012) / t_gglv::mpsh;
    map_imu_error_models[imu_type].AcceBiasProcNoisePSD = Eigen::Vector3d(9.000e-010, 9.000e-010, 9.000e-010) / t_gglv::mgpsh;


    imu_type = IMU_TYPE::NovAtel_SPAN_CPT;
    map_imu_error_models.insert(make_pair(imu_type, ErrorModel()));
    map_imu_error_models[imu_type].AttInitialSTD = Eigen::Vector3d(0.5, 0.5, 5);
    map_imu_error_models[imu_type].VelInitialSTD = Eigen::Vector3d(0.5, 0.5, 0.5);
    map_imu_error_models[imu_type].PosInitialSTD = Eigen::Vector3d(3, 3, 3);

    map_imu_error_models[imu_type].GyroBiasInitialSTD = Eigen::Vector3d(1.99999800e+001, 1.99999800e+001, 1.99999800e+001);
    map_imu_error_models[imu_type].AcceBiasInitialSTD = Eigen::Vector3d(5.00000000e-002, 5.00000000e-002, 5.00000000e-002) / t_gglv::mg;

    map_imu_error_models[imu_type].AttProcNoisePSD = Eigen::Vector3d(9.80100000e+001, 9.80100000e+001, 9.80100000e+001) * t_gglv::mpsh;
    map_imu_error_models[imu_type].VelProcNoisePSD = Eigen::Vector3d(4.68722500e-008, 4.68722500e-008, 4.68722500e-008) / t_gglv::mgpsHz;
    map_imu_error_models[imu_type].PosProcNoisePSD = Eigen::Vector3d(1.00000000e-008, 1.00000000e-008, 1.00000000e-008) / t_gglv::mpsh;
    map_imu_error_models[imu_type].GyroBiasProcNoisePSD = Eigen::Vector3d(7.71595062e-012, 7.71595062e-012, 7.71595062e-012) / t_gglv::mpsh;
    map_imu_error_models[imu_type].AcceBiasProcNoisePSD = Eigen::Vector3d(4.38906250e-009, 4.38906250e-009, 4.38906250e-009) / t_gglv::mgpsh;


    imu_type = IMU_TYPE::NovAtel_SPAN_uIRS;
    map_imu_error_models.insert(make_pair(imu_type, ErrorModel()));
    map_imu_error_models[imu_type].AttInitialSTD = Eigen::Vector3d(0.5, 0.5, 5);
    map_imu_error_models[imu_type].VelInitialSTD = Eigen::Vector3d(0.5, 0.5, 0.5);
    map_imu_error_models[imu_type].PosInitialSTD = Eigen::Vector3d(3, 3, 3);

    map_imu_error_models[imu_type].GyroBiasInitialSTD = Eigen::Vector3d(3.60000000e-003, 3.60000000e-003, 3.60000000e-003);
    map_imu_error_models[imu_type].AcceBiasInitialSTD = Eigen::Vector3d(3.00000000e-004, 3.00000000e-004, 3.00000000e-004) / t_gglv::mg;

    map_imu_error_models[imu_type].AttProcNoisePSD = Eigen::Vector3d(1.00000080e-001, 1.00000080e-001, 1.00000080e-001) * t_gglv::mpsh;
    map_imu_error_models[imu_type].VelProcNoisePSD = Eigen::Vector3d(1.00000000e-006, 1.00000000e-006, 1.00000000e-006) / t_gglv::mgpsHz;
    map_imu_error_models[imu_type].PosProcNoisePSD = Eigen::Vector3d(1.00000000e-006, 1.00000000e-006, 1.00000000e-006) / t_gglv::mpsh;
    map_imu_error_models[imu_type].GyroBiasProcNoisePSD = Eigen::Vector3d(9.99998950e-012, 9.99998950e-012, 9.99998950e-012) / t_gglv::mpsh;
    map_imu_error_models[imu_type].AcceBiasProcNoisePSD = Eigen::Vector3d(9.99799104e-014, 9.99799104e-014, 9.99799104e-014) / t_gglv::mgpsh;


    imu_type = IMU_TYPE::NovAtel_SPAN_LCI100C;
    map_imu_error_models.insert(make_pair(imu_type, ErrorModel()));
    map_imu_error_models[imu_type].AttInitialSTD = Eigen::Vector3d(0.5, 0.5, 5);
    map_imu_error_models[imu_type].VelInitialSTD = Eigen::Vector3d(0.5, 0.5, 0.5);
    map_imu_error_models[imu_type].PosInitialSTD = Eigen::Vector3d(3, 3, 3);

    map_imu_error_models[imu_type].GyroBiasInitialSTD = Eigen::Vector3d(3.24000000e-001, 3.24000000e-001, 3.24000000e-001);
    map_imu_error_models[imu_type].AcceBiasInitialSTD = Eigen::Vector3d(2.00000000e-002, 2.00000000e-002, 2.00000000e-002) / t_gglv::mg;

    map_imu_error_models[imu_type].AttProcNoisePSD = Eigen::Vector3d(1.00000160e+000, 1.00000160e+000, 1.00000160e+000) * t_gglv::mpsh;
    map_imu_error_models[imu_type].VelProcNoisePSD = Eigen::Vector3d(1.00000000e-006, 1.00000000e-006, 1.00000000e-006) / t_gglv::mgpsHz;
    map_imu_error_models[imu_type].PosProcNoisePSD = Eigen::Vector3d(1.00000000e-008, 1.00000000e-008, 1.00000000e-008) / t_gglv::mpsh;
    map_imu_error_models[imu_type].GyroBiasProcNoisePSD = Eigen::Vector3d(9.99998950e-016, 9.99998950e-016, 9.99998950e-016) / t_gglv::mpsh;
    map_imu_error_models[imu_type].AcceBiasProcNoisePSD = Eigen::Vector3d(2.50000000e-009, 2.50000000e-009, 2.50000000e-009) / t_gglv::mgpsh;

    imu_type = IMU_TYPE::Navigation_Grade;
    map_imu_error_models.insert(make_pair(imu_type, ErrorModel()));
    map_imu_error_models[imu_type].AttInitialSTD = Eigen::Vector3d(0.5, 0.5, 5);
    map_imu_error_models[imu_type].VelInitialSTD = Eigen::Vector3d(0.5, 0.5, 0.5);
    map_imu_error_models[imu_type].PosInitialSTD = Eigen::Vector3d(3, 3, 3);

    map_imu_error_models[imu_type].GyroBiasInitialSTD = Eigen::Vector3d(0.01, 0.01, 0.01);
    map_imu_error_models[imu_type].AcceBiasInitialSTD = Eigen::Vector3d(0.0003, 0.0003, 0.0003) / t_gglv::mg;

    map_imu_error_models[imu_type].AttProcNoisePSD = Eigen::Vector3d(1.0, 1.0, 1.0) * t_gglv::mpsh;
    map_imu_error_models[imu_type].VelProcNoisePSD = Eigen::Vector3d(1.0e-6, 1.0e-6, 1.0e-6) / t_gglv::mgpsHz;
    map_imu_error_models[imu_type].PosProcNoisePSD = Eigen::Vector3d(1.00000000e-006, 1.00000000e-006, 1.00000000e-006) / t_gglv::mpsh;
    map_imu_error_models[imu_type].GyroBiasProcNoisePSD = Eigen::Vector3d(1.0e-11, 1.0e-11, 1.0e-11) / t_gglv::mpsh;
    map_imu_error_models[imu_type].AcceBiasProcNoisePSD = Eigen::Vector3d(1.0e-11, 1.0e-11, 1.0e-11) / t_gglv::mgpsh;


    imu_type = IMU_TYPE::Tactical_Grade;
    map_imu_error_models.insert(make_pair(imu_type, ErrorModel()));
    map_imu_error_models[imu_type].AttInitialSTD = Eigen::Vector3d(0.5, 0.5, 5);
    map_imu_error_models[imu_type].VelInitialSTD = Eigen::Vector3d(0.5, 0.5, 0.5);
    map_imu_error_models[imu_type].PosInitialSTD = Eigen::Vector3d(3, 3, 3);

    map_imu_error_models[imu_type].GyroBiasInitialSTD = Eigen::Vector3d(0.1, 0.1, 0.1);
    map_imu_error_models[imu_type].AcceBiasInitialSTD = Eigen::Vector3d(0.003, 0.003, 0.003) / t_gglv::mg;

    map_imu_error_models[imu_type].AttProcNoisePSD = Eigen::Vector3d(20.0, 20.0, 20.0) * t_gglv::mpsh;
    map_imu_error_models[imu_type].VelProcNoisePSD = Eigen::Vector3d(1.0e-6, 1.0e-6, 1.0e-6) / t_gglv::mgpsHz;
    map_imu_error_models[imu_type].PosProcNoisePSD = Eigen::Vector3d(1.00000000e-006, 1.00000000e-006, 1.00000000e-006) / t_gglv::mpsh;
    map_imu_error_models[imu_type].GyroBiasProcNoisePSD = Eigen::Vector3d(1.0e-7, 1.0e-7, 1.0e-7) / t_gglv::mpsh;
    map_imu_error_models[imu_type].AcceBiasProcNoisePSD = Eigen::Vector3d(1.0e-7, 1.0e-7, 1.0e-7) / t_gglv::mgpsh;


    imu_type = IMU_TYPE::MEMS_Grade;
    map_imu_error_models.insert(make_pair(imu_type, ErrorModel()));
    map_imu_error_models[imu_type].AttInitialSTD = Eigen::Vector3d(1, 1, 10);
    map_imu_error_models[imu_type].VelInitialSTD = Eigen::Vector3d(0.5, 0.5, 0.5);
    map_imu_error_models[imu_type].PosInitialSTD = Eigen::Vector3d(3, 3, 3);

    map_imu_error_models[imu_type].GyroBiasInitialSTD = Eigen::Vector3d(100, 100, 100);
    map_imu_error_models[imu_type].AcceBiasInitialSTD = Eigen::Vector3d(0.3, 0.3, 0.3) / t_gglv::mg;

    map_imu_error_models[imu_type].AttProcNoisePSD = Eigen::Vector3d(100000, 100000, 100000) * t_gglv::mpsh;
    map_imu_error_models[imu_type].VelProcNoisePSD = Eigen::Vector3d(1.0e-6, 1.0e-6, 1.0e-6) / t_gglv::mgpsHz;
    map_imu_error_models[imu_type].PosProcNoisePSD = Eigen::Vector3d(1.00000000e-006, 1.00000000e-006, 1.00000000e-006) / t_gglv::mpsh;
    map_imu_error_models[imu_type].GyroBiasProcNoisePSD = Eigen::Vector3d(1.0e-4, 1.0e-4, 1.0e-4) / t_gglv::mpsh;
    map_imu_error_models[imu_type].AcceBiasProcNoisePSD = Eigen::Vector3d(1.0e-4, 1.0e-4, 1.0e-4) / t_gglv::mgpsh;


    imu_type = IMU_TYPE::NovAtel_SPAN_ADIS16488;
    map_imu_error_models.insert(make_pair(imu_type, ErrorModel()));
    map_imu_error_models[imu_type].AttInitialSTD = Eigen::Vector3d(0.5, 0.5, 5);
    map_imu_error_models[imu_type].VelInitialSTD = Eigen::Vector3d(0.5, 0.5, 0.5);
    map_imu_error_models[imu_type].PosInitialSTD = Eigen::Vector3d(3, 3, 3);

    map_imu_error_models[imu_type].GyroBiasInitialSTD = Eigen::Vector3d(7.20000000e+002, 7.20000000e+002, 7.20000000e+002);
    map_imu_error_models[imu_type].AcceBiasInitialSTD = Eigen::Vector3d(5.00000000e-002, 5.00000000e-002, 5.00000000e-002) / t_gglv::mg;

    map_imu_error_models[imu_type].AttProcNoisePSD = Eigen::Vector3d(3.24000000e+002, 3.24000000e+002, 3.24000000e+002) * t_gglv::mpsh;
    map_imu_error_models[imu_type].VelProcNoisePSD = Eigen::Vector3d(4.68722500e-006, 4.68722500e-006, 4.68722500e-006) / t_gglv::mgpsHz;
    map_imu_error_models[imu_type].PosProcNoisePSD = Eigen::Vector3d(1.00000000e-004, 1.00000000e-004, 1.00000000e-004) / t_gglv::mpsh;
    map_imu_error_models[imu_type].GyroBiasProcNoisePSD = Eigen::Vector3d(9.56014618e-003, 9.56014618e-003, 9.56014618e-003) / t_gglv::mpsh;
    map_imu_error_models[imu_type].AcceBiasProcNoisePSD = Eigen::Vector3d(3.16406250e-009, 3.16406250e-009, 3.16406250e-009) / t_gglv::mgpsh;


    imu_type = IMU_TYPE::ADIS16470;
    map_imu_error_models.insert(make_pair(imu_type, ErrorModel()));
    map_imu_error_models[imu_type].AttInitialSTD = Eigen::Vector3d(0.5, 0.5, 5);
    map_imu_error_models[imu_type].VelInitialSTD = Eigen::Vector3d(0.5, 0.5, 0.5);
    map_imu_error_models[imu_type].PosInitialSTD = Eigen::Vector3d(3, 3, 3);

    map_imu_error_models[imu_type].GyroBiasInitialSTD = Eigen::Vector3d(1.20000000e+002, 1.20000000e+002, 1.20000000e+002);
    map_imu_error_models[imu_type].AcceBiasInitialSTD = Eigen::Vector3d(5.00000000e-002, 5.00000000e-002, 5.00000000e-002) / t_gglv::mg;

    map_imu_error_models[imu_type].AttProcNoisePSD = Eigen::Vector3d(1.74000000e+001, 1.74000000e+001, 1.74000000e+001) * t_gglv::mpsh;
    map_imu_error_models[imu_type].VelProcNoisePSD = Eigen::Vector3d(0.68722500e-003, 0.68722500e-003, 0.68722500e-003) / t_gglv::mgpsHz;
    map_imu_error_models[imu_type].PosProcNoisePSD = Eigen::Vector3d(1.00000000e-004, 1.00000000e-004, 1.00000000e-004) / t_gglv::mpsh;
    map_imu_error_models[imu_type].GyroBiasProcNoisePSD = Eigen::Vector3d(8.56014618e-002, 8.56014618e-002, 8.56014618e-002) / t_gglv::mpsh;
    map_imu_error_models[imu_type].AcceBiasProcNoisePSD = Eigen::Vector3d(3.16406250e-004, 3.16406250e-004, 3.16406250e-004) / t_gglv::mgpsh;


    imu_type = IMU_TYPE::StarNeto;
    map_imu_error_models.insert(make_pair(imu_type, ErrorModel()));
    map_imu_error_models[imu_type].AttInitialSTD = Eigen::Vector3d(0.5, 0.5, 5);
    map_imu_error_models[imu_type].VelInitialSTD = Eigen::Vector3d(0.5, 0.5, 0.5);
    map_imu_error_models[imu_type].PosInitialSTD = Eigen::Vector3d(3, 3, 3);

    map_imu_error_models[imu_type].GyroBiasInitialSTD = Eigen::Vector3d(1, 1, 1);
    map_imu_error_models[imu_type].AcceBiasInitialSTD = Eigen::Vector3d(1.000e-002, 1.000e-002, 1.000e-002) / t_gglv::mg;

    map_imu_error_models[imu_type].AttProcNoisePSD = Eigen::Vector3d(16.350, 16.350, 16.350) * t_gglv::mpsh;
    map_imu_error_models[imu_type].VelProcNoisePSD = Eigen::Vector3d(1.600e-003, 1.600e-003, 1.600e-003) / t_gglv::mgpsHz;
    map_imu_error_models[imu_type].PosProcNoisePSD = Eigen::Vector3d(1.000e-006, 1.000e-006, 1.000e-006) / t_gglv::mpsh;
    map_imu_error_models[imu_type].GyroBiasProcNoisePSD = Eigen::Vector3d(3.240e-012, 3.240e-012, 3.240e-012) / t_gglv::mpsh;
    map_imu_error_models[imu_type].AcceBiasProcNoisePSD = Eigen::Vector3d(9.000e-010, 9.000e-010, 9.000e-010) / t_gglv::mgpsh;

    return map_imu_error_models;
}