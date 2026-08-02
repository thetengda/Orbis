/**
 * @file         gsetign.cpp
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        control set from XML for MSF
 * @version      1.0
 * @date         2025-01-01
 *
 * @copyright Copyright (c) 2025, Wuhan University. All rights reserved.
 *
 */

#include "gset/gsetign.h"

using namespace Eigen;

string great::meas2str(MEAS_TYPE type)
{
    string res;

    switch (type)
    {
    case great::GNSS_MEAS:
        res = "GNSS";
        break;
    case great::POS_MEAS:
        res = "GNSS_POS";
        break;
    case great::VEL_MEAS:
        res = "GNSS_VEL";
        break;
    case great::POS_VEL_MEAS:
        res = "GNSS";
        break;
    case great::MOTION_MEAS:
        break;
    case great::ZUPT_MEAS:
        res = "ZUPT";
        break;
    case great::ZIHR_MEAS:
        res = "ZIHR";
        break;
    case great::NHC_MEAS:
        res = "NHC";
        break;
    case great::YAW_MEAS:
        res = "YAW";
        break;
    case great::ODO_MEAS:
        res = "ODO";
        break;
    case great::Hgt_MEAS:
        res = "Hgt";
        break;
    case great::ATT_MEAS:
        res = "ATT";
        break;
    case great::ZUPT_POS_MEAS:
        res = "ZUPT";
        break;
    case great::NO_MEAS:
        break;
    case great::VIS_MEAS:
        res = "VIS";
        break;
    case great::DEFAULT_MEAS:
        break;
    case great::LIDAR_MEAS:
        res = "LIDAR";
        break;
    case great::UWB_MEAS:
        res = "UWB";
        break;
    default:
        break;
    }
    return res;
}

IGN_TYPE great::str2ign(const string& s)
{
    string tmp = s;
    transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
    if (tmp == "PURE_INS")
        return PURE_INS;
    if (tmp == "LCI")
        return LCI;
    if (tmp == "TCI")
        return TCI;
    if (tmp == "STCI")
        return STCI;
    if (tmp == "VIO")
        return VIO;
    if (tmp == "VIS_LCI")
        return VIS_LCI;
    if (tmp == "VIS_TCI")
        return VIS_TCI;
    if (tmp == "VIS_STCI")
        return VIS_STCI;
    if (tmp == "LIO")
        return LIO;
    if (tmp == "LIDAR_LCI")
        return LIDAR_LCI;
    if (tmp == "LIDAR_TCI")
        return LIDAR_TCI;
    if (tmp == "LIDAR_STCI")
        return LIDAR_STCI;
    if (tmp == "VLO")
        return VLO;
    if (tmp == "MULTIGN_LCI")
        return MULTIGN_LCI;
    if (tmp == "MULTIGN_TCI")
        return MULTIGN_TCI;
    if (tmp == "MULTIGN_STCI")
        return MULTIGN_STCI;
    if (tmp == "UWB_LCI")
        return UWB_LCI;
    if (tmp == "UWB_TCI")
        return UWB_TCI;
    return PURE_INS;
}

IMU_TYPE great::str2imu(const string& s)
{
    string tmp = s;
    //transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
    if (tmp == "SPAN FSAS")
        return IMU_TYPE::NovAtel_SPAN_FSAS;
    if (tmp == "SPAN CPT")
        return IMU_TYPE::NovAtel_SPAN_CPT;
    if (tmp == "SPAN uIRS")
        return IMU_TYPE::NovAtel_SPAN_uIRS;
    if (tmp == "SPAN LCI100C")
        return IMU_TYPE::NovAtel_SPAN_LCI100C;
    if (tmp == "Navigation Grade")
        return IMU_TYPE::Navigation_Grade;
    if (tmp == "Tactical Grade")
        return IMU_TYPE::Tactical_Grade;
    if (tmp == "MEMS Grade")
        return IMU_TYPE::MEMS_Grade;
    if (tmp == "ADIS 16470")
        return IMU_TYPE::ADIS16470;
    if (tmp == "ADIS 16488")
        return IMU_TYPE::NovAtel_SPAN_ADIS16488;
    if (tmp == "StarNeto")
        return IMU_TYPE::StarNeto;
    if (tmp == "Customize")
        return IMU_TYPE::Customize;
    if (tmp == "MEMS Grade")
        return IMU_TYPE::MEMS_Grade;
    if (tmp == "MEMS Grade")
        return IMU_TYPE::MEMS_Grade;
    return IMU_TYPE::Customize;
}

great::t_gsetign::t_gsetign() : t_gsetbase()
{
    _set.insert(XMLKEY_IGN);

    _initial_misalignment_std = Vector3d(1, 1, 10);
    _initial_vel_std = 0.1 * Vector3d::Ones();
    _initial_pos_std = 3 * Vector3d::Ones();
    _initial_gyro_std = 10 * Vector3d::Ones();
    _initial_acce_std = 5 * Vector3d::Ones();
    _initial_gyro_scale_std = 0.01 * Vector3d::Ones();
    _initial_acce_scale_std = 0.01 * Vector3d::Ones();
    _initial_odoscale_std = 0.1;

    _min_misalignment_std = 0.1 * Vector3d::Ones();
    _min_vel_std = 0.005 * Vector3d::Ones();
    _min_pos_std = 0.01 * Vector3d::Ones();
    _min_gyro_std = 0.1 * Vector3d::Ones();
    _min_acce_std = 0.1 * Vector3d::Ones();
    _odo_std = 0.5;

    _misalignment_psd = 1 * Vector3d::Ones();
    _vel_psd = 0.1 * Vector3d::Ones();
    _pos_psd = 0 * Vector3d::Ones();
    _gyro_psd = 0 * Vector3d::Ones();
    _acce_psd = 0 * Vector3d::Ones();
    _gyro_scale_psd = 0 * Vector3d::Ones();
    _acce_scale_psd = 0 * Vector3d::Ones();
    _odo_psd = 0.1;

    _meas_vel_noise = Vector3d(0.5, 0.5, 0.5);
    _meas_pos_noise = Vector3d(0.5, 0.5, 0.5);

    _nq = 15;
    _nr = 6;
    _lever = _odo_lever = _uwb_lever = Vector3d::Zero();
    _max_pdop = 8.0;
    _min_sat = 5;
    _order = "XYZXYZ";
    _TS = 1.0;
    _delay_t = 0.0;
    _delay_odo = 0.001;
    _fltmode = FORWARD;
    _odo = "OFF";
    _NHC = _ZUPT = false;
    //_imu_type = IMU_TYPE::StarNeto;
}

great::t_gsetign::~t_gsetign()
{
}

void great::t_gsetign::check()
{
    _gmutex.lock();

    _imu_type = str2imu(_doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).child("Estimator").child("IMUErrorModel").attribute("Type").value());
    _map_imu_error_models = IMUErrorModels();

    xml_node node = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN);

    string str = to_string(_initial_misalignment_std(0)) + ' ' + to_string(_initial_misalignment_std(1)) + ' ' + to_string(_initial_misalignment_std(2));
    _default_attr(node, "initial_misalignment_std", str);

    str.clear();
    str = "";
    str = to_string(_initial_vel_std(0)) + ' ' + to_string(_initial_vel_std(1)) + ' ' + to_string(_initial_vel_std(2));
    _default_attr(node, "initial_vel_std", str);

    str.clear();
    str = "";
    str = to_string(_initial_pos_std(0)) + ' ' + to_string(_initial_pos_std(1)) + ' ' + to_string(_initial_pos_std(2));
    _default_attr(node, "initial_pos_std", str);

    _default_attr(node, "initial_gyro_std", _initial_gyro_std(0));
    _default_attr(node, "initial_acce_std", _initial_acce_std(0));
    _default_attr(node, "initial_gyro_scale_std", _initial_gyro_scale_std(0));
    _default_attr(node, "initial_acce_scale_std", _initial_acce_scale_std(0));
    _default_attr(node, "initial_odo_k_std", _initial_odoscale_std);
    _default_attr(node, "odo_std", _odo_std);
    _default_attr(node, "NHC_std", _odo_std);
    _default_attr(node, "ZUPT_std", _odo_std);
    _default_attr(node, "ZAUPT_std", _odo_std);

    _default_attr(node, "min_misalignment_std", _min_misalignment_std(0));
    _default_attr(node, "min_vel_std", _min_vel_std(0));
    _default_attr(node, "min_pos_std", _min_pos_std(0));
    _default_attr(node, "min_gyro_std", _min_gyro_std(0));
    _default_attr(node, "min_acce_std", _min_acce_std(0));

    _default_attr(node, "misalignment_psd", _misalignment_psd(0));
    _default_attr(node, "vel_psd", _vel_psd(0));
    _default_attr(node, "pos_psd", _pos_psd(0));
    _default_attr(node, "gyro_psd", _gyro_psd(0));
    _default_attr(node, "acce_psd", _acce_psd(0));
    _default_attr(node, "gyro_scale_psd", _gyro_scale_psd(0));
    _default_attr(node, "acce_scale_psd", _acce_scale_psd(0));
    _default_attr(node, "odo_psd", _odo_psd);

    _default_attr(node, "hor_vel_noise", _meas_vel_noise(0));
    _default_attr(node, "ver_vel_noise", _meas_vel_noise(2));
    _default_attr(node, "hor_pos_noise", _meas_pos_noise(0));
    _default_attr(node, "ver_pos_noise", _meas_pos_noise(2));

    str.clear();
    str = "";
    str = to_string(_lever(0)) + ' ' + to_string(_lever(1)) + ' ' + to_string(_lever(2));
    _default_attr(node, "lever", str);
    _default_attr(node, "odo_lever", str);
    _default_attr(node, "uwb_lever", str);

    _default_attr(node, "nq", _nq);
    _default_attr(node, "nr", _nr);
    _default_attr(node, "order", _order);
    _default_attr(node, "TS", _TS);
    _default_attr(node, "delay_t", _delay_t);
    _default_attr(node, "delay_odo", _delay_odo);
    _default_attr(node, "filtermode", _fltmode);
    _default_attr(node, "max_pdop", _max_pdop);
    _default_attr(node, "min_satnum", _min_sat);
    _default_attr(node, "odometry", _odo);
    _default_attr(node, "NHC", _NHC);
    _default_attr(node, "ZUPT", _ZUPT);
    _default_attr(node, "UWB", false);
    _default_attr(node, "Hgt", false);

    _gmutex.unlock();
}

void great::t_gsetign::help()
{
    _gmutex.lock();

    // no Implement

    _gmutex.unlock();
}


Eigen::Vector3d great::t_gsetign::initial_misalignment_std()
{

    _gmutex.lock();
    //string str = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).attribute("initial_misalignment_std").value();
    xml_node AttNode = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).child("Estimator").child("Attitude");
    string str = AttNode.attribute("InitialSTD").value();
    for (int i = 0; i < str.size(); i++)
    {
        if (str[i] == ',')
            str[i] = ' ';
    }
    stringstream ss(str);
    double X, Y, Z;
    ss >> X >> Y >> Z;

    _gmutex.unlock();
    if (_imu_type != IMU_TYPE::Customize)
        return _map_imu_error_models[_imu_type].AttInitialSTD;
    return Eigen::Vector3d(X, Y, Z);
}

Eigen::Vector3d great::t_gsetign::initial_vel_std()
{

    _gmutex.lock();
    //string str = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).attribute("initial_vel_std").value();    
    xml_node VelNode = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).child("Estimator").child("Velocity");
    string str = VelNode.attribute("InitialSTD").value();
    for (int i = 0; i < str.size(); i++)
    {
        if (str[i] == ',')
            str[i] = ' ';
    }
    stringstream ss(str);
    double X, Y, Z;
    ss >> X >> Y >> Z;
    _gmutex.unlock();
    if (_imu_type != IMU_TYPE::Customize)
        return _map_imu_error_models[_imu_type].VelInitialSTD;
    return Eigen::Vector3d(X, Y, Z);
}

Eigen::Vector3d great::t_gsetign::initial_pos_std()
{
    _gmutex.lock();
    //string str = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).attribute("initial_pos_std").value();
    xml_node PosNode = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).child("Estimator").child("Position");
    string str = PosNode.attribute("InitialSTD").value();
    for (int i = 0; i < str.size(); i++)
    {
        if (str[i] == ',')
            str[i] = ' ';
    }
    stringstream ss(str);
    double X, Y, Z;
    ss >> X >> Y >> Z;
    _gmutex.unlock();
    if (_imu_type != IMU_TYPE::Customize)
        return _map_imu_error_models[_imu_type].PosInitialSTD;
    return Eigen::Vector3d(X, Y, Z);
}

Eigen::Vector3d great::t_gsetign::initial_gyro_std()
{
    _gmutex.lock();
    //double x = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).attribute("initial_gyro_std").as_double();
    xml_node GyroNode = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).child("Estimator").child("GyroBias");
    string str = GyroNode.attribute("InitialSTD").value();
    for (int i = 0; i < str.size(); i++)
    {
        if (str[i] == ',')
            str[i] = ' ';
    }
    stringstream ss(str);
    double X, Y, Z;
    ss >> X >> Y >> Z;
    _gmutex.unlock();
    if (_imu_type != IMU_TYPE::Customize)
        return _map_imu_error_models[_imu_type].GyroBiasInitialSTD;
    return Eigen::Vector3d(X, Y, Z);
}

Eigen::Vector3d great::t_gsetign::initial_gyro_scale_std()
{
    _gmutex.lock();
    xml_node ExtraStatesNode = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).child("ExtraStates");
    string str = ExtraStatesNode.child("IMUScale").child("GyroScale").attribute("InitialSTD").value();
    for (int i = 0; i < str.size(); i++)
    {
        if (str[i] == ',')
            str[i] = ' ';
    }
    stringstream ss(str);
    double X, Y, Z;
    ss >> X >> Y >> Z;
    _gmutex.unlock();
    return Eigen::Vector3d(X, Y, Z);
}

Eigen::Vector3d great::t_gsetign::initial_acce_std()
{
    _gmutex.lock();
    //double x = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).attribute("initial_acce_std").as_double();
    xml_node AcceNode = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).child("Estimator").child("AcceBias");
    string str = AcceNode.attribute("InitialSTD").value();
    for (int i = 0; i < str.size(); i++)
    {
        if (str[i] == ',')
            str[i] = ' ';
    }
    stringstream ss(str);
    double X, Y, Z;
    ss >> X >> Y >> Z;
    _gmutex.unlock();
    if (_imu_type != IMU_TYPE::Customize)
        return _map_imu_error_models[_imu_type].AcceBiasInitialSTD;
    return Eigen::Vector3d(X, Y, Z);
}

Eigen::Vector3d great::t_gsetign::initial_acce_scale_std()
{
    _gmutex.lock();
    xml_node ExtraStatesNode = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).child("ExtraStates");
    string str = ExtraStatesNode.child("IMUScale").child("AcceScale").attribute("InitialSTD").value();
    for (int i = 0; i < str.size(); i++)
    {
        if (str[i] == ',')
            str[i] = ' ';
    }
    stringstream ss(str);
    double X, Y, Z;
    ss >> X >> Y >> Z;
    _gmutex.unlock();
    return Eigen::Vector3d(X, Y, Z);
}

Eigen::Vector3d great::t_gsetign::initial_imu_installation_att_std()
{
    _gmutex.lock();
    xml_node ExtraStatesNode = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).child("ExtraStates");
    string str = ExtraStatesNode.child("IMUInstallation").child("Rotation").attribute("InitialSTD").value();
    for (int i = 0; i < str.size(); i++)
    {
        if (str[i] == ',')
            str[i] = ' ';
    }
    stringstream ss(str);
    double X, Y, Z;
    ss >> X >> Y >> Z;
    _gmutex.unlock();
    return Eigen::Vector3d(X, Y, Z);
}

double great::t_gsetign::initial_odoscale_std()
{
    _gmutex.lock();
    xml_node ExtraStatesNode = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).child("ExtraStates");
    double res = ExtraStatesNode.child("OdometerScale").attribute("InitialSTD").as_double();
    _gmutex.unlock();
    return res;
}

Eigen::Vector3d great::t_gsetign::min_misalignment_std()
{
    _gmutex.lock();
    string str = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).attribute("min_misalignment_std").value();
    for (int i = 0; i < str.size(); i++)
    {
        if (str[i] == ',')
            str[i] = ' ';
    }
    stringstream ss(str);
    double X, Y, Z;
    ss >> X >> Y >> Z;
    _gmutex.unlock();
    return Eigen::Vector3d(X, Y, Z);
}

Eigen::Vector3d great::t_gsetign::min_vel_std()
{
    _gmutex.lock();
    string str = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).attribute("min_vel_std").value();
    for (int i = 0; i < str.size(); i++)
    {
        if (str[i] == ',')
            str[i] = ' ';
    }
    stringstream ss(str);
    double X, Y, Z;
    ss >> X >> Y >> Z;
    _gmutex.unlock();
    return Eigen::Vector3d(X, Y, Z);
}

Eigen::Vector3d great::t_gsetign::min_pos_std()
{
    _gmutex.lock();
    string str = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).attribute("min_pos_std").value();
    for (int i = 0; i < str.size(); i++)
    {
        if (str[i] == ',')
            str[i] = ' ';
    }
    stringstream ss(str);
    double X, Y, Z;
    ss >> X >> Y >> Z;
    _gmutex.unlock();
    return Eigen::Vector3d(X, Y, Z);
}

Eigen::Vector3d great::t_gsetign::min_gyro_std()
{
    _gmutex.lock();
    string str = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).attribute("min_gyro_std").value();
    for (int i = 0; i < str.size(); i++)
    {
        if (str[i] == ',')
            str[i] = ' ';
    }
    stringstream ss(str);
    double X, Y, Z;
    ss >> X >> Y >> Z;
    _gmutex.unlock();
    return Eigen::Vector3d(X, Y, Z);
}

Eigen::Vector3d great::t_gsetign::min_acce_std()
{
    _gmutex.lock();
    string str = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).attribute("min_acce_std").value();
    for (int i = 0; i < str.size(); i++)
    {
        if (str[i] == ',')
            str[i] = ' ';
    }
    stringstream ss(str);
    double X, Y, Z;
    ss >> X >> Y >> Z;
    _gmutex.unlock();
    return Eigen::Vector3d(X, Y, Z);
}

double great::t_gsetign::min_odo_std()
{
    _gmutex.lock();
    double x = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).attribute("min_odo_k_std").as_double();
    _gmutex.unlock();
    return x;
}

Eigen::Vector3d great::t_gsetign::misalignment_psd()
{
    _gmutex.lock();
    //double x = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).attribute("misalignment_psd").as_double();
    xml_node AttNode = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).child("Estimator").child("Attitude");
    string str = AttNode.attribute("ProcNoiseSD").value();
    for (int i = 0; i < str.size(); i++)
    {
        if (str[i] == ',')
            str[i] = ' ';
    }
    stringstream ss(str);
    double X, Y, Z;
    ss >> X >> Y >> Z;

    _gmutex.unlock();
    if (_imu_type != IMU_TYPE::Customize)
        return _map_imu_error_models[_imu_type].AttProcNoisePSD;
    return Eigen::Vector3d(X, Y, Z);
}

Eigen::Vector3d great::t_gsetign::vel_psd()
{
    _gmutex.lock();
    //double x = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).attribute("vel_psd").as_double();
    xml_node VelNode = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).child("Estimator").child("Velocity");
    string str = VelNode.attribute("ProcNoiseSD").value();
    for (int i = 0; i < str.size(); i++)
    {
        if (str[i] == ',')
            str[i] = ' ';
    }
    stringstream ss(str);
    double X, Y, Z;
    ss >> X >> Y >> Z;
    _gmutex.unlock();
    if (_imu_type != IMU_TYPE::Customize)
        return _map_imu_error_models[_imu_type].VelProcNoisePSD;
    return Eigen::Vector3d(X, Y, Z);
}

Eigen::Vector3d great::t_gsetign::pos_psd()
{
    _gmutex.lock();
    //double x = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).attribute("pos_psd").as_double();
    xml_node PosNode = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).child("Estimator").child("Position");
    string str = PosNode.attribute("ProcNoiseSD").value();
    for (int i = 0; i < str.size(); i++)
    {
        if (str[i] == ',')
            str[i] = ' ';
    }
    stringstream ss(str);
    double X, Y, Z;
    ss >> X >> Y >> Z;
    _gmutex.unlock();
    if (_imu_type != IMU_TYPE::Customize)
        return _map_imu_error_models[_imu_type].PosProcNoisePSD;
    return Eigen::Vector3d(X, Y, Z);
}

Eigen::Vector3d great::t_gsetign::gyro_psd()
{
    _gmutex.lock();
    //double x = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).attribute("gyro_psd").as_double();
    xml_node GyroNode = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).child("Estimator").child("GyroBias");
    string str = GyroNode.attribute("ProcNoiseSD").value();
    for (int i = 0; i < str.size(); i++)
    {
        if (str[i] == ',')
            str[i] = ' ';
    }
    stringstream ss(str);
    double X, Y, Z;
    ss >> X >> Y >> Z;
    _gmutex.unlock();
    if (_imu_type != IMU_TYPE::Customize)
        return _map_imu_error_models[_imu_type].GyroBiasProcNoisePSD;
    return Eigen::Vector3d(X, Y, Z);
}

Eigen::Vector3d great::t_gsetign::acce_psd()
{
    _gmutex.lock();
    //double x = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).attribute("acce_psd").as_double();
    xml_node AcceNode = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).child("Estimator").child("AcceBias");
    string str = AcceNode.attribute("ProcNoiseSD").value();
    for (int i = 0; i < str.size(); i++)
    {
        if (str[i] == ',')
            str[i] = ' ';
    }
    stringstream ss(str);
    double X, Y, Z;
    ss >> X >> Y >> Z;
    _gmutex.unlock();
    if (_imu_type != IMU_TYPE::Customize)
        return _map_imu_error_models[_imu_type].AcceBiasProcNoisePSD;
    return Eigen::Vector3d(X, Y, Z);
}

Eigen::Vector3d great::t_gsetign::gyro_scale_psd()
{
    _gmutex.lock();
    xml_node ExtraStatesNode = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).child("ExtraStates");
    string str = ExtraStatesNode.child("IMUScale").child("GyroScale").attribute("ProcNoiseSD").value();
    for (int i = 0; i < str.size(); i++)
    {
        if (str[i] == ',')
            str[i] = ' ';
    }
    stringstream ss(str);
    double X, Y, Z;
    ss >> X >> Y >> Z;
    _gmutex.unlock();
    return Eigen::Vector3d(X, Y, Z);
}

Eigen::Vector3d great::t_gsetign::acce_scale_psd()
{
    _gmutex.lock();
    xml_node ExtraStatesNode = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).child("ExtraStates");
    string str = ExtraStatesNode.child("IMUScale").child("AcceScale").attribute("ProcNoiseSD").value();
    for (int i = 0; i < str.size(); i++)
    {
        if (str[i] == ',')
            str[i] = ' ';
    }
    stringstream ss(str);
    double X, Y, Z;
    ss >> X >> Y >> Z;
    _gmutex.unlock();
    return Eigen::Vector3d(X, Y, Z);
}

Eigen::Vector3d great::t_gsetign::imu_inst_trans_psd()
{
    _gmutex.lock();
    xml_node ExtraStatesNode = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).child("ExtraStates");
    string str = ExtraStatesNode.child("IMUInstallation").child("Translation").attribute("ProcNoiseSD").value();
    for (int i = 0; i < str.size(); i++)
    {
        if (str[i] == ',')
            str[i] = ' ';
    }
    stringstream ss(str);
    double X, Y, Z;
    ss >> X >> Y >> Z;
    _gmutex.unlock();
    return Eigen::Vector3d(X, Y, Z);
}

Eigen::Vector3d great::t_gsetign::imu_inst_rot_psd()
{
    _gmutex.lock();
    xml_node ExtraStatesNode = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).child("ExtraStates");
    string str = ExtraStatesNode.child("IMUInstallation").child("Rotation").attribute("ProcNoiseSD").value();
    for (int i = 0; i < str.size(); i++)
    {
        if (str[i] == ',')
            str[i] = ' ';
    }
    stringstream ss(str);
    double X, Y, Z;
    ss >> X >> Y >> Z;
    _gmutex.unlock();
    return Eigen::Vector3d(X, Y, Z);
}

double great::t_gsetign::odo_scale()
{
    _gmutex.lock();
    xml_node OdometerNode = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).child("Odometer");
    double x = OdometerNode.child("Scale").attribute("Value").as_double();
    _gmutex.unlock();
    return x;
}

double great::t_gsetign::odo_psd()
{
    _gmutex.lock();
    xml_node ExtraStatesNode = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).child("ExtraStates");
    double x = ExtraStatesNode.child("OdometerScale").attribute("ProcNoiseSD").as_double();
    _gmutex.unlock();
    return x;
}

double great::t_gsetign::odo_std()
{
    _gmutex.lock();
    xml_node OdometerNode = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).child("Odometer");
    double x = OdometerNode.child("NoiseSD").attribute("Value").as_double();
    _gmutex.unlock();
    return x;
}

Eigen::Vector3d great::t_gsetign::NHC_std()
{
    _gmutex.lock();
    xml_node NHCNode = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).child("NHC");
    double x = NHCNode.child("NoiseSD").attribute("Value").as_double();
    _gmutex.unlock();
    return Eigen::Vector3d(x, x, x);
}

Eigen::Vector3d great::t_gsetign::ZUPT_std()
{
    _gmutex.lock();
    xml_node ZUPTNode = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).child("ZUPT");
    double x = ZUPTNode.child("VelNoiseSD").attribute("Value").as_double();
    _gmutex.unlock();
    return Eigen::Vector3d(x, x, x);
}
double great::t_gsetign::ZIHR_std()
{
    _gmutex.lock();
    xml_node ZUPTNode = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).child("ZUPT");
    double x = ZUPTNode.child("ZIHRNoiseSD").attribute("Value").as_double();
    _gmutex.unlock();
    return x * t_gglv::deg;
}

double great::t_gsetign::Yaw_std()
{
    _gmutex.lock();
    xml_node ZUPTNode = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).child("ZUPT");
    double x = ZUPTNode.child("YawNoiseSD").attribute("Value").as_double();
    _gmutex.unlock();
    return x * t_gglv::deg;
}

Eigen::Vector3d great::t_gsetign::ATT_std()
{
    _gmutex.lock();
    xml_node AttNode = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).child("Attitude");
    string x = AttNode.child("NoiseSD").attribute("Value").value();
    for (int i = 0; i < x.size(); i++)
    {
        if (x[i] == ',')x[i] = ' ';
    }
    stringstream ss(x);
    double p, r, y;
    ss >> p >> r >> y;
    _gmutex.unlock();
    return Eigen::Vector3d(p, r, y) * t_gglv::deg;
}

int great::t_gsetign::nq()
{
    _gmutex.lock();
    double x = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).attribute("nq").as_double();
    _gmutex.unlock();
    return x;
}

int great::t_gsetign::nr()
{
    _gmutex.lock();
    double x = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).attribute("nr").as_double();
    _gmutex.unlock();
    return x;
}

double great::t_gsetign::delay_odo()
{
    _gmutex.lock();
    double x = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).child("Odometer").child("DelayTime").attribute("Value").as_double();
    _gmutex.unlock();
    return x;
}

Eigen::Vector3d great::t_gsetign::odo_lever()
{
    _gmutex.lock();
    // IMU lever from imu center to wheel
    xml_node lever = _doc.child(XMLKEY_ROOT).child(XMLKEY_INS).child("Installation").child("Lever");
    double lx = lever.attribute("X").as_double();
    double ly = lever.attribute("Y").as_double();
    double lz = lever.attribute("Z").as_double();
    _gmutex.unlock();
    return Vector3d(lx, ly, lz);
}

Eigen::Vector3d great::t_gsetign::uwb_lever()
{
    _gmutex.lock();
    xml_node lever = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).child("UltraWideBand").child("AntennaLever");
    double lx = lever.attribute("X").as_double();
    double ly = lever.attribute("Y").as_double();
    double lz = lever.attribute("Z").as_double();
    _gmutex.unlock();
    return Vector3d(lx, ly, lz);
}

IGN_TYPE great::t_gsetign::ign_type()
{
    _gmutex.lock();
    string x = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).child("GNSS").attribute("Type").value();
    _gmutex.unlock();
    return great::str2ign(x);
}

map<double, int> great::t_gsetign::sim_gnss_outages()
{
    map<double, int> outages;
    xml_node SimOutages = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).child("GNSS").child("SimOutages");
    string beg = SimOutages.attribute("beg").value();
    string len = SimOutages.attribute("length").value();
    stringstream begss(beg), lenss(len);
    double begt;
    int length;
    while (begss >> begt)
    {
        lenss >> length;
        outages.insert(make_pair(begt, length));
    }
    return outages;
}

IMU_TYPE great::t_gsetign::imu_type()
{
    _gmutex.lock();
    string x = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).child("GNSS").attribute("Type").value();
    _gmutex.unlock();
    return great::str2imu(x);
}

Eigen::Vector3d great::t_gsetign::lever()
{
    _gmutex.lock();
    xml_node LeverNode = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).child("GNSS").child("AntennaLever");
    string type = LeverNode.attribute("Type").value();
    double X = LeverNode.attribute("X").as_double();
    double Y = LeverNode.attribute("Y").as_double();
    double Z = LeverNode.attribute("Z").as_double();
    if (type == "RFU")
    {
        // X; Y; Z;
    }

    _gmutex.unlock();
    return Vector3d(X, Y, Z);
}

double great::t_gsetign::delay_t()
{
    _gmutex.lock();
    double x = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).child("GNSS").child("DelayTime").attribute("Value").as_double();
    _gmutex.unlock();
    return x;
}

double great::t_gsetign::max_pdop()
{
    _gmutex.lock();
    double x = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).child("GNSS").child("LCISetting").attribute("MaxPDOP").as_double();
    _gmutex.unlock();
    return x;
}

int great::t_gsetign::min_sat()
{
    _gmutex.lock();
    //int x = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).attribute("min_satnum").as_int();
    double x = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).child("GNSS").child("LCISetting").attribute("MinSat").as_double();
    _gmutex.unlock();
    return x;
}

string great::t_gsetign::odo()
{
    _gmutex.lock();
    string str = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).child("Odometer").attribute("Type").value();
    _gmutex.unlock();
    return str;
}

string great::t_gsetign::odo_inst()
{
    _gmutex.lock();
    xml_node OdometerNode = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).child("Odometer");
    string str = OdometerNode.child("Installation").attribute("Type").value();
    _gmutex.unlock();
    return str;
}

bool great::t_gsetign::NHC()
{
    _gmutex.lock();
    string str = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).child("NHC").attribute("Type").value();
    bool b = false;
    if (str == "ON")b = true;
    _gmutex.unlock();
    return b;
}

bool great::t_gsetign::ZUPT()
{
    _gmutex.lock();
    string str = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).child("ZUPT").attribute("Type").value();
    bool b = false;
    if (str == "ON")b = true;
    _gmutex.unlock();
    return b;
}

bool great::t_gsetign::Attitude()
{
    _gmutex.lock();
    string str = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).child("Attitude").attribute("Type").value();
    bool b = false;
    if (str == "ON")b = true;
    _gmutex.unlock();
    return b;
}

bool great::t_gsetign::UWB()
{
    _gmutex.lock();
    string str = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).child("UltraWideBand").attribute("Type").value();
    bool b = false;
    if (str == "ON")b = true;
    _gmutex.unlock();
    return b;
}

bool great::t_gsetign::Hgt()
{
    _gmutex.lock();
    string str = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).child("Hgt").attribute("Type").value();
    bool b = false;
    if (str == "ON")b = true;
    _gmutex.unlock();
    return b;
}

bool great::t_gsetign::imu_scale()
{
    _gmutex.lock();
    xml_node ExtraStatesNode = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).child("ExtraStates");
    string str = ExtraStatesNode.child("IMUScale").attribute("Type").value();
    bool b = false;
    if (str == "ON")b = true;
    _gmutex.unlock();
    return b;
}

bool great::t_gsetign::imu_inst_rot()
{
    _gmutex.lock();
    xml_node ExtraStatesNode = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).child("ExtraStates");
    string str = ExtraStatesNode.child("IMUInstallation").attribute("Type").value();
    bool b = false;
    if (str == "Rotation" || str == "TransRot")b = true;
    _gmutex.unlock();
    return b;
}

bool great::t_gsetign::imu_inst_trans()
{
    _gmutex.lock();
    xml_node ExtraStatesNode = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).child("ExtraStates");
    string str = ExtraStatesNode.child("IMUInstallation").attribute("Type").value();
    bool b = false;
    if (str == "Translation" || str == "TransRot")b = true;
    _gmutex.unlock();
    return b;
}

double great::t_gsetign::wheelraduis()
{
    _gmutex.lock();
    xml_node OdometerNode = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).child("Odometer");
    double d = OdometerNode.child("WheelRadius").attribute("Value").as_double();
    _gmutex.unlock();
    return d;
}

double great::t_gsetign::UWB_start()
{
    _gmutex.lock();
    double res = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).child("UltraWideBand").child("Time").attribute("Start").as_double();

    _gmutex.unlock();
    return res;
}

double great::t_gsetign::UWB_end()
{
    _gmutex.lock();
    double res = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).child("UltraWideBand").child("Time").attribute("End").as_double();

    _gmutex.unlock();
    return res;
}

double great::t_gsetign::Hgt_start()
{
    _gmutex.lock();
    double res = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).child("Hgt").child("Time").attribute("Start").as_double();

    _gmutex.unlock();
    return res;
}

double great::t_gsetign::Hgt_end()
{
    _gmutex.lock();
    double res = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).child("Hgt").child("Time").attribute("End").as_double();

    _gmutex.unlock();
    return res;
}

double great::t_gsetign::Hgt_std()
{
    _gmutex.lock();
    double res = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).child("Hgt").child("STD").attribute("Value").as_double();

    _gmutex.unlock();
    return res;
}

double great::t_gsetign::Hgt_info()
{
    _gmutex.lock();
    double res = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).child("Hgt").attribute("Value").as_double();

    _gmutex.unlock();
    return res;
}

map<string, Eigen::Vector3d> great::t_gsetign::MultiAntLever()
{
    _gmutex.lock();
    int i = 1;
    map<string, Eigen::Vector3d> levers;
    set<string> sites = t_gsetbase::_setvals("gen", "rover");
    set<string> rover_sites = t_gsetbase::_setvals("gen", "rover");
    if (rover_sites.size() > 0)
    {
        sites = set<string>(rover_sites.begin(), rover_sites.end());
    }
    for (auto it : sites)
    {
        string Ant = "AntennaLever" + to_string(i);
        xml_node LeverNode = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).child("GNSS").child(Ant.c_str());
        string name = LeverNode.attribute("Name").value();
        if (name != it)
        {
            _gmutex.unlock();
            return map<string, Eigen::Vector3d>();
        }
        string type = LeverNode.attribute("Type").value();
        if (type == "RFU")
        {
            // X; Y; Z;
        }
        double X = LeverNode.attribute("X").as_double();
        double Y = LeverNode.attribute("Y").as_double();
        double Z = LeverNode.attribute("Z").as_double();
        levers.insert(make_pair(name, Eigen::Vector3d(X, Y, Z)));
        ++i;
    }
    _gmutex.unlock();
    return levers;
}

map<MEAS_TYPE, double> great::t_gsetign::max_norm()
{
    _gmutex.lock();
    map<MEAS_TYPE, double> map_norm;
    xml_node ign = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN);
    xml_node::iterator itnode = ign.begin();
    while (itnode != ign.end())
    {
        string name = itnode->name();
        MEAS_TYPE meas = xmlname2meas(name);
        if (meas != MEAS_TYPE::DEFAULT_MEAS)
        {
            double norm = 3.0;
            xml_node node = itnode->child("MaxNorm");
            if (!node.empty())
                norm = node.attribute("Value").as_double();
            if (meas == MEAS_TYPE::GNSS_MEAS)
            {
                norm = itnode->child("LCISetting").attribute("MaxNorm").as_double();
                map_norm.insert(make_pair(MEAS_TYPE::POS_MEAS, norm));
                map_norm.insert(make_pair(MEAS_TYPE::POS_VEL_MEAS, norm));

                ++itnode;
                continue;
            }
            map_norm.insert(make_pair(meas, norm));

        }

        ++itnode;
    }
    _gmutex.unlock();
    return map_norm;
}

MEAS_TYPE great::t_gsetign::xmlname2meas(string str)
{
    if (str == "GNSS")
        return MEAS_TYPE::GNSS_MEAS;
    else if (str == "Odometer")
        return MEAS_TYPE::ODO_MEAS;
    else if (str == "NHC")
        return MEAS_TYPE::NHC_MEAS;
    else if (str == "ZUPT")
        return MEAS_TYPE::ZUPT_MEAS;
    else if (str == "Attitude")
        return MEAS_TYPE::ATT_MEAS;
    else
        return MEAS_TYPE::DEFAULT_MEAS;
}
