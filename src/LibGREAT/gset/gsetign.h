/**
 * @file         gsetign.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        control set from XML for MSF
 * @version      1.0
 * @date         2025-01-01
 *
 * @copyright Copyright (c) 2025, Wuhan University. All rights reserved.
 *
 */

#ifndef GSETIGN_H
#define GSETIGN_H

#define XMLKEY_IGN "integration"

#include "gexport/ExportLibGREAT.h"
#include "gexport/ExportLibGREAT.h"
#include "gset/gsetbase.h"
#include "gset/gsetins.h"
#include "gutils/gtime.h"
using namespace gnut;

namespace great
{
    /**
    * @enum MEAS_TYPE
    * @brief Measurement mode,included GNSS,ZUPT,Motion Constrain,Yaw et al
    */
    enum MEAS_TYPE
    {
        GNSS_MEAS = 1,        ///< GNSS measurement
        POS_MEAS,             ///< position measurement
        VEL_MEAS,             ///< velocity measurement
        POS_VEL_MEAS,         ///< TODO

        MOTION_MEAS = 20, ///< motion measurement
        ZUPT_MEAS,        ///< zero velocity update measurement
        ZIHR_MEAS,        ///< zero integrated heading rate update measurement
        ODO_MEAS,         ///< odometer measurement
        NHC_MEAS,         ///< nonholonomic constraints measurement
        YAW_MEAS,         ///< yaw measurement
        Hgt_MEAS,         ///< height measurement
        ATT_MEAS,         ///< attitude measurement
        ZUPT_POS_MEAS,    ///< zero velocity update position measurement

        OTHER_MEAS = 34, ///< other measurement
        VIS_MEAS,        ///< vision measurement
        LIDAR_MEAS,      ///< LiDAR measurement

        Indoor_MEAS = 50, ///< indoor measurement
        UWB_MEAS,         ///< UWB measurement
        WIFI_MEAS,        ///< WiFi measurement
        BLE_MEAS,         ///< BLE measurement
        IMGTAG_MEAS,      ///< image tag measurement

        DEFAULT_MEAS = 100, ///< default measurement
        NO_MEAS             ///< none
    };

    /**
    * @enum IGN_TYPE
    * @brief GNSS/INS integration option.
    * PURE_INS is ins process only;
    * LCI is loosely coupled integration;
    * TCI is tightly coupled intergration;
    * STCI is semi-tightly coupled integration.
    * MULTIGN_LC is GNSS/INS/VISION loosely coupled integrartion
    * MULTIGN_TC is GNSS/INS/VISION tightly coupled integrartion
    * MULTIGN_STC is GNSS/INS/VISION semi-tightly coupled integrartion
    * MULTIGN_LIDAR_LC is GNSS/INS/LIDAR loosely coupled integrartion
    * MULTIGN_LIDAR_TC is GNSS/INS/LIDAR tightly coupled integrartion
    * MULTIGN_LIDAR_STC is GNSS/INS/LIDAR semi-tightly coupled integrartion
    * LIDAR_INS is INS/LIDAR tightly coupled integration
    * VIO is INS/VISION tightly coupled integration
    */
    enum IGN_TYPE
    {
        PURE_INS,     ///< ins process only
        LCI,          ///< loosely coupled integration
        TCI,          ///< tightly coupled intergration
        STCI,         ///< semi-tightly coupled integration
        VIS_LCI,      ///< INS/VISION coupled integration
        VIS_TCI,      ///< INS/VISION tightly integration
        VIS_STCI,     ///< INS/VISION semi-tightly integration
        VIO,          ///< INS/VISION coupled integration
        LIDAR_LCI,    ///< INS/LIDAR loosely coupled integration
        LIDAR_TCI,    ///< INS/LIDAR tightly coupled integration
        LIDAR_STCI,   ///< INS/LIDAR semi-tightly coupled integration
        LIO,          ///< INS/LIDAR coupled integration
        MULTIGN_LCI,  ///< GNSS/INS/VISION loosely coupled integrartion
        MULTIGN_TCI,  ///< GNSS/INS/VISION tightly coupled integrartion
        MULTIGN_STCI, ///< GNSS/INS/VISION semi-tightly coupled integrartion
        VLO,          ///< INS/VISION/LIDAR coupled integration
        UWB_LCI,      ///< INS/UWB loosely coupled integration
        UWB_TCI,      ///< INS/UWB tightly coupled integration
        IGN_DEFAULT   ///< defalut coupled integration
    };

    

    
    /**
    * @brief change MEAS_TYPE to string.
    * @param[in]    type        type of measurement
    * @return        string        type of measurement in string form
    */
    LibGREAT_LIBRARY_EXPORT string meas2str(MEAS_TYPE type);

    /**
    * @brief change string to IGN_TYPE.
    * @param[in]    s            type of integration
    * @return        IGN_TYPE    type of integration in string form
    */
    LibGREAT_LIBRARY_EXPORT IGN_TYPE str2ign(const string& s);
    LibGREAT_LIBRARY_EXPORT IMU_TYPE str2imu(const string& s);
    /**
    * @class t_gsetign
    * @brief coupled integration settings
    * @author ShenZhiheng
    */
    class LibGREAT_LIBRARY_EXPORT t_gsetign : public virtual t_gsetbase
    {
    public:
        /** @brief default constructor. */
        t_gsetign();

        /** @brief default destructor. */
        ~t_gsetign();

        /**
        * @brief settings check.
        */
        void check();

        /**
        * @brief settings help.
        */
        void help();

        /**
        * @brief    get initial standard convariance of misalignment.
        * @return    Eigen::Vector3d     initial standard convariance of misalignment
        */
        Eigen::Vector3d initial_misalignment_std();

        /**
        * @brief    get initial standard convariance of velocity.
        * @return    Eigen::Vector3d     initial standard convariance of velocity
        */
        Eigen::Vector3d initial_vel_std();

        /**
        * @brief    get initial standard convariance of position.
        * @return    Eigen::Vector3d     initial standard convariance of position
        */
        Eigen::Vector3d initial_pos_std();

        /**
        * @brief    get initial standard convariance of gyro.
        * @return    Eigen::Vector3d     initial standard convariance of gyro
        */
        Eigen::Vector3d initial_gyro_std();

        /**
         * @brief
         *
         * @return Eigen::Vector3d
         */
        Eigen::Vector3d initial_gyro_scale_std();

        /**
        * @brief    get initial standard convariance of accelerometer.
        * @return    Eigen::Vector3d     initial standard convariance of accelerometer
        */
        Eigen::Vector3d initial_acce_std();

        /**
         * @brief
         *
         * @return Eigen::Vector3d
         */
        Eigen::Vector3d initial_acce_scale_std();

        /**
         * @brief
         *
         * @return Eigen::Vector3d
         */
        Eigen::Vector3d initial_imu_installation_att_std();

        /**
         * @brief
         *
         * @return double
         */
        double initial_odoscale_std(); ///< TODO

        /**
        * @brief    get minimum standard convariance of misalignment.
        * @return    Eigen::Vector3d     minimum standard convariance of misalignment
        */
        Eigen::Vector3d min_misalignment_std();

        /**
        * @brief    get minimum standard convariance of velocity.
        * @return    Eigen::Vector3d     minimum standard convariance of velocity
        */
        Eigen::Vector3d min_vel_std();

        /**
        * @brief    get minimum standard convariance of position.
        * @return    Eigen::Vector3d     minimum standard convariance of position
        */
        Eigen::Vector3d min_pos_std();

        /**
        * @brief    get minimum standard convariance of gyro.
        * @return    Eigen::Vector3d     minimum standard convariance of gyro
        */
        Eigen::Vector3d min_gyro_std();

        /**
        * @brief    get minimum standard convariance of accelerometer.
        * @return    Eigen::Vector3d     minimum standard convariance of accelerometer
        */
        Eigen::Vector3d min_acce_std();

        /**
        * @brief    get minimum standard convariance of odometer.
        * @return    double     minimum standard convariance of odometer
        */
        double min_odo_std();

        /**
        * @brief    get noise power spectral density of misalignment.
        * @return    Eigen::Vector3d     noise power spectral density of misalignment
        */
        Eigen::Vector3d misalignment_psd();

        /**
        * @brief    get noise power spectral density of velocity.
        * @return    Eigen::Vector3d     noise power spectral density of velocity
        */
        Eigen::Vector3d vel_psd();

        /**
        * @brief    get noise power spectral density of position.
        * @return    Eigen::Vector3d     noise power spectral density of position
        */
        Eigen::Vector3d pos_psd();

        /**
        * @brief    get noise power spectral density of gyro.
        * @return    Eigen::Vector3d     noise power spectral density of gyro
        */
        Eigen::Vector3d gyro_psd();

        /**
        * @brief    get noise power spectral density of accelerometer.
        * @return    Eigen::Vector3d     noise power spectral density of accelerometer
        */
        Eigen::Vector3d acce_psd();

        /**
         * @brief
         *
         * @return Eigen::Vector3d
         */
        Eigen::Vector3d gyro_scale_psd();

        /**
         * @brief
         *
         * @return Eigen::Vector3d
         */
        Eigen::Vector3d acce_scale_psd();

        /**
         * @brief
         *
         * @return Eigen::Vector3d
         */
        Eigen::Vector3d imu_inst_rot_psd(); // imu installation attitude psd

        /**
         * @brief
         *
         * @return Eigen::Vector3d
         */
        Eigen::Vector3d imu_inst_trans_psd(); // imu translation attitude psd

        /**
        * @brief    get noise power spectral density of odometer.
        * @return    double     noise power spectral density of odometer
        */
        double odo_scale();

        /**
         * @brief
         *
         * @return double
         */
        double odo_psd();

        /**
         * @brief
         *
         * @return double
         */
        double odo_std();

        /**
         * @brief
         *
         * @return Eigen::Vector3d
         */
        Eigen::Vector3d NHC_std();

        /**
         * @brief
         *
         * @return Eigen::Vector3d
         */
        Eigen::Vector3d ZUPT_std();

        /**
         * @brief
         *
         * @return double
         */
        double ZIHR_std();

        /**
         * @brief
         *
         * @return double
         */
        double Yaw_std();

        /**
         * @brief
         *
         * @return Eigen::Vector3d
         */
        Eigen::Vector3d ATT_std();

        /**
        * @brief    get number of ins state vector.
        * @return    int     number of ins state vector
        */
        int nq();

        /**
         * @brief
         * @return int
         */
        int nr(); ///< TODO

        /**
        * @brief    get maximum time delay.
        * @return    double     maximum time delay
        */
        double delay_t();

        /**
         * @brief
         *
         * @return double
         */
        double delay_odo(); ///< TODO

        /**
        * @brief    get maximum of PDOP.
        * @return    double     maximum of PDOP
        */
        double max_pdop();

        /**
        * @brief    get minimum number of satellite.
        * @return    int     minimum number of satellite
        */
        int min_sat();

        /**
        * @brief    get lever arm.
        * @return    Eigen::Vector3d     lever arm
        */
        Eigen::Vector3d lever();

        /**
         * @brief
         *
         * @return Eigen::Vector3d
         */
        Eigen::Vector3d odo_lever();

        /**
         * @brief
         *
         * @return Eigen::Vector3d
         */
        Eigen::Vector3d uwb_lever();

        /**
        * @brief    get integration type.
        * @return    IGN_TYPE     integration type
        */
        IGN_TYPE ign_type();


        map<double, int> sim_gnss_outages();

        /**
         * @brief
         *
         * @return IMU_TYPE
         */
        IMU_TYPE imu_type();


        /**
         * @brief get the format of ODO file
         *
         * @return string
         * @return "OFF/Raw/Pulse/Velocity"
         */
        string odo(); 

        /**
         * @brief get the installation of ODO
         *
         * @return string
         * @return "Left/Right"
         */
        string odo_inst(); ///< TODO

        /**
         * @brief
         *
         * @return true
         * @return false
         */
        bool NHC();

        /**
         * @brief
         *
         * @return true
         * @return false
         */
        bool ZUPT();

        /**
         * @brief
         *
         * @return true
         * @return false
         */
        bool Attitude();

        /**
         * @brief
         * @return true
         * @return false
         */
        bool UWB();

        /**
         * @brief
         * @return true
         * @return false
         */
        bool Hgt();

        /**
         * @brief
         * @return true
         * @return false
         */
        bool imu_scale();

        /**
         * @brief
         * @return true
         * @return false
         */
        bool imu_inst_rot();

        /**
         * @brief
         * @return true
         * @return false
         */
        bool imu_inst_trans();

        /**
         * @brief
         * @return double
         */
        double wheelraduis();

        double UWB_start();

        double UWB_end();

        double Hgt_start();

        double Hgt_end();

        double Hgt_std();

        double Hgt_info();

        map<string, Eigen::Vector3d> MultiAntLever();

        map<MEAS_TYPE, double> max_norm();
        MEAS_TYPE xmlname2meas(string str);

    protected:
        Eigen::Vector3d _initial_misalignment_std; ///< Standard deviation of initial misalignment angle [deg].
        Eigen::Vector3d _initial_vel_std;          ///< Standard deviation of initial velocity [m/s].
        Eigen::Vector3d _initial_pos_std;          ///< Standard deviation of initial position [m].
        Eigen::Vector3d _initial_gyro_std;         ///< Standard deviation of initial gyro drift [deg/h].
        Eigen::Vector3d _initial_acce_std;         ///< Standard deviation of initial accelerator bias [mg].
        Eigen::Vector3d _initial_gyro_scale_std;   ///< Standard deviation of initial gyro scale [#].
        Eigen::Vector3d _initial_acce_scale_std;   ///< Standard deviation of initial accelerator scale [#].
        double _initial_odoscale_std;              ///< TODO

        Eigen::Vector3d _min_misalignment_std; ///< Standard deviation lower bound of misalignment angle [deg].
        Eigen::Vector3d _min_vel_std;          ///< Standard deviation lower bound of velocity [m/s].
        Eigen::Vector3d _min_pos_std;          ///< Standard deviation lower bound of position [m].
        Eigen::Vector3d _min_gyro_std;         ///< Standard deviation lower bound of gyro drift [deg/h].
        Eigen::Vector3d _min_acce_std;         ///< Standard deviation lower bound of accelerator bias [mg].
        double _odo_std;

        Eigen::Vector3d _misalignment_psd; ///< Spectral density of misalignment angle [deg/sqrt(h)].
        Eigen::Vector3d _vel_psd;          ///< Spectral density of velocity [mg/sqrt(Hz)].
        Eigen::Vector3d _pos_psd;          ///< Spectral density of position
        Eigen::Vector3d _gyro_psd;         ///< Spectral density of gyro bias
        Eigen::Vector3d _acce_psd;         ///< Spectral density of accelerometer bias
        Eigen::Vector3d _gyro_scale_psd;   ///< Spectral density of gyro scale
        Eigen::Vector3d _acce_scale_psd;   ///< Spectral density of accelerometer scale
        double _odo_psd;

        Eigen::Vector3d _meas_vel_noise; ///< Measurement noise of velocity [m/s]
        Eigen::Vector3d _meas_pos_noise; ///< Measurement noise of position [m]

        Eigen::Vector3d _lever, _odo_lever, _uwb_lever; ///< GNSS Receiver center's lever relative to imu center on b Coor
        int _nq, _nr;                                   ///< parameter dimension and measurement dimension
        FLT_TYPE _fltmode;                       ///< filter mode (forward, forward and backward smooth, RTS)
        double _TS;                                     ///< measurement data interval
        string _order;                                  ///< the order for pos file (size is 6)
        double _delay_t;                                ///< the GNSS misalignment time diff with imu
        double _delay_odo;                              ///< the ODO misalignment time diff with imu
        double _max_pdop;                               ///< the max PDOP threshold for LCI
        int _min_sat;                                   ///< the minimum sat num threshold.
        IGN_TYPE _ign_type;                             ///< the integrated mode
        IMU_TYPE _imu_type;                             ///< IMU Error Models
        map<IMU_TYPE, ErrorModel> _map_imu_error_models;
        string _odo;                                    ///< odometry
        bool _NHC, _ZUPT;                                ///< motion constraints
    };
}

#endif
