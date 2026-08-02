/**
 * @file         gsetins.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        control set from XML for INS
 * @version      1.0
 * @date         2025-01-01
 *
 * @copyright Copyright (c) 2025, Wuhan University. All rights reserved.
 *
 */

#ifndef GSETINS_H
#define GSETINS_H

#define XMLKEY_INS "ins"

#include "gexport/ExportLibGREAT.h"
#include "gset/gsetbase.h"
#include "gins/gbase.h"
#include <Eigen/Eigen>

using namespace gnut;

namespace great
{
    /**
    * @class t_gsetins
    * @brief ins configure
    * @author ShenZhiheng
    */
    //class LibSet_LIBRARY_EXPORT t_gglv;

   



    class LibGREAT_LIBRARY_EXPORT t_gsetins : public virtual t_gsetbase
    {
    public:
        /** @brief default constructor. */
        t_gsetins();

        /** @brief default destructor. */
        virtual ~t_gsetins();

        /**
        * @brief settings check.
        */
        void check();


       
        /**
        * @brief settings help.
        */
        void help();

        double ts(); ///<TODO

        /**
        * @brief    get frequency of IMU.
        * @return    int     frequency of IMU
        */
        int freq();

        /**
        * @brief    get resampled frequency of IMU.
        * @return    int
        */
        int resampled_freq();
        /**
        * @brief    get    unit of gyro.
        * @return    UNIT_TYPE    unit of gyro
        */
        UNIT_TYPE GyroUnit();

        /**
        * @brief    get    unit of accelerometer.
        * @return    UNIT_TYPE    unit of accelerometer
        */
        UNIT_TYPE AcceUnit();

        /**
        * @brief    get    unit of attitude.
        * @return    UNIT_TYPE    unit of attitude
        */
        UNIT_TYPE AttUnit();

        /**
         * @brief
         *
         * @return UNIT_TYPE
         */
        UNIT_TYPE MagUnit(); ///< TODO

        /**
        * @brief    get    imu data order.
        * @return    string    imu data order
        */
        string order();

        /**
        * @brief    get    initial position.
        * @return    Eigen::Vector3d    initial position
        */
        Eigen::Vector3d pos();

        /**
        * @brief    get    initial velocity.
        * @return    Eigen::Vector3d    initial velocity
        */
        Eigen::Vector3d vel();

        /**
        * @brief    get    initial attitude.
        * @return    Eigen::Vector3d    initial attitude
        */
        Eigen::Vector3d att();

        /**
        * @brief    check whether use alignment.
        * @return
            @retval true    use alignment
            @retval false    do not use alignment
        */
        bool align();

        /**
         * @brief
         *
         * @return true
         * @return false
         */
        bool out_intsec(); ///< TODO

        /**
         * @brief
         *
         * @return int
         */
        int out_freq(); ///< TODO

        /**
         * @brief
         *
         * @return int
         */
        int cps(); ///< TODO

        /**
        * @brief    get    subsample number.
        * @return    int    subsample number
        */
        int subsample();

        /**
        * @brief    get    start time.
        * @return    double    start time
        */
        double start();

        /**
        * @brief    get    end time.
        * @return    double    end time
        */
        double end();

        /**
        * @brief    get    time of alignment.
        * @return    double    time of alignment
        */
        double align_time();

        /**
        * @brief    get    displacement of alignment.
        * @return    double    displacement of alignment
        */
        double pos_dist(); //"pos" align threshold

        /**
        * @brief    get    velocity of alignment.
        * @return    double    velocity of alignment
        */
        double vel_norm(); //"vel" align threshold

        /**
        * @brief    get    INS alignment method.
        * @return    ALIGN_TYPE    INS alignment method
        */
        ALIGN_TYPE align_type();

        /**
        * @brief    get    motion mode.
        * @return    vector<string>    motion mode
        */
        vector<string> motion();

        /**
        * @brief    get    constant accelerator bias. for simulation use.
        * @return    Eigen::Vector3d    constant accelerator bias. for simulation use
        */
        Eigen::Vector3d acce_bias(); /**< constant accelerator bias. for simulation use. */

        /**
        * @brief    get    constant gyroscope bias. for simulation use.
        * @return    Eigen::Vector3d    constant gyroscope bias. for simulation use
        */
        Eigen::Vector3d gyro_bias(); /**< constant gyroscope bias. for simulation use. */

        /**
         * @brief
         *
         * @return Eigen::Vector3d
         */
        Eigen::Vector3d imu_installation_rotation();

        /**
         * @brief
         *
         * @return Eigen::Vector3d
         */
        Eigen::Vector3d imu_installation_translation();

    protected:
        double _ts;                       ///< sample interval.
        int _freq;                        ///< sample freq.
        string _cps;                      ///< compensation mode.
        int _subsample;                   ///< subsample number.
        string _order;                    ///< imu data order.
        string _out_order;                ///< output data format
        int _out_freq;                    ///< output data format
        bool _int_sec;                    ///< output data format
        string _GyroUnit;                 ///< Gyro data unit.
        string _AcceUnit;                 ///< Acce data unit.
        Eigen::Vector3d _pos, _vel, _att; ///< pos is XYZ,vel is Vxyz,att is P,R,Y(deg).
        double _start, _end;              ///< start time and end time.
        double _align_time;               ///< coarse align time.
        string _align_type;               ///< coarse align type.
        Eigen::Vector3d _acce_bias;       ///< constant accelerator bias. for simulation use.
        Eigen::Vector3d _gyro_bias;       ///< constant gyroscope bias. for simulation use.
    };


    enum class IMU_TYPE
    {
        NovAtel_SPAN_FSAS,
        NovAtel_SPAN_CPT,
        NovAtel_SPAN_uIRS,
        NovAtel_SPAN_HG9900,
        NovAtel_SPAN_AG11or58,
        NovAtel_SPAN_CPTorKVH,
        NovAtel_SPAN_HG1900,
        NovAtel_SPAN_HG1930,
        NovAtel_SPAN_HG4930,
        NovAtel_SPAN_ADIS16488,
        NovAtel_SPAN_LCI100C,
        NovAtel_SPAN_STIM300,
        NovAtel_SPAN_KVH1750,
        NovAtel_SPAN_uIMU,

        Navigation_Grade,
        Tactical_Grade,
        MEMS_Grade,

        HG1900,
        HG1930,

        ADIS16470,
        StarNeto,
        Customize
    };
    /**
    * @brief change string to UNIT_TYPE.
    * @param[in]    s            imu data unit in string form
    * @return        UNIT_TYPE    imu data unit
    */
    LibGREAT_LIBRARY_EXPORT UNIT_TYPE str2Unit(const string& s);

    /**
    * @brief change string to ALIGN_TYPE.
    * @param[in]    s            INS align method format in string form
    * @return        ALIGN_TYPE    INS align method format
    */
    LibGREAT_LIBRARY_EXPORT ALIGN_TYPE str2align(const string& s);

    /**
    * @brief change string to CPS_TYPE.
    * @param[in]    s            compensation mode in string form
    * @return        CPS_TYPE    compensation mode format
    */
    LibGREAT_LIBRARY_EXPORT CPS_TYPE str2cps(const string& s);

    struct ErrorModel
    {
        Eigen::Vector3d AttInitialSTD;      // deg;
        Eigen::Vector3d VelInitialSTD;      // m/s
        Eigen::Vector3d PosInitialSTD;      // m
        Eigen::Vector3d GyroBiasInitialSTD; // deg/h
        Eigen::Vector3d AcceBiasInitialSTD; // mg

        Eigen::Vector3d AttProcNoisePSD;      // dpsh
        Eigen::Vector3d VelProcNoisePSD;      // mgpsHz
        Eigen::Vector3d PosProcNoisePSD;      // mpsh
        Eigen::Vector3d GyroBiasProcNoisePSD; // dphpsh
        Eigen::Vector3d AcceBiasProcNoisePSD; // mgpsh
    };

    LibGREAT_LIBRARY_EXPORT map<IMU_TYPE, ErrorModel> IMUErrorModels();

}

#endif