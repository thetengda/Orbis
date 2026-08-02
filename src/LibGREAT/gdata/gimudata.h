/**
 * @file         gimudata.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        IMU data structure for storing imu data
 * @version      1.0
 * @date         2025-01-01
 *
 * @copyright Copyright (c) 2025, Wuhan University. All rights reserved.
 *
 */

#ifndef GIMUDATA_H
#define GIMUDATA_H

#include <queue>
#include <stack>
#include "gdata/gdata.h"
#include "gins/gutility.h"
#include "gexport/ExportLibGREAT.h"
using namespace gnut;

namespace great
{

    /**
    * @brief t_gimudata Class for storing the IMU data
    *
    * t_gimudata is used for sins mechanical arrangement.
    * The data contains information about gyro and acce observation.
    * The gdata t_gimudata corresponding to the gcoder t_imufile.
    */
    class LibGREAT_LIBRARY_EXPORT t_gimudata : public t_gdata
    {
    public:
        /** 
        * @brief default constructor. 
        */
        explicit t_gimudata();

        /** 
        * @brief default destructor. 
        */
        explicit t_gimudata(t_spdlog spdlog);

        virtual ~t_gimudata() {}

        /**
        * @brief add one data into vector
        * @note data includes time and wm and vm
        * @param[in] 't' is time
        * @param[in] 'wm' is Angular increment(rad)
        * @param[in] 'vm' is Velocity increment(m/s)
        * @return value reprents success or failure of add
        */
        int add_IMU(const double& t, const Eigen::Vector3d& wm, const Eigen::Vector3d& vm);

        int add_IMU(const double& t, const Eigen::Vector3d& wm, const Eigen::Vector3d& vm, const Eigen::Vector3d& mm);

        /**
        * @brief set imu data interval
        * @param[in] 'ts' is interval
        * @return void
        */
        void set_ts(double ts);

        /**
        * @brief set backward mechanical arrangement
        * @param[in] 'b' is bool represents whether to perform backward mechanical arrangement
        * @return void
        */
        void set_backward(bool b) { _backward = b; }

        /**
        * @brief load data into wm and vm and t accroding to subsamples 
        * @param[in] 't' is last time
        * @param[in] 'nSamples' is subsamples
        * @param[in] '_beg_end' is the direction
        * @param[out] 't' is current time
        * @param[out] 'wm' is Angular increment
        * @param[out] 'vm' is Velocity increment
        * @param[out] 'ts' is computing intervals
        * @return bool
        */
        bool load(vector<Eigen::Vector3d>& wm, vector<Eigen::Vector3d>& vm, double& t, double& ts, int nSamples,bool _beg_end = true);

        /**
        * @brief erase imu data before t
        * @param[in] 't' is current time
        * @param[in] '_beg_end' is the direction
        * return t_gtime represent the most recent undeleted epoch time
        */
        t_gtime erase_bef(t_gtime t, bool _beg_end = true);


        /**
        * @brief return imu data vector size
        * @param[in] '_beg_end' is the direction
        * @return int represent the size of imu data
        */
        int size(bool _beg_end = true);

        /**
        * @brief judge whether IMU data is available
        * @param[in]  now         current epoch
        * @param[in] '_beg_end' is the direction
        * @return
            @retval true   available
            @retval false  unavailable
        */
        virtual bool available(const t_gtime& now, bool _beg_end = true);

        /**
        * @brief get the begin time of imu data
        * @param[in] '_beg_end' is the direction
        * @return double  first epoch
        */
        virtual double beg_obs(bool _beg_end = true); 

        /**
        * @brief get the end time of imu data
        * @param[in] '_beg_end' is the direction
        * @return double  end epoch
        */
        virtual double end_obs(bool _beg_end = true);
     
        virtual int interpolate(const double& intv);

    private:
        /**
        * @struct dataIMU
        * @brief describe imu data information
        */
        struct dataIMU
        {
            double t; /**< time */
            Eigen::Vector3d wm, vm, mm;
        };

        double _ts;                     /// imu interval 
        deque<dataIMU> _imu_back;    /// imu date stored for backward processing
        deque<dataIMU> _imu_forward;         /// imu date stored for backward processing
        int first, end;                 /// imu starting and ending time 
        bool _backward = false;
    };

}



#endif