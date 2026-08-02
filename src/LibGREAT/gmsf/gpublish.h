/**
 * @file         gpublish.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        visualization
 * @version      1.0
 * @date         2025-01-01
 *
 * @copyright Copyright (c) 2025, Wuhan University. All rights reserved.
 *
 */

#ifndef GPUBLISH_H
#define GPUBLISH_H

#include "gviewer.h"
#include "gins/gutility.h"
#include "gins/gbase.h"
#include "gexport/ExportLibGREAT.h"
#include <Eigen/Eigen>
using namespace std;

namespace great
{

    /**
    * @class t_gpublish
    * @brief t_gpublish Class is used for trajectory display
    */
    class LibGREAT_LIBRARY_EXPORT t_gpublish
    {
    public:
        /** @brief default constructor. */
        t_gpublish();

        /** @brief default destructor. */
        ~t_gpublish();

        /**
        * @brief initialize the viewer for display
        */
        void Initialize();

        /**
        * @brief update the parameters for trajectory display
        *
        * @param[in] imu_state        store the imu information to show on the viewer
        * @param[in] init_campos    the initial value of the camera position
        */
        void UpdateNewState(const IMUState& imu_state);
        void AddMapPoints(const vector<Eigen::Vector3d> & map_points);
        void UpdateVisualPoints(const vector<Eigen::Vector3d> & map_points);
        void UpdatePlanePoints(const vector<Eigen::Vector3d> & pcs, const vector<vector<Eigen::Vector3d>>& _near_points, Eigen::Matrix3d R_l_e, Eigen::Vector3d t_l_e);
        
    private:
        
        gviewer *viewer=nullptr;                // the core class of trajectory display
        bool firstFlag_imu = true;        // A indicator to determine the value of _init_imupos
        bool firstFlag_cam = true;        // A indicator to determine the value of _init_campos,
        Eigen::Vector3d _init_campos;    // transfer the zero of the IMU coordinate
        Eigen::Vector3d _init_imupos;    // transfer the zero of the IMU coordinate
        Eigen::Matrix3d R_e_n;          // rotation from ECEF frame to ENU frame for first imudata

        
    };

}
#endif