/**
 * @file         gimu.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        IMU core algorithm
 * @version      1.0
 * @date         2025-01-01
 * @note         Acknowledgement to Prof. Gongmin Yan, Northwestern Polytechnical University
 * @cite         PSINS, Precise Strapdown Inertial Navigation System (https://psins.org.cn/)
 * 
 * @copyright Copyright (c) 2025, Wuhan University. All rights reserved.
 *
 */

#ifndef GIMU_H
#define GIMU_H
#include "gbase.h"

namespace great
{

    /**
    * @class t_gimu
    * @brief t_gimu Class is IMU core class for imu updating.
    *
    * t_gimu is used for simu core algorithms,including cone,poly error compensation
    * and the main function is imu info updating.
    */
    class LibGREAT_LIBRARY_EXPORT t_gimu
    {
    public:

        /** @brief default constructor. */
        explicit t_gimu();

        /**
        * @brief Phi cone compensation
        * @param[in] wm                angle increments
        * @return 
            result of cone error.
        */
        Eigen::Vector3d phi_cone_crt(const vector<Eigen::Vector3d>& wm);


        /**
        * @brief Phi poly compensation
        * @param[in] wm                angle increments
        * @return 
            result of poly error.
        */
        Eigen::Vector3d phi_poly_crt(const vector<Eigen::Vector3d>& wm);

        /**
        * @brief Vm cone compensation
        * @param[in] wm                angle increments
        * @param[in] vm                velocity increments
        * @return 
            result of cone error.
        */
        Eigen::Vector3d vm_cone_crt(const vector<Eigen::Vector3d>& wm, const vector<Eigen::Vector3d>& vm);


        /**
        * @brief Vm poly compensation
        * @param[in] wm                angle increments
        * @param[in] vm                velocity increments
        * @return
            result of poly error.
        */
        Eigen::Vector3d vm_poly_crt(const vector<Eigen::Vector3d>& wm, const vector<Eigen::Vector3d>& vm);


        /**
        * @brief IMU error compensation
        *
        * @param[in] wm                angle increments
        * @param[in] vm                velocity increments
        * @param[in] scm            scheme of processing
        */
        void Update(const vector<Eigen::Vector3d>& wm, const vector<Eigen::Vector3d>& vm,const t_scheme& scm);

        
        Eigen::Vector3d phim, dvbm;        /// increment because of poly or cone.

    private:
        Eigen::Vector3d wm_1, vm_1;        /// last increment.
        double pf1, pf2;                /// flag of first time.

    };

}



#endif
