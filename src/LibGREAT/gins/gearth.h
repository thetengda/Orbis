/**
 * @file         gearth.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        Class for Earth parameter updating
 * @version      1.0
 * @date         2025-01-01
 * @note         Acknowledgement to Prof. Gongmin Yan, Northwestern Polytechnical University
 * @cite         PSINS, Precise Strapdown Inertial Navigation System (https://psins.org.cn/)
 *
 * @copyright Copyright (c) 2025, Wuhan University. All rights reserved.
 *
 */

#ifndef GEARTH_H
#define GEARTH_H

#include "gbase.h"
#include "gset/gsetins.h"
#include "gexport/ExportLibGREAT.h"
//#include "gutility.h"

using namespace std;

namespace great
{

    /**
    * @class t_gearth
    * @brief t_gearth Class for Earth parameter updating.
    *
    * The Earth Parameter,such as e,a,b et al.
    * The main function is Earth Updating.
    */
    class LibGREAT_LIBRARY_EXPORT t_gearth
    {
    public:

        /**
        * @brief constructor.
        *
        * @param[in]  a0   major semi axis  
        * @param[in]  f0   flattening
        * @param[in]  g0   gravity
        */
        explicit t_gearth(double a0 = Aell, double f0 = 1.0 / Finv, double g0 = G_EQUA);


        /**
        * @brief The Earth Parameter Update
        * The Earth Parameter mainly includes sb, sb2, sb4, cb, tb, sl, cl.
        *
        * @param[in]  pos        position of body
        * @param[in]  vn         velocity of body
        */
        void Update(const Eigen::Vector3d &pos, const Eigen::Vector3d &vn);


        /**
        * @brief Compute dpos according to vn and ts
        *
        * @param[in]  vn         velocity of body
        * @param[in]  ts         interval
        * @return 
            Position Increments
        */
        Eigen::Vector3d v2dp(const Eigen::Vector3d& vn, double ts);

    public:
        double a, b;                                         /// major & Minor Axissemi axis 
        double f, e, e2, wie;                                 /// flattening,eccentricity,..,angular of rotation
        double sb, sb2, sb4, cb, tb, sl, cl;                 /// sinB,..,..,cosB,tanB,sinL,cosL
        double RMh, RNh, cbRNh, f_RMh, f_RNh, f_cbRNh;         /// RM+H,RN+H,RMH*cosB
        Eigen::Vector3d pos, vn;                             /// position,velocity
        Eigen::Vector3d wnie, wnen, wnin, weie, gn, gcc;
        Eigen::Matrix3d Cen, Cne;                             /// DCM from n-fram to e-frame

    };


}



#endif