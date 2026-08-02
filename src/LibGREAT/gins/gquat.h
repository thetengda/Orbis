/**
 * @file         gquat.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        Class for Quaternion
 * @version      1.0
 * @date         2025-01-01

 * @note         Acknowledgement to Prof. Gongmin Yan, Northwestern Polytechnical University
 * @cite         PSINS, Precise Strapdown Inertial Navigation System (https://psins.org.cn/)
 *
 * @copyright Copyright (c) 2025, Wuhan University. All rights reserved.
 *
 */

#ifndef GQUAT_H
#define GQUAT_H

#include "gutility.h"
#include <string>
#include <fstream>
#include <vector>
#include "gexport/ExportLibGREAT.h"

using namespace std;

namespace great
{

    /**
    * @class t_gquat
    * @brief t_gquat Class for Quaternion.
    * 
    * t_gquat is Hamilton Quaternion Class
    * the main function is operator overriding,such as +,-,*.
    * and including some special algorithm,such as normlize,conjugation and so on.
    */
    class LibGREAT_LIBRARY_EXPORT t_gquat
    {
    public:

        /**
        * @brief Constructor
        * @note default value (1.0,0.0,0.0,0.0)
        */
        t_gquat(double q0 = 1.0, double q1 = 0.0, double q2 = 0.0, double q3 = 0.0);
        
        /**
        * @brief constructor.
        * @param[in]  m  Eigen:Vector4d initialize
        */
        t_gquat(const Eigen::Vector4d& m);

        /**
        * @brief + operator overloading
        * @param[in] q
        * @return t_gquat
        */
        t_gquat operator+(const t_gquat &q) const;


        /**
        * @brief + operator overloading
        * @note Quaternion add misalignment angle phi to true Quaternion
        * @param[in] phi    misalignment angle
        * @return t_gquat   true Quaternion
        */
        t_gquat operator+(const Eigen::Vector3d& phi)const;


        /**
        * @brief - operator overloading
        * @note Quaternion delete misalignment angle phi
        * @param[in] phi      misalignment angle phi
        * @return t_gquat
        */
        t_gquat operator-(const Eigen::Vector3d &phi) const;


        /**
        * @brief * operator overloading
        */
        t_gquat operator*(const t_gquat &q) const;
        Eigen::Vector3d operator-(const t_gquat &quat) const;

        /**
        * @brief * operator overloading
        * @note Quaternion multuplies angle 
        * @param[in] v
        * @return Eigen::Vector3d
        */
        Eigen::Vector3d operator*(const Eigen::Vector3d &v) const;


        /**
        * @brief *= operator overloading
        */
        t_gquat& operator*=(const t_gquat &q);

        /**
        * @brief Normlize
        * @return
        */
        static void normlize(t_gquat& q);

        /**
        * @brief Conjugation
        * @note
        * @param[in] q
        * @return t_gquat
        */


        static t_gquat conj(const t_gquat& q);

        Eigen::Matrix4d left();

        Eigen::Matrix4d right();

        Eigen::Vector3d vec();


    public:
        double q0, q1, q2, q3;        /// Quaternion value
    };

}



#endif