/**
 * @file         gbase.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        base function
 * @version      1.0
 * @date         2025-01-01
 * @note 		 pos is defined by lattitude(rad),longitude(rad),height(m).
		 		 velocity is defined by Eastern Velocity(m/s),Northern Velocity(m/s),Upper Velocity(m/s),
		 		 attitude is defined by Pitch(rad),Roll(rad),Yaw(rad),
         		 Pitch:[-pi/2,pi/2], Roll:[-pi,pi], Yaw:[-pi,pi], and the Yaw is positive from North to West.
		 		 local Coordinate System is East,North,Up.
		 		 b Coordinate System is Right,Forward,Up.
		 		 The DCM rotation order is Yaw-->Pitch-->Roll.
 *
 * @copyright Copyright (c) 2025, Wuhan University. All rights reserved.
 *
 */

#ifndef GINSBASE_H
#define GINSBASE_H

#include "gutility.h"
#include "gins/gquat.h"
#include "gins/gearth.h"
#include <string>
#include <fstream>
#include <vector>
#include "gexport/ExportLibGREAT.h"
using namespace std;

namespace great
{

	/**
	* @class t_gbase
	* @brief t_gbase Class for some base algorithms
	*
	* t_gbase is used for sins base algorithms.
	* the main function is convert about different class.such as quat,DCM,Eular angle.
	* and including some special transform,such as symmetry,askew and so on.
	*
	*/
	class LibGREAT_LIBRARY_EXPORT t_gbase
	{
	public:
		
		/**
		* @brief att to DCM
		* @attention att's unit is rad
		* @param[in] att
		* @return Eigen::Matrix3d DCM
		*/
		static Eigen::Matrix3d a2mat(const Eigen::Vector3d& att);


		/**
		* @brief DCM to att
		* @attention att's unit is rad
		* @param[in] m
		* @return Eigen::Vector3d Eular angle
		*/
		static Eigen::Vector3d m2att(const Eigen::Matrix3d& m);


		/**
		* @brief att to quaternion
		* @attention att's unit is rad
		* @param[in] att
		* @return t_gquat quat
		*/
		static t_gquat a2qua(const Eigen::Vector3d& att);
		
		/**
		* @brief quaternion to att
		* @attention att's unit is rad
		* @param[in] qnb
		* @return Eigen::Vector3d Eular angle
		*/
		static Eigen::Vector3d q2att(const t_gquat& qnb);


		/**
		* @brief rotation vector to quaternion
		* @param[in] rv   rotation vector
		* @return t_gquat quat
		*/
		static t_gquat rv2q(const Eigen::Vector3d& rv);

		/**
		* @brief quaternion to rotation vector
		* @param[in] q quaternion
		* @return Eigen::Vector3d
		*/
		static Eigen::Vector3d q2rv(const t_gquat& q);

		/**
		* @brief DCM to quaternion
		* @param[in] DCM
		* @return t_gquat
		*/
		static t_gquat m2qua(const Eigen::Matrix3d& Cnb);

		/**
		* @brief quaternion to  DCM
		* @param[in] qnb
		* @return Eigen::Matrix3d
		*/
		static Eigen::Matrix3d q2mat(const t_gquat& qnb);

		/**
		* @brief rotation vector to  DCM
		* @param[in] rv
		* @return Eigen::Matrix3d
		*/
		static Eigen::Matrix3d rv2m(const Eigen::Vector3d& rv);

		/**
		* @brief double vectors determine DCM
		* @param[in] vb1 vb2 vn1 vn2
		* @return Cnb
		*/
		static Eigen::Matrix3d dv2mat(const Eigen::Vector3d& vb1, const Eigen::Vector3d& vb2, 
			const Eigen::Vector3d& vn1, const Eigen::Vector3d& vn2);

		/**
		* @brief Antisymmetric matrix
		* @param[in] v
		* @return Eigen::Matrix3d
		*/
		static Eigen::Matrix3d askew(const Eigen::Vector3d& v);

		/**
		* @brief 4 dimensions Antisymmetric matrix for 3 dimensions vector
		* @param[in] v
		* @return Eigen::Matrix4d
		*/
		static Eigen::Matrix4d m2m4(const Eigen::Vector3d& v);

		/**
		* @brief 4 dimensions Antisymmetric matrix for 3 dimensions vector
		* @note m2m4 another expression
		* @param[in] v
		* @return Eigen::Matrix4d
		*/
		static Eigen::Matrix4d m2m4_(const Eigen::Vector3d& v);	

		/**
		* @brief DCM from n to e
		* @param[in] pos
		* @return Eigen::Matrix4d
		*/
		static Eigen::Matrix3d Cen(const Eigen::Vector3d& pos);

		/**
		* @brief symmetry
		* @param[in] pos
		* @param[out] pos
		* @return
		*/
		static void symmetry(Eigen::MatrixXd& m);

		/**
		* @brief Vector3d products Matrix3d
		* @param[in] vec
		* @param[in] mat
		* @return Eigen::Vector3d
		*/
		static Eigen::Vector3d product(const Eigen::Vector3d& vec, const Eigen::Matrix3d& mat);

		/**
		* @brief Delete matrix i row and col
		* @param[in] M
		* @param[out] M
		* @return
		*/
		static void delrowcol(Eigen::MatrixXd& M, int i);

		/**
		* @brief Delete matrix i row
		* @param[in] M
		* @param[out] M
		* @return
		*/
		static void delrow(Eigen::MatrixXd& M, int i);

		/**
		* @brief Delete matrix i col
		* @param[in] M
		* @param[out] M
		* @return
		*/
		static void delcol(Eigen::MatrixXd& M, int i);

		/**
		* @brief Delete VectorXd i row
		* @param[in] V
		* @param[out] V
		* @return
		*/
		static void delrow(Eigen::VectorXd& V, int i);

	};


}



#endif