#ifndef GUTILITY_H
#define GUTILITY_H

/**
 * @file         gfgo_para.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        common algorithm
 * @version      1.0
 * @date         2025-01-01
 *
 * @copyright Copyright (c) 2025, Wuhan University. All rights reserved.
 *
 */

#include <ctime>
#include <chrono>
#include <cstdlib>
#include <unordered_map>
#include <thread>
#include <mutex>
#include <queue>
#include <iostream>
#include <string>
#include <ceres/ceres.h>
#include <ceres/rotation.h>
#include <Eigen/Dense>
#include "gexport/ExportLibGREAT.h"


#define M_PI 3.1415926535897932384626433832795

namespace gfgo
{

	/**
	* @enum parameter block size
	* @brief different parameter has different size
	*/
	enum SIZE_PARAMETERIZATION
	{
		SIZE_POSE = 7,
		SIZE_SPEEDBIAS = 9,
		SIZE_FEATURE = 1,
		SIZE_AMB=1,
		SIZE_CRD =3,
		SIZE_SPEED =3,

		// for PPP IF:
		//add byzhang
		SIZE_CLK = 1,
		SIZE_TRP = 1,
		SIZE_ISB = 1
	};

	/**
	* @enum state order 
	* @brief  state order in basic imu state vector (15 dimension)
	*/
	enum StateOrder
	{
		O_P = 0,
		O_R = 3,
		O_V = 6,
		O_BA = 9,
		O_BG = 12
	};

	/**
	* @enum noise order
	* @brief  noise order in Q matrix
	*/
	enum NoiseOrder
	{
		O_AN = 0,
		O_GN = 3,
		O_AW = 6,
		O_GW = 9
	};

	/**
	*@brief utility class for factor graph optimization
	*/
	class LibGREAT_LIBRARY_EXPORT t_gfgo_utility
	{
	
	public:
		template <typename Derived>
		static Eigen::Quaternion<typename Derived::Scalar> deltaQ(const Eigen::MatrixBase<Derived> &theta)
		{
			typedef typename Derived::Scalar Scalar_t;

			Eigen::Quaternion<Scalar_t> dq;
			Eigen::Matrix<Scalar_t, 3, 1> half_theta = theta;
			half_theta /= static_cast<Scalar_t>(2.0);
			dq.w() = static_cast<Scalar_t>(1.0);
			dq.x() = half_theta.x();
			dq.y() = half_theta.y();
			dq.z() = half_theta.z();
			return dq;
		}

		template <typename Derived>
		static Eigen::Matrix<typename Derived::Scalar, 3, 3> skewSymmetric(const Eigen::MatrixBase<Derived> &q)
		{
			Eigen::Matrix<typename Derived::Scalar, 3, 3> ans;
			ans << typename Derived::Scalar(0), -q(2), q(1),
				q(2), typename Derived::Scalar(0), -q(0),
				-q(1), q(0), typename Derived::Scalar(0);
			return ans;
		}


		template <typename Derived>
		static Eigen::Quaternion<typename Derived::Scalar> positify(const Eigen::QuaternionBase<Derived> &q)
		{
			//printf("a: %f %f %f %f", q.w(), q.x(), q.y(), q.z());
			//Eigen::Quaternion<typename Derived::Scalar> p(-q.w(), -q.x(), -q.y(), -q.z());
			//printf("b: %f %f %f %f", p.w(), p.x(), p.y(), p.z());
			//return q.template w() >= (typename Derived::Scalar)(0.0) ? q : Eigen::Quaternion<typename Derived::Scalar>(-q.w(), -q.x(), -q.y(), -q.z());
			return q;
		}

		template <typename Derived>
		static Eigen::Matrix<typename Derived::Scalar, 4, 4> Qleft(const Eigen::QuaternionBase<Derived> &q)
		{
			Eigen::Quaternion<typename Derived::Scalar> qq = positify(q);
			Eigen::Matrix<typename Derived::Scalar, 4, 4> ans;
			ans(0, 0) = qq.w(), ans.template block<1, 3>(0, 1) = -qq.vec().transpose();
			ans.template block<3, 1>(1, 0) = qq.vec(), ans.template block<3, 3>(1, 1) = qq.w() * Eigen::Matrix<typename Derived::Scalar, 3, 3>::Identity() + skewSymmetric(qq.vec());
			return ans;
		}

		template <typename Derived>
		static Eigen::Matrix<typename Derived::Scalar, 4, 4> Qright(const Eigen::QuaternionBase<Derived> &p)
		{
			Eigen::Quaternion<typename Derived::Scalar> pp = positify(p);
			Eigen::Matrix<typename Derived::Scalar, 4, 4> ans;
			ans(0, 0) = pp.w(), ans.template block<1, 3>(0, 1) = -pp.vec().transpose();
			ans.template block<3, 1>(1, 0) = pp.vec(), ans.template block<3, 3>(1, 1) = pp.w() * Eigen::Matrix<typename Derived::Scalar, 3, 3>::Identity() - skewSymmetric(pp.vec());
			return ans;
		}


	};

	/**
	*@brief  class for calculating the time interval
	*/
	class LibGREAT_LIBRARY_EXPORT t_tictoc
	{
	public:
		t_tictoc()
		{
			tic();
		}

		void tic()
		{
			start = std::chrono::system_clock::now();
		}

		double toc()
		{
			end = std::chrono::system_clock::now();
			std::chrono::duration<double> elapsed_seconds = end - start;
			return elapsed_seconds.count() * 1000;
		}

	private:
		std::chrono::time_point<std::chrono::system_clock> start, end;
	};

	
}

#endif