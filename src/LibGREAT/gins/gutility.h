/**
 * @file         gutility.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        implement common algorithm and define some enums
 * @version      1.0
 * @date         2025-01-01
 *
 * @copyright Copyright (c) 2025, Wuhan University. All rights reserved.
 *
 */

#ifndef UTILITY_H
#define UTILITY_H

#include "gexport/ExportLibGREAT.h"
#include <string>
#include <Eigen/Eigen>
#include "gutils/gsysconv.h"

#define UNICORE_ACC_SCALE 400 / pow(2, 31)				///< m/s/LSB
#define UNICORE_GYRO_SCALE 2160 / pow(2, 31)			///< deg/LSB
#define GREAT_ACC_SCALE          400 / pow(2,31)	    ///< m/s/LSB
#define GREAT_GYRO_SCALE         2160 / pow(2,31)		///< deg/LSB
#define STARNETO_ACC_SCALE 0.05 / pow(2, 15)
#define STARNETO_GYRO_SCALE 0.1 / (3600 * 256.0)
#define STARNETO_G 9.80144145
#define IMR_FSAS_GYRO_SCALE 0.1 / pow(2, 8)				///< arcsec/LSB
#define IMR_FSAS_ACCE_SCALE 0.05 / pow(2, 15)			///< m/s/LSB

namespace great
{

	typedef long long int IMUStateIDType;
	/**
	* @struct IMUState
	* @brief some information only used in visualization
	*/
	struct IMUState
	{
		IMUStateIDType id;                // id of imu state,not use
		double time;                    // timstamp of imu state,not use
		Eigen::Quaterniond orientation; // attitude of imu state
		Eigen::Vector3d position;        // position of imu state

		IMUState() : id(0), time(0),
			orientation(Eigen::Quaterniond::Identity()),
			position(Eigen::Vector3d::Zero()) {}
	};

	/**
	* @class t_gglv
	* @brief define some global const values,such as PI, deg et al
	*/
	//class LibGREAT_LIBRARY_EXPORT t_gglv;

	////LibSet_LIBRARY_EXPORT extern const t_gglv glv;
	//

	//LibGREAT_LIBRARY_EXPORT extern const t_gglv glv;
	class t_gglv
	{
	public:
		//inline t_gglv(double Re = 6378137.0, double f = (1.0 / 298.257), double wie0 = 7.2921151467e-5, double g0 = 9.7803267714);
		
		//double Re, f, wie, g0, mg, ug;	///< earth ellipsoid parameters 
		static const double Re;              ///< semi major axes [m]
		static const double f;               ///< inverse flattening [-]
		static const double wie;             ///< self-rotation
		static const double g0;              ///< equatorial gravity
		static const double mg;              ///< mg
		static const double ug;              ///< ug
		static const double PI;              ///< PI
		static const double EPS;             ///< Epsilon
		static const double INF;             ///< Infinity
		static const double deg;             ///< 1 degree in radians
		static const double min;             ///< 1 minute in radians
		static const double sec;             ///< 1 second in radians
		static const double ppm;             ///< part per million
		static const double hur;             ///< 1 h = 3600.0 s
		static const double dps;             ///< deg per second
		static const double dpss;            ///< deg per sqrt(second)
		static const double dph;             ///< deg per hour
		static const double dpsh;            ///< deg per sqrt(hour)
		static const double dphpsh;          ///< dph per sqrt(hour)
		static const double mpsh;            ///< m per hour
		static const double mpspsh;          ///< m per sqrt(hour)
		static const double ppmpsh;          ///< ppm per sqrt(hour)
		static const double secpsh;          ///< second per sqrt(hour)
		static const double mgpsHz;
		static const double ugpsHz;
		static const double mgpsh;
		static const double ugpsh;
		//double mgpsHz, ugpsHz, mgpsh, ugpsh;	///<TODO
	};
	/**
	* @enum FLT_TYPE
	* @brief filter mode,included forward,backward,forward and backward,RTS
	*/
	enum FLT_TYPE
	{
		FORWARD,  ///< forward mode
		BACKWARD, ///< backward mode
		FBS,      ///< forward and backward mode
		RTS       ///< RTS mode
	};

	/**
	* @enum CPS_TYPE
	* @brief Compensation mode,included Cone,Poly,One plus previous
	*/
	enum CPS_TYPE
	{
		CONE,         ///< Cone mode
		ONE_PLUS_PRE, ///< One plus previous mdoe
		POLY          ///< Poly mode
	};

	/**
	* @enum MOTION_TYPE
	* @brief motion mode
	* @note be used for imu simulation,includes different motion mode
	*/
	enum MOTION_TYPE
	{
		m_keep,       ///< keep mode
		m_accelerate, ///< accelerate mode
		m_decelerate, ///< decelerate mode
		m_yawleft,    ///< yaw left mode
		m_yawright,   ///< yaw right mode
		m_pitchup,    ///< pitch uo mode
		m_pitchdown,  ///< pitch down mode
		m_rollleft,   ///< roll left mode
		m_rollright,  ///< roll right mode
		m_turnleft,   ///< turn left mode
		m_turnright,  ///< turn right mode
		m_climb,      ///< climb mode
		m_descent,    ///< descent mode
		m_s,          ///< S-shaped trajectory mode
		m_8,          ///< 8-shaped trajectory mode
		m_straight,   ///< straight mode
		m_static,     ///< static mode
		m_default     ///< default mode
	};

	/**
	* @enum UNIT_TYPE
	* @brief imu data unit
	*/
	enum UNIT_TYPE
	{
		RAD,  ///< rad
		DEG,  ///< deg
		RPS,  ///< TODO
		DPS,  ///< TODO
		RPH,  ///< TODO
		DPH,  ///< TODO
		MPS,  ///< TODO
		MPS2, ///< TODO
		UNDF  ///< undefined
	};

	/**
	* @enum ALIGN_TYPE
	* @brief INS align method, STC_AGN is static base alignment, MA is moving base alignment
	* and the POS_AGN is position vector alignment, and the VEL_AGN is velocity vector alignment
	*/
	enum ALIGN_TYPE
	{
		AUTO,               ///< auto alignment
		STC_AGN,            ///< static base alignment
		MA,                 ///< moving base alignment
		POS_AGN,            ///< position vector alignment
		VEL_AGN,            ///< velocity vector alignment
		Multi_Ant,
		ALIGN_DEFAULT = 999 ///< default alignment
	};

	/**
	* @struct t_scheme
	* @brief  processing scheme
	*/
	struct LibGREAT_LIBRARY_EXPORT t_scheme
	{
		double t = 0.0;                 ///< current time
		double ts = 0.0;                ///< sampling inteval
		int freq = 0;                   ///< IMU Freq
		double start = 0.0;             ///< IMU start time
		double end = 0.0;               ///< IMU end time
		double delay = 0.0;             ///< the GNSS misalignment time diff with imu.
		double delay_odo = 0.0;         ///< the ODO misalignment time diff with imu.
		double max_pdop = 0.0;          ///< the max PDOP threshold for LC.
		int min_sat = 0;                ///< the minimum satellite number threshold for LC.
		int align_time = 0;             ///< coarse align time.
		int nSamples = 0;               ///< subsample number
		int Cps = 0;                    ///< compensation mode
		bool align = false;             ///< whether self-aligned or not
		bool _imu_scale = false;        ///< whether to estimate IMU factor
		bool _imu_inst_rot = false;     ///< whether to estimate IMU installation rotation
		bool _imu_inst_trans = false;   ///< whether to estimate IMU installation translation
		bool _ZUPT = false;             ///< whether to enable ZUPT
		bool _NHC = false;              ///< whether to enable NHC
		bool _Hgt = false;              ///< whether to enable Height constraint

		t_scheme();
		/**
		* @brief Constructor
		* set parameter value by set
		*/
		t_scheme(gnut::t_gsetbase* set);
	};

	int sign(double d);

	/**
	* @brief Cart to Geod
	* overwrite gnut's transform function
	*/
	Eigen::Vector3d Cart2Geod(const Eigen::Vector3d& CartPos, bool b);

	/**
	* @brief Geod to Cart
	* overwrite gnut's transform function
	*/
	Eigen::Vector3d Geod2Cart(const Eigen::Vector3d& GeodPos, bool b);

	/**
	* @brief Geod to Cart
	* overwrite gnut's transform function
	*/
	Eigen::Vector3d XYZ2ENU(const Eigen::Vector3d& xyz, const Eigen::Vector3d& xyz_ref);

	Eigen::MatrixXd NewMat2Eigen(const Matrix& newmat);
	Matrix Eigen2newMat(const Eigen::MatrixXd& eigen);
	SymmetricMatrix Eigen2BaseMatrix(const Eigen::MatrixXd& eigen);
	Eigen::MatrixXd BaseMatrix2Eigen(const SymmetricMatrix& newmat);
	Eigen::VectorXd Columns2VectorXd(const ColumnVector& newmat);
	ColumnVector VectorXd2Columns(const Eigen::VectorXd& eigen);

}

#endif
