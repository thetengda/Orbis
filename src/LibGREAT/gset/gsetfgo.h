/**
 * @file         gsetfgo.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        control set from XML for FGO
 * @version      1.0
 * @date         2025-11-04
 *
 * @copyright Copyright (c) 2024, Wuhan University. All rights reserved.
 *
 */

#ifndef GSETFGO_H
#define GSETFGO_H

#include "gset/gsetbase.h"
#include "gexport/ExportLibGREAT.h"
#include <Eigen/Eigen>

#define XMLKEY_FGO "fgo"

#define WINDOW_SIZE 10
#define NUM_OF_F 10000

using namespace std;
using namespace gnut;

namespace  gfgo
{
	/**
	* @brief	class for set fgo xml
	*/
	class LibGREAT_LIBRARY_EXPORT t_gsetfgo : public virtual t_gsetbase
	{
	public:

		/** @brief default constructor. */
		t_gsetfgo();

		/** @brief default destructor. */
		~t_gsetfgo();

		/**
		* @brief settings check.
		*/
		void check();

		/**
		* @brief settings help.
		*/
		void help();

		/**
		* @brief  get accelerometer measurement noise standard deviation.
		* @return	double	 accelerometer measurement noise standard deviation
		*/
		double acc_n();

		/**
		* @brief  get accelerometer bias random work noise std.
		* @return	double	 accelerometer measurement noise std
		*/
		double acc_w();

		/**
		* @brief  get gyroscope measurement noise std.
		* @return	double	 gyroscope measurement noise std
		*/
		double gyr_n();

		/**
		* @brief  get gyroscope bias random work noise std.
		* @return	double	 gyroscope bias random work noise std
		*/
		double gyr_w();

		/**
		* @brief  get gravity magnitude.
		* @return	Eigen::Vector3d	 gravity magnitude
		*/
		Eigen::Vector3d gravity();

		/**
		* @brief  get pixel noise std.
		* @return	double	 pixel noise std
		*/
		double pixel_error();

		/**
		* @brief  get laser cloud noise std.
		* @return	double	 laser cloud noise std
		*/
		double laser_cloud_error();

		//control
		/**
		* @brief  check whether use imu data.
		* @return
			@retval 1	use imu data
			@retval 0	do not use imu data
		*/
		int imu_enable();

		/**
		* @brief  check whether use lidar data.
		* @return
			@retval 1	use lidar data
			@retval 0	do not use lidar data
		*/
		int lidar_enable();

		/**
		* @brief  check whether use gnss data.
		* @return
			@retval 1	use gnss data
			@retval 0	do not use gnss data
		*/
		int gnss_enable();	

		/**
		* @brief  check whether use loop closure.
		* @return
			@retval 1	use loop closure
			@retval 0	do not use loop closure
		*/
		int loop_enable();

		/**
		* @brief  check whether use camera data and get the number of cameras.
		* @return
			@retval >0	number of cameras
			@retval 0	do not use camera data
		*/
		int num_of_cam();

		/**
		* @brief  get size of the slide window.
		* @return	int	size of the slide window
		*/
		int window_size();

		/**
	    * @brief  get size of the slide window for FGO-GNSS.
	    * @return    int    size of the slide window
	    */
		int gwindow_size();

		/**
		* @brief  get size of the slide window for FGO-GINS.
		* @return    int    size of the slide window
		*/
		int gins_window_size();

		/**
		* @brief  check whether estimate cam->IMU extrinsic parameters.
		* @return
			@retval 1	estimate cam->IMU extrinsic parameters
			@retval 0	do not estimate cam->IMU extrinsic parameters
		*/
		int estimate_extrinsic();

		/**
		* @brief  check whether estimate extrinsic time delay.
		* @return
			@retval 1	estimate extrinsic time delay
			@retval 0	do not estimate extrinsic time delay
		*/
		int estimate_td();

		/**
		* @brief  get max time of slover.
		* @return	double max time of slover
		*/
		double max_solver_time();

		/**
		* @brief  get max number of iterations.
		* @return	int max number of iterations
		*/
		int max_num_iterations();

		/**
		* @brief  get relative position transform variance.
		* @return	double relative position transform variance
		*/
		double relative_pos_var();

		/**
		* @brief  get relative rotation transform variance.
		* @return	double relative rotation transform variance
		*/
		double relative_rot_var();

	};

}
#endif