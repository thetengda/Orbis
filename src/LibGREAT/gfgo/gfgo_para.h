#ifndef GFO_PARA_H
#define GFO_PARA_H
/**
 * @file         gfgo_para.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        Definition of  common parameters for factor graph optimization.
 * @version      1.0
 * @date         2025-01-01
 *
 * @copyright Copyright (c) 2025, Wuhan University. All rights reserved.
 *
 */
#include"gset/gsetbase.h"
#include "gset/gsetfgo.h"
#include <Eigen/Dense>
#include "gexport/ExportLibGREAT.h"
#define WINDOW_SIZE 20
#define NUM_OF_F 1000
#define GWINDOW_SIZE 20  //for GNSS epoch sliding window
#define GINS_WINDOW_SIZE 20  //for GNSS/INS  sliding window
#define NUM_OF_ARC 10000  //for ambiguity arc 

namespace gfgo
{
	/**
	*@brief parameter class for factor graph optimization
	*/
	class LibGREAT_LIBRARY_EXPORT t_gfgo_para
	{
	public:

		explicit t_gfgo_para(gnut::t_gsetbase* gset);
		~t_gfgo_para();
		
	protected:
		int _num_of_cam=2;					///< number of cameras
		int _stereo=0;						///< sign of whether it is a binocular camera  
		int _imu_enable=1;					///< flags that enable the use of imu data
		int _lidar_enable=0;				///< flags that enable the use of lidar data
		int _gnss_enable=0;					///< flags that enable the use of gnss data
		int _loop_enable=0;					///< flags that enable the use of loop closure

		int _estimate_extrinsic=1;				///< estimated extrinsic parameters
		int _estimate_td=0;						///< estimated time delay
		double _max_solver_time=0.08;			///< max sloution time
		int _max_num_iterations=10;				///< max number of iteration
		double _pixel_error=1.0;				///< error of pixel
		double _laser_cloud_error=1.0;			///< error of laser cloud

		double _acc_n=0.05;						///< accelerometer measurement noise
		double _acc_w=0.001;						///< accelerometer bias random work noise
		double _gyr_n=0.005;						///< gyroscope measurement noise
		double _gyr_w=0.0001;						///< gyroscope bias random work noise
		Eigen::Vector3d _gravity;				///< gravity

		double _relative_pos_var=0.1;			///< relative position variance
		double _relative_rot_var=0.001;			///< relative rotation variance		

		int window_size = 3;
		int gwindow_size = 3;
		int gins_window_size = 3;

	};

}
#endif
