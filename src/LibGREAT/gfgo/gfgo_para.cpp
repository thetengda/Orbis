/**
 * @file         gfgo_para.cpp
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        Definition of  common parameters for factor graph optimization.
 * @version      1.0
 * @date         2025-01-01
 *
 * @copyright Copyright (c) 2025, Wuhan University. All rights reserved.
 *
 */

#include "gfgo_para.h"

namespace gfgo
{

	t_gfgo_para::t_gfgo_para(gnut::t_gsetbase * gset)
	{
		_num_of_cam = dynamic_cast<t_gsetfgo*>(gset)->num_of_cam();		
		_imu_enable = dynamic_cast<t_gsetfgo*>(gset)->imu_enable();
		_lidar_enable = dynamic_cast<t_gsetfgo*>(gset)->lidar_enable();
		_gnss_enable = dynamic_cast<t_gsetfgo*>(gset)->gnss_enable();
		_loop_enable = dynamic_cast<t_gsetfgo*>(gset)->loop_enable();
		_estimate_extrinsic = dynamic_cast<t_gsetfgo*>(gset)->estimate_extrinsic();
		_estimate_td = dynamic_cast<t_gsetfgo*>(gset)->estimate_td();
		_max_solver_time = dynamic_cast<t_gsetfgo*>(gset)->max_solver_time();
		_max_num_iterations = dynamic_cast<t_gsetfgo*>(gset)->max_num_iterations();
		_pixel_error = dynamic_cast<t_gsetfgo*>(gset)->pixel_error();
		_laser_cloud_error = dynamic_cast<t_gsetfgo*>(gset)->laser_cloud_error();
		if (_num_of_cam == 2) _stereo = 1;
		else _stereo = 0;
		
		_acc_n = dynamic_cast<t_gsetfgo*>(gset)->acc_n();
		_acc_w = dynamic_cast<t_gsetfgo*>(gset)->acc_w();
		_gyr_n = dynamic_cast<t_gsetfgo*>(gset)->gyr_n();
		_gyr_w = dynamic_cast<t_gsetfgo*>(gset)->gyr_w();
		_gravity.setZero();
		_gravity = dynamic_cast<t_gsetfgo*>(gset)->gravity();

		_relative_pos_var = dynamic_cast<t_gsetfgo*>(gset)->relative_pos_var();
		_relative_rot_var = dynamic_cast<t_gsetfgo*>(gset)->relative_rot_var();	

		window_size = dynamic_cast<t_gsetfgo*>(gset)->window_size();
		gwindow_size = dynamic_cast<t_gsetfgo*>(gset)->gwindow_size();
		gins_window_size= dynamic_cast<t_gsetfgo*>(gset)->gins_window_size();
		
	}

	t_gfgo_para::~t_gfgo_para()
	{

	}

}

