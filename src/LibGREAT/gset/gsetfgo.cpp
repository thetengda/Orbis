/**
 * @file         gsetfgo.cpp
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        control set from XML for FGO
 * @version      1.0
 * @date         2025-11-04
 *
 * @copyright Copyright (c) 2024, Wuhan University. All rights reserved.
 *
 */

#include "gsetfgo.h"

namespace gfgo
{
	t_gsetfgo::t_gsetfgo()
	{
		_set.insert(XMLKEY_FGO);
	}

	t_gsetfgo::~t_gsetfgo()
	{
	}

	void t_gsetfgo::check()
	{
	}

	void t_gsetfgo::help()
	{
	}

	double t_gsetfgo::acc_n()
	{
		_gmutex.lock();
		string x = _doc.child(XMLKEY_ROOT).child(XMLKEY_FGO).child_value("acc_n");
		str_erase(x);
		double x_double = 0.0;
		if (x != "")  x_double = std::stod(x);

		if (double_eq(x_double, 0.0))
		{
			_gmutex.unlock();
			return 0.05;	//default value
		}
		else
		{
			_gmutex.unlock();
			return x_double;
		}
	}

	double t_gsetfgo::acc_w()
	{
		_gmutex.lock();
		string x = _doc.child(XMLKEY_ROOT).child(XMLKEY_FGO).child_value("acc_w");
		str_erase(x);
		double x_double = 0.0;
		if (x != "")  x_double = std::stod(x);

		if (double_eq(x_double, 0.0))
		{
			_gmutex.unlock();
			return 0.001;	//default value
		}
		else
		{
			_gmutex.unlock();
			return x_double;
		}
	}

	double t_gsetfgo::gyr_n()
	{
		_gmutex.lock();
		string x = _doc.child(XMLKEY_ROOT).child(XMLKEY_FGO).child_value("gyr_n");
		str_erase(x);
		double x_double = 0.0;
		if (x != "")  x_double = std::stod(x);

		if (double_eq(x_double, 0.0))
		{
			_gmutex.unlock();
			return 0.005;	//default value
		}
		else
		{
			_gmutex.unlock();
			return x_double;
		}
	}

	double t_gsetfgo::gyr_w()
	{
		_gmutex.lock();
		string x = _doc.child(XMLKEY_ROOT).child(XMLKEY_FGO).child_value("gyr_w");
		str_erase(x);
		double x_double = 0.0;
		if (x != "")  x_double = std::stod(x);

		if (double_eq(x_double, 0.0))
		{
			_gmutex.unlock();
			return 0.0001;	//default value
		}
		else
		{
			_gmutex.unlock();
			return x_double;
		}
	}

	Eigen::Vector3d t_gsetfgo::gravity()
	{
		_gmutex.lock();
		string x = _doc.child(XMLKEY_ROOT).child(XMLKEY_FGO).child_value("gravity_norm");
		str_erase(x);
		double x_double = 0.0;
		if (x != "")  x_double = std::stod(x);

		if (double_eq(x_double, 0.0))
		{
			_gmutex.unlock();
			return Eigen::Vector3d(0, 0, 9.81);	//default value
		}
		else
		{
			_gmutex.unlock();
			return Eigen::Vector3d(0, 0, x_double);
		}
	}

	double t_gsetfgo::pixel_error()
	{
		_gmutex.lock();
		string x = _doc.child(XMLKEY_ROOT).child(XMLKEY_FGO).child_value("pixel_error");
		str_erase(x);
		double x_double = 0.0;
		if (x != "")  x_double = std::stod(x);

		if (double_eq(x_double, 0.0))
		{
			_gmutex.unlock();
			return 15.0;	//default value
		}
		else
		{
			_gmutex.unlock();
			return x_double;
		}
	}
	double t_gsetfgo::laser_cloud_error()
	{
		_gmutex.lock();
		string x = _doc.child(XMLKEY_ROOT).child(XMLKEY_FGO).child_value("laser_cloud_error");
		str_erase(x);
		double x_double = 0.0;
		if (x != "")  x_double = std::stod(x);

		if (double_eq(x_double, 0.0))
		{
			_gmutex.unlock();
			return 1.0;	//default value
		}
		else
		{
			_gmutex.unlock();
			return x_double;
		}
	}
	int t_gsetfgo::imu_enable()
	{
		_gmutex.lock();
		string x = _doc.child(XMLKEY_ROOT).child(XMLKEY_FGO).child_value("imu_enable");
		str_erase(x);
		int x_int = 0;
		if (x != "")  x_int = std::stoi(x);

		if (x_int != 0)
		{
			_gmutex.unlock();
			return 1;
		}
		else
		{
			_gmutex.unlock();
			return 0;
		}
	}

	int t_gsetfgo::lidar_enable()
	{
		_gmutex.lock();
		string x = _doc.child(XMLKEY_ROOT).child(XMLKEY_FGO).child_value("lidar_enable");
		str_erase(x);
		int x_int = 0;
		if (x != "")  x_int = std::stoi(x);

		if (x_int != 0)
		{
			_gmutex.unlock();
			return 1;
		}
		else
		{
			_gmutex.unlock();
			return 0;
		}
	}

	int t_gsetfgo::gnss_enable()
	{
		_gmutex.lock();
		string x = _doc.child(XMLKEY_ROOT).child(XMLKEY_FGO).child_value("gnss_enable");
		str_erase(x);
		int x_int = 0;
		if (x != "")  x_int = std::stoi(x);

		if (x_int != 0)
		{
			_gmutex.unlock();
			return 1;
		}
		else
		{
			_gmutex.unlock();
			return 0;
		}
	}

	int t_gsetfgo::loop_enable()
	{
		_gmutex.lock();
		string x = _doc.child(XMLKEY_ROOT).child(XMLKEY_FGO).child_value("loop_enable");
		str_erase(x);
		int x_int = 0;
		if (x != "")  x_int = std::stoi(x);

		if (x_int != 0)
		{
			_gmutex.unlock();
			return 1;
		}
		else
		{
			_gmutex.unlock();
			return 0;
		}
	}

	int t_gsetfgo::num_of_cam()
	{
		_gmutex.lock();
		string x = _doc.child(XMLKEY_ROOT).child(XMLKEY_FGO).child_value("num_of_cam");
		str_erase(x);
		int x_int = 0;
		if (x != "")  x_int = std::stoi(x);
		_gmutex.unlock();
		return x_int;
	}

	int t_gsetfgo::window_size()
	{
		_gmutex.lock();
		string x = _doc.child(XMLKEY_ROOT).child(XMLKEY_FGO).child_value("window_size");
		str_erase(x);
		int x_int = 0;
		if (x != "")  x_int = std::stoi(x);

		if (x_int <= 0)
		{
			_gmutex.unlock();
			return 10;	//default value
		}
		else
		{
			_gmutex.unlock();
			return x_int;
		}
	}

	int t_gsetfgo::gwindow_size()
	{
		_gmutex.lock();
		string x = _doc.child(XMLKEY_ROOT).child(XMLKEY_FGO).child_value("gnss_window_size");
		str_erase(x);
		int x_int = 0;
		if (x != "")  x_int = std::stoi(x);

		if (x_int <= 0)
		{
			_gmutex.unlock();
			return 10;	//default value
		}
		else
		{
			_gmutex.unlock();
			return x_int;
		}
	}

	int t_gsetfgo::gins_window_size()
	{
		_gmutex.lock();
		string x = _doc.child(XMLKEY_ROOT).child(XMLKEY_FGO).child_value("gins_window_size");
		str_erase(x);
		int x_int = 0;
		if (x != "")  x_int = std::stoi(x);

		if (x_int <= 0)
		{
			_gmutex.unlock();
			return 10;	//default value
		}
		else
		{
			_gmutex.unlock();
			return x_int;
		}
	}

	int t_gsetfgo::estimate_extrinsic()
	{
		_gmutex.lock();
		string x = _doc.child(XMLKEY_ROOT).child(XMLKEY_FGO).child_value("estimate_extrinsic");
		str_erase(x);
		int x_int = 0;
		if (x != "")  x_int = std::stoi(x);

		if (x_int != 0)
		{
			_gmutex.unlock();
			return 1;
		}
		else
		{
			_gmutex.unlock();
			return 0;
		}
	}

	int t_gsetfgo::estimate_td()
	{
		_gmutex.lock();
		string x = _doc.child(XMLKEY_ROOT).child(XMLKEY_FGO).child_value("estimate_td");
		str_erase(x);
		int x_int = 0;
		if (x != "")  x_int = std::stoi(x);

		if (x_int != 0)
		{
			_gmutex.unlock();
			return 1;
		}
		else
		{
			_gmutex.unlock();
			return 0;
		}
	}

	double t_gsetfgo::max_solver_time()
	{
		_gmutex.lock();
		string x = _doc.child(XMLKEY_ROOT).child(XMLKEY_FGO).child_value("max_solver_time");
		str_erase(x);
		double x_double = 0.0;
		if (x != "")  x_double = std::stod(x);

		if (double_eq(x_double, 0.0))
		{
			_gmutex.unlock();
			return 0.08;	//default value
		}
		else
		{
			_gmutex.unlock();
			return x_double;
		}
	}

	int t_gsetfgo::max_num_iterations()
	{
		_gmutex.lock();
		string x = _doc.child(XMLKEY_ROOT).child(XMLKEY_FGO).child_value("max_num_iterations");
		str_erase(x);
		int x_int = 0;
		if (x != "")  x_int = std::stoi(x);

		if (x_int <= 0)
		{
			_gmutex.unlock();
			return 10;	//default value
		}
		else
		{
			_gmutex.unlock();
			return x_int;
		}
	}

	double t_gsetfgo::relative_pos_var()
	{
		_gmutex.lock();
		string x = _doc.child(XMLKEY_ROOT).child(XMLKEY_FGO).child_value("relative_pos_var");
		str_erase(x);
		double x_double = 0.0;
		if (x != "")  x_double = std::stod(x);

		if (double_eq(x_double, 0.0))
		{
			_gmutex.unlock();
			return 0.1;	//default value
		}
		else
		{
			_gmutex.unlock();
			return x_double;
		}
	}

	double t_gsetfgo::relative_rot_var()
	{
		_gmutex.lock();
		string x = _doc.child(XMLKEY_ROOT).child(XMLKEY_FGO).child_value("relative_rot_var");
		str_erase(x);
		double x_double = 0.0;
		if (x != "")  x_double = std::stod(x);

		if (double_eq(x_double, 0.0))
		{
			_gmutex.unlock();
			return 0.001;	//default value
		}
		else
		{
			_gmutex.unlock();
			return x_double;
		}
	}	
}