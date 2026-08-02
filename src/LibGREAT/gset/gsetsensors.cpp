/**
 * @file         gsetsensors.cpp
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        control set from XML for MSF based on factor graph optimization
 * @version      1.0
 * @date         2025-11-04
 *
 * @copyright Copyright (c) 2024, Wuhan University. All rights reserved.
 *
 */
#include "gsetsensors.h"
namespace gfgomsf
{
	MSF_TYPE str2msf(const string &s)
	{
		string tmp = s;
		transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);

		if (tmp == "GINS_TC")
			return MSF_TYPE::GINS_TC_MODE;
		if (tmp == "GNSS")
			return MSF_TYPE::GNSS_MODE;
		return GINS_TC_MODE;
	}
	


	void t_gset_gomsf_sensors::check()
	{
		

	}

	void t_gset_gomsf_sensors::help()
	{


	}
	MSF_TYPE t_gset_gomsf_sensors::msf_type()
	{
		_gmutex.lock();
		string res = _doc.child(XMLKEY_ROOT).child(XMLKEY_GOMSF_SENSORS).child_value("msf_type");
		str_erase(res);
		_gmutex.unlock();
		return str2msf(res);	//default GINS_TC_MODE

	}

	string t_gset_gomsf_sensors::result_path()
	{
		_gmutex.lock();
		string x = _doc.child(XMLKEY_ROOT).child(XMLKEY_GOMSF_SENSORS).child_value("result_path");
		str_erase(x);
		_gmutex.unlock();
		return x;
	}
}