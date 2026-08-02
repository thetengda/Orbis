/**
 * @file         gimudata.cpp
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        IMU data structure for storing imu data
 * @version      1.0
 * @date         2025-01-01
 *
 * @copyright Copyright (c) 2025, Wuhan University. All rights reserved.
 *
 */

#include "gimudata.h"
#include "gmodels/ginterp.h"
using namespace Eigen;

great::t_gimudata::t_gimudata():t_gdata()
{
    first = end = 1;
    id_type(ID_TYPE::IMUDATA);
}

great::t_gimudata::t_gimudata(t_spdlog spdlog) : t_gdata(spdlog)
{
    first = end = 1;
    id_type(ID_TYPE::IMUDATA);
}
int great::t_gimudata::add_IMU(const double& t, const Vector3d& wm, const Vector3d& vm)
{

#ifdef BMUTEX   
    boost::mutex::scoped_lock lock(_mutex);
#endif
    _gmutex.lock();

    dataIMU tmp = { t,wm,vm };
    _imu_forward.push_back(tmp);
    if (_backward)_imu_back.push_back(tmp);

    _gmutex.unlock();
    return 0;
}

int great::t_gimudata::add_IMU(const double& t, const Vector3d& wm, const Vector3d& vm, const Vector3d& mm)
{

#ifdef BMUTEX   
    boost::mutex::scoped_lock lock(_mutex);
#endif
    _gmutex.lock();

    dataIMU tmp = { t,wm,vm,mm };
    _imu_forward.push_back(tmp);
    if (_backward)_imu_back.push_back(tmp);
    _gmutex.unlock();
    return 0;
}

void great::t_gimudata::set_ts(double ts)
{
    _ts = ts;
}

bool great::t_gimudata::load(vector<Eigen::Vector3d>& wm, vector<Eigen::Vector3d>& vm, double & t, double & ts, int nSamples,bool _beg_end)
{

#ifdef BMUTEX
    boost::mutex::scoped_lock lock(_mutex);
#endif
    _gmutex.lock();
    wm.clear(); vm.clear();

    if ((_beg_end&&_imu_forward.size() < nSamples) || (!_beg_end&&_imu_back.size() < nSamples)) {
        t = 0; 
        _gmutex.unlock(); 
        return false;
    }
    double _pre_time = t;
    //_pre_time = t;
    for (int j = 0; j < nSamples; j++)
    {
        if (_beg_end)
        {
            t = _imu_forward.front().t;
            wm.push_back(_imu_forward.front().wm);
            vm.push_back(_imu_forward.front().vm);
            _imu_forward.pop_front();
        }
        else
        {
            t = _imu_back.back().t;
            wm.push_back(-_imu_back.back().wm);
            vm.push_back(_imu_back.back().vm);
            _imu_back.pop_back();
        }

    }
    if (_beg_end && first == 1 ) {
        ts = nSamples * _ts;
        first = 0;
    }
    else if (!_beg_end && end == 1) {
        ts = nSamples * _ts;
        end = 0;
    }
    else ts = fabs(t - _pre_time);

    _gmutex.unlock();
    return true;
}

t_gtime great::t_gimudata::erase_bef(t_gtime t, bool _beg_end)
{
	double t_ins;
	double t_gnss = t.sow() + t.dsec();
	t_ins = _imu_forward.front().t;
	while (t_ins < t_gnss)
	{
		_imu_forward.pop_front();
		t_ins = _imu_forward.front().t;
	}
	return t_gtime(t.gwk(), t_ins);
}


int great::t_gimudata::size(bool _beg_end)
{
    if(_beg_end) return _imu_forward.size();
    else return _imu_back.size();
}

bool great::t_gimudata::available(const t_gtime & now, bool _beg_end)
{
    try {
        double t = now.sow() + now.dsec();

        while (_imu_forward.size() < 5) t_gtime::gmsleep(500);
        if ((_imu_forward.size() > 0 && _beg_end && t < _imu_forward.back().t) || (_imu_back.size() > 0 && !_beg_end && t > _imu_back.front().t))
            return true;
        else
            return false;
    }

    catch(...)
    { 
        return false;
    }
}

double great::t_gimudata::beg_obs(bool _beg_end)
{
    if (_beg_end)
    {
        if (_imu_forward.size() == 0)return 0;
        return _imu_forward.front().t;
    }
    else
    {
        if (_imu_back.size() == 0)return 0;
        return _imu_back.front().t;
    }
}

double great::t_gimudata::end_obs(bool _beg_end)
{
    if (_beg_end)
    {
        if (_imu_forward.size() == 0)return 0;
        return _imu_forward.back().t;
    }
    else
    {
        if (_imu_back.size() == 0)return 0;
        return _imu_back.back().t;
    }
}

int great::t_gimudata::interpolate(const double& intv)
{
	if (intv > 0.1)return -1;

	deque<dataIMU> resampled_imu;
	double target_t = int(_imu_forward[0].t) + 1.0;
	double last_orig_t = _imu_forward[0].t;
	double error = 0.0;
	for (int i = 0; i < _imu_forward.size() - 1; i++)
	{
		double next_orig_t = _imu_forward[i + 1].t;
		while (target_t <= next_orig_t)
		{
			dataIMU imu01;
			t_ginterp interp;
			for (auto j = 0; j < 3; j++)
			{
				map<double, double> orig;
				orig.insert(make_pair(_imu_forward[i].t, _imu_forward[i].wm[j] / _ts));
				orig.insert(make_pair(_imu_forward[i + 1].t, _imu_forward[i + 1].wm[j] / _ts));
				if (interp.linear(orig, target_t, imu01.wm[j]) < 0)
				{
					imu01.wm[j] = _imu_forward[i].wm[j] / _ts;
					if (_spdlog)
						SPDLOG_LOGGER_ERROR(_spdlog, "imu interpolation failed between " + to_string(_imu_forward[i].t)
							+ " and " + to_string(_imu_forward[i + 1].t));
				}
				orig.clear();
				orig.insert(make_pair(_imu_forward[i].t, _imu_forward[i].vm[j] / _ts));
				orig.insert(make_pair(_imu_forward[i + 1].t, _imu_forward[i + 1].vm[j] / _ts));
				if (interp.linear(orig, target_t, imu01.vm[j]) < 0)
				{
					imu01.vm[j] = _imu_forward[i].vm[j] / _ts;
					if (_spdlog)
						SPDLOG_LOGGER_ERROR(_spdlog, "imu interpolation failed between " + to_string(_imu_forward[i].t)
							+ " and " + to_string(_imu_forward[i + 1].t));
				}
			}
			imu01.t = target_t; imu01.wm *= intv; imu01.vm *= intv;
			resampled_imu.push_back(imu01);
			double y = intv - error;
			double tmp = target_t + y;
			error = (tmp - target_t) - y;
			target_t = tmp;
		}
		last_orig_t = next_orig_t;
	}
	_imu_forward = resampled_imu;

	return 1;
}
