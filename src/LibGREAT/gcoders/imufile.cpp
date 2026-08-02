/**
 * @file         imufile.cpp
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        read and decode imu file
 * @version      1.0
 * @date         2025-01-01
 *
 * @copyright Copyright (c) 2025, Wuhan University. All rights reserved.
 *
 */

#include "imufile.h"
#include "gdata/gimudata.h"
#include "gset/gsetins.h"
using namespace Eigen;

great::t_imufile::t_imufile(t_gsetbase* s, string version, int sz) :t_gcoder(s, sz),
_ts(0.005),
_GyroUnit(DPS),
_AcceUnit(MPS2),
_complete(false),
_order("garfu")
{
    _ts = dynamic_cast<t_gsetins*>(s)->ts();
    _freq = dynamic_cast<t_gsetins*>(s)->freq();
    _GyroUnit = dynamic_cast<t_gsetins*>(s)->GyroUnit();
    _AcceUnit = dynamic_cast<t_gsetins*>(s)->AcceUnit();
    _MagUnit = dynamic_cast<t_gsetins*>(s)->MagUnit();
    _order = dynamic_cast<t_gsetins*>(s)->order();
}

int great::t_imufile::decode_head(char* buff, int sz, vector<string>& errmsg)
{
#ifdef BMUTEX
    boost::mutex::scoped_lock lock(_mutex);
#endif
    _mutex.lock();

    // no header expected, but fill the buffer
    t_gcoder::_add2buffer(buff, sz);
    _mutex.unlock(); return -1;
}

int great::t_imufile::decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
{
	_mutex.lock();

	if (t_gcoder::_add2buffer(buff, sz) == 0) { _mutex.unlock(); return 0; };

	int tmpsize = 0;
	string line;
	while (!_complete && (tmpsize = t_gcoder::_getline(line, 0)) >= 0)
	{
		if ((_order == "unicore" || _order == "great") && line[0] != '#')
		{
			t_gcoder::_consume(tmpsize);
			continue;
		}
		if ((line[0] == '%' && _order != "imr_fsas") || (line[0] == '#' && _order != "unicore" && _order != "great") || line[0] == 'T')
		{
			t_gcoder::_consume(tmpsize);
			continue;
		}
		for (int i = 0; i < line.size(); i++)
		{
			if (line[i] == ',' || line[i] == '*' || line[i] == ';')line[i] = ' ';
		}
		double t, gx, gy, gz, ax, ay, az, mx, my, mz;
		Vector3d g_tmp, a_tmp, m_tmp;
		Vector3d wtmp, vtmp, mtmp;
		stringstream ss(line);
		if (_order[0] == 'g' && _order[1] == 'a') // g->a
		{
			ss >> t >> gx >> gy >> gz >> ax >> ay >> az;
			g_tmp = Vector3d(gx, gy, gz); a_tmp = Vector3d(ax, ay, az);
		}
		else if (_order[0] == 'a' && _order[1] == 'g') // a->g
		{
			ss >> t >> ax >> ay >> az >> gx >> gy >> gz;
			g_tmp = Vector3d(gx, gy, gz); a_tmp = Vector3d(ax, ay, az);
		}
		else if (_order == "starneto")
		{
			vector<double> v;
			if (this->decode_starneto(line, t, v) > 0)
			{
				ax = v[3]; ay = v[4]; az = v[5];
				gx = v[0]; gy = v[1]; gz = v[2];
				g_tmp = Vector3d(gx, gy, gz); a_tmp = Vector3d(ax, ay, az);
				wtmp = g_tmp; vtmp = a_tmp;
			}
			else
			{
				if (_spdlog)
					SPDLOG_LOGGER_ERROR(_spdlog, "decode unicore data failed!");
				t_gcoder::_consume(tmpsize);
				continue;
			}
		}
		else
		{
			if (_spdlog) SPDLOG_LOGGER_ERROR(_spdlog, "The imu data format is not exist!");
		}

		if (_order.substr(_order.size() - 3, 3) == "rfu")  // right -> forward -> up
		{
			wtmp = g_tmp; vtmp = a_tmp;
		}
		else if (_order.substr(_order.size() - 3, 3) == "flu") // forward -> left -> up
		{
			Matrix3d Rotation; Rotation << 0, -1, 0, 1, 0, 0, 0, 0, 1;
			wtmp = Rotation * g_tmp; vtmp = Rotation * a_tmp;
		}
		else if (_order.substr(_order.size() - 3, 3) == "frd") // forward -> right -> down
		{
			Matrix3d Rotation; Rotation << 0, 1, 0, 1, 0, 0, 0, 0, -1;
			wtmp = Rotation * g_tmp; vtmp = Rotation * a_tmp;
		}
		else if (_order.substr(_order.size() - 3, 3) == "rbd") // right -> behind -> down
		{
			Matrix3d Rotation; Rotation << 1, 0, 0, 0, -1, 0, 0, 0, -1;
			wtmp = Rotation * g_tmp; vtmp = Rotation * a_tmp;
		}
		else if (_order.substr(_order.size() - 3, 3) == "lbu") // left -> behind -> up
		{
			Matrix3d Rotation; Rotation << -1, 0, 0, 0, -1, 0, 0, 0, 1;
			wtmp = Rotation * g_tmp; vtmp = Rotation * a_tmp;
		}
		else if (_order.substr(_order.size() - 3, 3) == "bru") // behind -> right -> up
		{
			Matrix3d Rotation; Rotation << 0, 1, 0, -1, 0, 0, 0, 0, 1;
			wtmp = Rotation * g_tmp; vtmp = Rotation * a_tmp;
		}
        else if (_order.substr(_order.size() - 3, 3) == "bld") // behind -> Left -> Down
        {
            Matrix3d Rotation;
            Rotation << 0, -1, 0, -1, 0, 0, 0, 0, -1;
            wtmp = Rotation * g_tmp; vtmp = Rotation * a_tmp;
        }
        else if (_order.substr(_order.size() - 3, 3) == "drb") // down -> right -> back
        {
            Matrix3d Rotation;
            Rotation << 0, 1, 0, 0, 0, -1, -1, 0, 0;
            wtmp = Rotation * g_tmp; vtmp = Rotation * a_tmp;
        }
        else if (_order.substr(_order.size() - 3, 3) == "dfr") // down -> forward -> right
        {
            Matrix3d Rotation;
            Rotation << 0, 0, 1, 0, 1, 0, -1, 0, 0;
            wtmp = Rotation * g_tmp; vtmp = Rotation * a_tmp;
        }
        else if (_order.substr(_order.size() - 3, 3) == "ufl") // Source: Up -> Forward -> Left
        {
            Matrix3d Rotation;
            Rotation << 0, 0, -1, 0, 1, 0, 1, 0, 0;
            wtmp = Rotation * g_tmp; vtmp = Rotation * a_tmp;
        }
        else if (_order.substr(_order.size() - 3, 3) == "urf") // Source: Up -> Right -> Forward
        {
            Matrix3d Rotation;
            Rotation << 0, 1, 0, 0, 0, 1, 1, 0, 0;
            wtmp = Rotation * g_tmp; vtmp = Rotation * a_tmp;
        }


		switch (_GyroUnit)
		{
		case RAD:
			break;
		case DEG:
			wtmp *= t_gglv::deg; break;
		case RPS:
			wtmp *= _ts; break;
		case DPS:
			wtmp *= t_gglv::dps * _ts; break;
		case RPH:
			wtmp /= t_gglv::hur * _ts; break;
		case DPH:
			wtmp *= t_gglv::dph * _ts; break;
		default:
			break;
		}
		switch (_AcceUnit)
		{
		case MPS:
			break;
		case MPS2:
			vtmp *= _ts; break;
		default:
			break;
		}

		if (t < dynamic_cast<t_gsetins*>(_set)->start())
		{
			t_gcoder::_consume(tmpsize);
			continue;
		}
		if (t > dynamic_cast<t_gsetins*>(_set)->end())
			_complete = true;

		_tt = t;

		map<string, t_gdata*>::iterator it = _data.begin();
		while (it != _data.end())
		{
			((t_gimudata*)it->second)->set_ts(_ts);
			if (it->second->id_type() == t_gdata::IMUDATA)
			{
				///modified by wh  read mag data
				if (_order[2] == 'm')
					((t_gimudata*)it->second)->add_IMU(t, wtmp, vtmp, mtmp);
				else
					((t_gimudata*)it->second)->add_IMU(t, wtmp, vtmp);
			}
			++it;
		}
		if (ss.fail())
		{
			if (_spdlog)
				SPDLOG_LOGGER_DEBUG(_spdlog, "imufile", "warning: incorrect IMU data record: " + ss.str());
			t_gcoder::_consume(tmpsize);
			_mutex.unlock(); return -1;
		}
		t_gcoder::_consume(tmpsize);
		cnt++;
	}
	_mutex.unlock();
	return 0;
}

int great::t_imufile::decode_starneto(const string & line, double & t, vector<double>& v)
{
    try
    {
        vector<string> ret;
        split(line, " ", ret);
        if (ret.size() != 11)  return -1;

        t = stod(ret[2]);

        for (int i = 3; i <= 5; i++)
        {
            v.push_back(stod(ret[i]) * t_gglv::deg * _ts);
        }
        for (int i = 6; i <= 8; i++)
            v.push_back(stod(ret[i]) * STARNETO_G * _ts);
        return 1;
    }
    catch (...)
    {
        return -1;
    }
}

int great::t_imufile::decode_starneto(const char* block, int sz, vector<pair<double, vector<double>>>& v)
{
    try
    {
        const int MARGIN = 50;
        int consume_size = 0;

        for (int i = 0; i < sz - MARGIN; i++)
        {
            if (block[i] == (char)0xAA && block[i + 1] == (char)0x55 && block[i + 2] == (char)0x05)
            {
                i += 3;
                int byten = 0;
                uint16_t    week;
                uint32_t    sow;
                double      gx;
                double      gy;
                double      gz;
                double      ax;
                double      ay;
                double      az;

                byten = sizeof(week); memcpy(&week, block + i, byten);  i += byten;
                byten = sizeof(sow);  memcpy(&sow, block + i, byten);  i += byten;
                byten = sizeof(gx);   memcpy(&gx, block + i, byten);  i += byten;
                byten = sizeof(gy);   memcpy(&gy, block + i, byten);  i += byten;
                byten = sizeof(gz);   memcpy(&gz, block + i, byten);  i += byten;
                byten = sizeof(ax);   memcpy(&ax, block + i, byten);  i += byten;
                byten = sizeof(ay);   memcpy(&ay, block + i, byten);  i += byten;
                byten = sizeof(az);   memcpy(&az, block + i, byten);  i += byten;

                gx *= t_gglv::deg * _ts;
                gy *= t_gglv::deg * _ts;
                gz *= t_gglv::deg * _ts;
                ax *= STARNETO_G * _ts;
                ay *= STARNETO_G * _ts;
                az *= STARNETO_G * _ts;

                v.push_back(make_pair(sow / 1000.0, vector<double>({ gx,gy,gz,ax,ay,az })));
            }
            else if (block[i] == (char)0xAA && block[i + 1] == (char)0x44 && block[i + 2] == (char)0x13)
            {
                int byten = 0;
                i += 6;

                uint16_t    week= 0;
                uint32_t    sow= 0;
                uint32_t    week2= 0;
                double      sow2= 0;
                int32_t     state= 0;
                int32_t     az= 0;
                int32_t     ay= 0;
                int32_t     ax= 0;
                int32_t     gz= 0;
                int32_t     gy= 0;
                int32_t     gx= 0;
                byten = sizeof(week);  memcpy(&week, block + i, byten); i += byten;
                byten = sizeof(sow);   memcpy(&sow, block + i, byten); i += byten;
                byten = sizeof(week2); memcpy(&week2, block + i, byten); i += byten;
                byten = sizeof(sow2);  memcpy(&sow2, block + i, byten); i += byten;
                byten = sizeof(state); memcpy(&state, block + i, byten); i += byten;
                byten = sizeof(az);    memcpy(&az, block + i, byten); i += byten;
                byten = sizeof(ay);    memcpy(&ay, block + i, byten); i += byten;
                byten = sizeof(ax);    memcpy(&ax, block + i, byten); i += byten;
                byten = sizeof(gz);    memcpy(&gz, block + i, byten); i += byten;
                byten = sizeof(gy);    memcpy(&gy, block + i, byten); i += byten;
                byten = sizeof(gx);    memcpy(&gx, block + i, byten); i += byten;

                double gyro[3] = { gx,-gy,gz };
                double accel[3] = { ax,-ay,az };

                for (int i = 0; i < 3; i++)  gyro[i] *= STARNETO_GYRO_SCALE;
                for (int i = 0; i < 3; i++) accel[i] *= STARNETO_ACC_SCALE;
                v.push_back(make_pair(sow / 1000.0, vector<double>({ gyro[0],gyro[1],gyro[2],accel[0],accel[1],accel[2] })));

            }
            consume_size = i;
        }

        return consume_size;
    }
    catch (...)
    {
        return -1;
    }
}

bool great::t_imufile::available(const t_gtime& now)
{
    double t = now.sow() + now.dsec();
    if (_tt > t)return true;
    else
    {
        return false;
    }
}

int great::t_imufile::encode_head(char* buff, int sz, vector<string>& errmsg)
{
#ifdef BMUTEX
    boost::mutex::scoped_lock lock(_mutex);
#endif
    _mutex.lock();

    _mutex.unlock(); return -1;
}
int great::t_imufile::encode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
{
#ifdef BMUTEX
    boost::mutex::scoped_lock lock(_mutex);
#endif
    _mutex.lock();
    try
    {
        if (_ss_position == 0)
        {
            t_gimudata* imudata = NULL;
            for (auto it = _data.begin(); it != _data.end(); ++it)
            {
                if (it->second->id_type() == t_gdata::IMUDATA)
                {
                    imudata = dynamic_cast<t_gimudata*>(it->second);
                }
            }
            if (imudata)
            {
                while (true)
                {
                    std::vector<Eigen::Vector3d> wm;
                    std::vector<Eigen::Vector3d> vm;
                    double t = -1.0, ts = -1.0;
                    imudata->load(wm, vm, t, ts, 1); // rfu
                    if (t <= 0)
                        break;
                    Eigen::Vector3d g_tmp = wm[0]; // RAD
                    Eigen::Vector3d a_tmp = vm[0]; // MPS

                    Eigen::Vector3d wtmp, vtmp;

                    if (_order.substr(2, 3) == "rfu")  // right -> forward -> up
                    {
                        wtmp = g_tmp; vtmp = a_tmp;
                    }
                    else if (_order.substr(2, 3) == "flu") // forward -> left -> up
                    {
                        Matrix3d Rotation; Rotation << 0, -1, 0, 1, 0, 0, 0, 0, 1;
                        wtmp = Rotation.transpose() * g_tmp; vtmp = Rotation.transpose() * a_tmp;
                    }
                    else if (_order.substr(2, 3) == "frd") // forward -> right -> down
                    {
                        Matrix3d Rotation; Rotation << 0, 1, 0, 1, 0, 0, 0, 0, -1;
                        wtmp = Rotation.transpose() * g_tmp; vtmp = Rotation.transpose() * a_tmp;
                    }
                    else if (_order.substr(2, 3) == "rbd") // right -> behind -> down
                    {
                        Matrix3d Rotation; Rotation << 1, 0, 0, 0, -1, 0, 0, 0, -1;
                        wtmp = Rotation.transpose() * g_tmp; vtmp = Rotation.transpose() * a_tmp;
                    }

                    switch (_GyroUnit)
                    {
                    case RAD:
                        break;
                    case DEG:
                        wtmp /= t_gglv::deg; break;
                    case RPS:
                        wtmp /= _ts; break;
                    case DPS:
                        wtmp /= t_gglv::dps * _ts; break;
                    case RPH:
                        wtmp *= t_gglv::hur * _ts; break;
                    case DPH:
                        wtmp /= t_gglv::dph * _ts; break;
                    default:
                        break;
                    }
                    switch (_AcceUnit)
                    {
                    case MPS:
                        break;
                    case MPS2:
                        vtmp /= _ts; break;
                    default:
                        break;
                    }

                    if (_order[0] == 'g' && _order[1] == 'a') // g->a
                    {
                        _ss << setw(15) << setiosflags(ios::fixed) << setprecision(4) << t << " " << resetiosflags(ios::fixed)
                            << setw(20) << setprecision(10) << wtmp(0) << " "
                            << setw(20) << setprecision(10) << wtmp(1) << " "
                            << setw(20) << setprecision(10) << wtmp(2) << " "
                            << setw(20) << setprecision(10) << vtmp(0) << " "
                            << setw(20) << setprecision(10) << vtmp(1) << " "
                            << setw(20) << setprecision(10) << vtmp(2) << " " << std::endl;
                    }
                    else if (_order[0] == 'a' && _order[1] == 'g') // a->g
                    {
                        _ss << setw(15) << setiosflags(ios::fixed) << setprecision(4) << t << " " << resetiosflags(ios::fixed)
                            << setw(20) << setprecision(10) << vtmp(0) << " "
                            << setw(20) << setprecision(10) << vtmp(1) << " "
                            << setw(20) << setprecision(10) << vtmp(2) << " "
                            << setw(20) << setprecision(10) << wtmp(0) << " "
                            << setw(20) << setprecision(10) << wtmp(1) << " "
                            << setw(20) << setprecision(10) << wtmp(2) << " " << std::endl;
                    }
                    else
                    {
                        if (_spdlog)
                            SPDLOG_LOGGER_ERROR(_spdlog, "ERROR : t_imufile::encode_data order bot supported.");
                        }
                        }
                    }
                }

            }
    catch (...)
    {
        if (_spdlog)
            SPDLOG_LOGGER_ERROR(_spdlog, "ERROR : t_imufile::encode_data throw exception");
        _mutex.unlock();
        return -1;
    }

    int size = _fill_buffer(buff, sz);
    _mutex.unlock();
    return size;
}



