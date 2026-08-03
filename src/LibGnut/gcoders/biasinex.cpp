/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 *
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  This file is part of the G-Nut C++ library.
-*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <cstring>
#include <memory>

#include "gcoders/biasinex.h"
#include "gutils/gtypeconv.h"
#include "gmodels/gbias.h"
#include "gall/gallbias.h"

using namespace std;

namespace gnut
{

    t_biasinex::t_biasinex(t_gsetbase *s, string version, int sz, string id)
        : t_sinex(s, version, sz, id)
    {

        _allbias = 0;
    }

    int t_biasinex::_decode_comm()
    {
        std::string::size_type idx;
        if ((idx = _line.find("PRN")) != string::npos)
            _mapidx["SAT"] = make_pair(idx, 3);
        if ((idx = _line.find("OBS1")) != string::npos)
            _mapidx["OBS1"] = make_pair(idx, 4);
        if ((idx = _line.find("OBS2")) != string::npos)
            _mapidx["OBS2"] = make_pair(idx, 4);

        if ((idx = _line.find("BIAS_START____")) != string::npos)
            _mapidx["BEG"] = make_pair(idx, 14);
        if ((idx = _line.find("BIAS_END______")) != string::npos)
            _mapidx["END"] = make_pair(idx, 14);

        if ((idx = _line.find("__ESTIMATED_VALUE____")) != string::npos)
            _mapidx["EST"] = make_pair(idx, 21);
        if ((idx = _line.find("_STD_DEV___")) != string::npos)
            _mapidx["STD"] = make_pair(idx, 11);

        if ((idx = _line.find("__ESTIMATED_SLOPE____")) != string::npos)
            _mapidx["SLOPE"] = make_pair(idx, 21);

        return 1;
    }

    int t_biasinex::_decode_block()
    {

        t_sinex::_decode_block();

        // -------- "BIAS/DESCRIPTION" --------
        if (_block.find("BIAS/DESCRIPTION") != string::npos)
        {

            if (_line.find(" OBSERVATION_SAMPLING ") != string::npos)
            {
                if (_spdlog)
                    SPDLOG_LOGGER_DEBUG(_spdlog, "Read BIAS/SMP: " + cut_crlf(_line.substr(31)));
            }
            else if (_line.find(" PARAMETER_SPACING ") != string::npos)
            {
                if (_spdlog)
                    SPDLOG_LOGGER_DEBUG(_spdlog, "Read BIAS/SPC: " + cut_crlf(_line.substr(31)));
            }
            else if (_line.find(" DETERMINATION_METHOD ") != string::npos)
            {
                if (_spdlog)
                    SPDLOG_LOGGER_DEBUG(_spdlog, "Read BIAS/MTD: " + cut_crlf(_line.substr(31)));
            }
            else if (_line.find(" BIAS_MODE ") != string::npos)
            {
                if (_spdlog)
                    SPDLOG_LOGGER_DEBUG(_spdlog, "Read BIAS/MOD: " + cut_crlf(_line.substr(31)));
            }
            else if (_line.find(" TIME_SYSTEM ") != string::npos)
            {
                string ts_str = trim(cut_crlf(_line.substr(31)));
                if (ts_str.size() == 1)
                {
                    switch (ts_str[0])
                    {
                    case 'G': _bias_tsys = t_gtime::GPS; break;
                    case 'C': _bias_tsys = t_gtime::BDS; break;
                    case 'R': _bias_tsys = t_gtime::GLO; break;
                    case 'E': _bias_tsys = t_gtime::GAL; break;
                    default:  _bias_tsys = t_gtime::TAI; break;
                    }
                }
                else
                {
                    _bias_tsys = t_gtime::str2tsys(ts_str);
                }
                if (_spdlog)
                    SPDLOG_LOGGER_DEBUG(_spdlog, "Read BIAS/TSY: " + ts_str);
            }
        }
        else if (_block.find("BIAS/SOLUTION") != string::npos)
        {

            string prn = "";
            GOBS gobs1, gobs2;
            gobs1 = gobs2 = X;
            t_gtime beg = FIRST_TIME;
            t_gtime end = LAST_TIME;
            beg.tsys(_bias_tsys);
            end.tsys(_bias_tsys);
            double dcb = 0.0;
            double slope_ns = 0.0;

            for (auto it = _mapidx.begin(); it != _mapidx.end(); it++)
            {
                size_t pos = it->second.first;
                size_t len = it->second.second;
                if (it->first == "SAT")
                    prn = _line.substr(pos, len);
                if (it->first == "OBS1")
                    gobs1 = str2gobs(_line.substr(pos, len));
                if (it->first == "OBS2")
                    gobs2 = str2gobs(_line.substr(pos, len));
                if (it->first == "BEG")
                    beg.from_str("%Y:%j:%s", _line.substr(pos, len));
                if (it->first == "END" && _line.substr(pos, len) != "0000:000:00000")
                    end.from_str("%Y:%j:%s", _line.substr(pos, len));
                if (it->first == "EST")
                    dcb = str2dbl(_line.substr(pos, len));
                if (it->first == "SLOPE")
                    slope_ns = str2dbl(_line.substr(pos, len));
            }
            shared_ptr<t_gbias> p_bias;

            if (_allbias)
            {
                p_bias = make_shared<t_gbias>(_spdlog);
                double val = dcb * 1e-9 * CLIGHT;
                double slope = slope_ns * 1e-9 * CLIGHT;
                p_bias->set(beg, end, val, slope, gobs1, gobs2);
                _allbias->add(_ac, beg, prn, p_bias);
            }
        }

        return 1;
    }

    void t_biasinex::_add_data(string id, t_gdata *pt_data)
    {

        if (pt_data == 0)
            return;

        // ALL OBJECTS
        if (pt_data->id_type() == t_gdata::ALLBIAS)
        {
            if (!_allbias)
            {
                _allbias = dynamic_cast<t_gallbias *>(pt_data);
            }
        }

        return;
    }

} // namespace
