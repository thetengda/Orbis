/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 *
   (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)

   (c) 2011-2017 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  This file is part of the G-Nut C++ library.

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License as
  published by the Free Software Foundation; either version 3 of
  the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, see <http://www.gnu.org/licenses>.

-*/

#include <iomanip>
#include <sstream>

#include "gcfg_ppp.h"

using namespace std;
using namespace pugi;

namespace gnut
{

  // Constructor
  // ----------
  t_gcfg_ppp::t_gcfg_ppp()
      : t_gsetbase(),
        t_gsetgen(),
        t_gsetinp(),
        t_gsetout(),
        t_gsetgnss(),
        t_gsetproc(),
        t_gsetflt(),
        t_gsetrec(),
        t_gsetfgo(),
        t_gset_gomsf_sensors()
  {
    _IFMT_supported.insert(IFMT::RINEXO_INP);
    _IFMT_supported.insert(IFMT::RINEXC_INP);
    _IFMT_supported.insert(IFMT::RINEXN_INP);
    _IFMT_supported.insert(IFMT::ATX_INP);
    _IFMT_supported.insert(IFMT::BLQ_INP);
    _IFMT_supported.insert(IFMT::SP3_INP);
    _IFMT_supported.insert(IFMT::BIAS_INP);
    _IFMT_supported.insert(IFMT::BIASINEX_INP);
    _IFMT_supported.insert(IFMT::DE_INP);
    _IFMT_supported.insert(IFMT::EOP_INP);
    _IFMT_supported.insert(IFMT::LEAPSECOND_INP);
    _IFMT_supported.insert(IFMT::UPD_INP);
    _IFMT_supported.insert(IFMT::IFCB_INP);
    _OFMT_supported.insert(LOG_OUT);
    _OFMT_supported.insert(PPP_OUT);
    _OFMT_supported.insert(FLT_OUT);
    _OFMT_supported.insert(FGO_OUT);
    _OFMT_supported.insert(FGO_AR_OUT);
    _OFMT_supported.insert(RATIO_OUT);
  }

  // Destructor
  // ----------
  t_gcfg_ppp::~t_gcfg_ppp()
  {
  }

  // settings check (FGO-only settings are validated only when <est>=FGO so a
  // pure FLT configuration runs exactly as before, with no FGO defaults
  // injected into the parsed XML).
  // ----------
  void t_gcfg_ppp::check()
  {
    t_gsetgen::check();
    t_gsetinp::check();
    t_gsetout::check();
    t_gsetrec::check();
    t_gsetflt::check();
    t_gsetproc::check();
    t_gsetgnss::check();
    t_gsetamb::check();
    if (estimator() == "FGO")
    {
      t_gsetfgo::check();
      t_gset_gomsf_sensors::check();
    }
  }

  // settings help
  // ----------
  void t_gcfg_ppp::help()
  {
    t_gsetbase::help_header();
    t_gsetgen::help();
    t_gsetinp::help();
    t_gsetout::help();
    t_gsetrec::help();
    t_gsetflt::help();
    t_gsetproc::help();
    t_gsetgnss::help();
    t_gsetfgo::help();
    t_gset_gomsf_sensors::help();
    t_gsetbase::help_footer();
    t_gsetamb::help_footer();
  }

} // namespace