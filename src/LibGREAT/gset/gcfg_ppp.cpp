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

namespace great
{

// Constructor
// ----------
t_gcfg_ppp::t_gcfg_ppp() 
 : t_gsetgen(),
   t_gsetinp(),
   t_gsetout(),
   t_gsetgnss(),
   t_gsetproc(),
   t_gsetflt(),
   t_gsetrec(),
   t_gsetamb()
{
  _IFMT_supported.insert(IFMT::RINEXO_INP);
  _IFMT_supported.insert(IFMT::RINEXC_INP);
  _IFMT_supported.insert(IFMT::RINEXN_INP);
  _IFMT_supported.insert(IFMT::ATX_INP);
  _IFMT_supported.insert(IFMT::BLQ_INP);
  _IFMT_supported.insert(IFMT::SP3_INP);
  _IFMT_supported.insert(IFMT::BIAS_INP);
  _IFMT_supported.insert(IFMT::BIASINEX_INP);
  _IFMT_supported.insert(IFMT::UPD_INP);
  _IFMT_supported.insert(IFMT::EOP_INP);
  _IFMT_supported.insert(IFMT::LEAPSECOND_INP);
  _IFMT_supported.insert(IFMT::DE_INP);
  _IFMT_supported.insert(IFMT::BIASINEX_INP);
  _IFMT_supported.insert(IFMT::IFCB_INP);
  _OFMT_supported.insert(LOG_OUT);
  _OFMT_supported.insert(PPP_OUT);
  _OFMT_supported.insert(FLT_OUT);
  _OFMT_supported.insert(KML_OUT);
}


// Destructor
// ----------
t_gcfg_ppp::~t_gcfg_ppp()
    {
    }

// settings check
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
  t_gsetbase::help_footer();
  t_gsetamb::help_footer();
}

} // namespace
