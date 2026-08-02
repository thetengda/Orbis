
#ifndef GCFG_PPP_H
#define GCFG_PPP_H

/**
 * @file         gcfg_ppp.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        control set from XML for main
 * @version      1.0
 * @date         2025-11-04
 *
 * @copyright Copyright (c) 2025, Wuhan University. All rights reserved.
 *
 */

#include <string>
#include <iostream>
#include <signal.h>
#include "gfgognss/gpvtfgo.h"
#include "gall/gallprec.h"
#include "gall/gallpcv.h"
#include "gall/gallobs.h"
#include "gall/gallotl.h"
#include "gall/gallbias.h"
#include "gproc/gpppflt.h"
#include "gproc/gpreproc.h"
#include "gio/gfile.h"
#include "gutils/gtime.h"
#include "gutils/gtypeconv.h"
#include "gdata/gdata.h"
#include "gcoders/gcoder.h"
#include "gcoders/rinexo.h"
#include "gcoders/rinexc.h"
#include "gcoders/rinexn.h"
#include "gcoders/biasinex.h"
#include "gcoders/biabernese.h"
#include "gcoders/sp3.h"
#include "gcoders/atx.h"
#include "gcoders/blq.h"
#include "gmodels/gpcv.h"
#include "gmodels/gbancroft.h"

#include "gset/gsetgen.h"
#include "gset/gsetinp.h"
#include "gset/gsetout.h"
#include "gset/gsetproc.h"
#include "gset/gsetgnss.h"
#include "gset/gsetflt.h"
#include "gset/gsetrec.h"

#include "gproc/gpvtflt.h"
#include "gio/gfile.h"
#include "gdata/gifcb.h"
#include "gcoders/upd.h"

#include "gcoders/sp3.h"
#include "gcoders/rinexo.h"
#include "gcoders/rinexn.h"
#include "gcoders/rinexc.h"
#include "gcoders/poleut1.h"
#include "gcoders/ifcb.h"
#include "gcoders/gcoder.h"
#include "gcoders/dvpteph405.h"
#include "gcoders/blq.h"
#include "gcoders/biasinex.h"
#include "gcoders/biabernese.h"
#include "gcoders/atx.h"
#include "gdata/gnavde.h"
#include "gset/gsetfgo.h"
#include "gset/gsetsensors.h"

using namespace std;
using namespace pugi;
using namespace gfgo;
using namespace gfgomsf;

namespace gnut
{

  class t_gcfg_ppp : public t_gsetgen,
                     public t_gsetinp,
                     public t_gsetout,
                     public t_gsetgnss,
                     public t_gsetproc,
                     public t_gsetflt,
                     public t_gsetrec,
					 public t_gsetamb,
					 public t_gsetfgo,
	                 public t_gset_gomsf_sensors
  {

  public:
    t_gcfg_ppp();
    ~t_gcfg_ppp();

    void check(); // settings check
    void help();  // settings help

  protected:
  private:
  };

} // namespace

#endif