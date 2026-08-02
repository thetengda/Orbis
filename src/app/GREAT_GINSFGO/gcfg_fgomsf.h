#ifndef GCFG_FGOMSF_H
#define GCFG_FGOMSF_H

/**
 * @file         gcfg_fgomsf.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        control set from XML for main
 * @version      1.0
 * @date         2025-11-04
 *
 * @copyright Copyright (c) 2025, Wuhan University. All rights reserved.
 *
 */

#include "gcoders/imufile.h"
#include "gset/gsetins.h"
#include "gset/gcfg_ppp.h"
#include "gset/gsetfgo.h"
#include "gset/gsetign.h"
#include "gset/gsetsensors.h"
#include "gio/gfile.h"
#include "gfgomsf/gfgo_gins.h"
#include "gdata/gifcb.h"
#include "gcoders/ifcb.h"

using namespace std;
using namespace gnut;
using namespace great;
using namespace gfgomsf;
using namespace gfgo;

class t_gcfg_fgomsf : 
	public t_gcfg_ppp,	
	public t_gsetins,	
	public t_gsetign, 
    public t_gsetfgo,
	public t_gset_gomsf_sensors
{

public:
	t_gcfg_fgomsf();
	~t_gcfg_fgomsf();

	void check();
	void help();


};



#endif