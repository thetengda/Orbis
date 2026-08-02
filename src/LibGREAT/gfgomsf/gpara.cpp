/**
 * @file         gpara.cpp
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        Definition of  common parameters.
 * @version      1.0
 * @date         2025-11-04
 *
 * @copyright Copyright (c) 2025, Wuhan University. All rights reserved.
 *
 */

#include "gpara.h"

gfgomsf::t_gpara::t_gpara(gnut::t_gsetbase * gset):
t_gfgo_para(gset)
{
	_nav_frame = E_F;
	_msf_type = dynamic_cast<t_gset_gomsf_sensors*>(gset)->msf_type();
	_imu_ts = dynamic_cast<t_gsetins*>(gset)->ts();
	_result_path = dynamic_cast<t_gset_gomsf_sensors*>(gset)->result_path();
}

gfgomsf::t_gpara::~t_gpara()
{


}
