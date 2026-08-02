#ifndef GPARA_H
#define GPARA_H
/**
 * @file         gpara.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        Definition of  common parameters.
 * @version      1.0
 * @date         2025-11-04
 *
 * @copyright Copyright (c) 2025, Wuhan University. All rights reserved.
 *
 */

#include"gfgo/gfgo_para.h"
#include "gset/gsetins.h"
#include"gset/gsetbase.h"
#include <Eigen/Dense>
#include "gexport/ExportLibGREAT.h"
#include <gset/gsetign.h>
#include "gset/gsetsensors.h"


using namespace gfgo;
using namespace great;
namespace gfgomsf
{
	/**
	*@brief t_gara Class for storing common parameters.
	*
	*/
	class LibGREAT_LIBRARY_EXPORT t_gpara: virtual public t_gfgo_para
	{
	public:

		/**
		 * @brief Navigation parameter initialization constructor
		 * Sets up frame type, integration mode, IMU specs and output paths
		 */
		explicit t_gpara(gnut::t_gsetbase* gset);			/// get Parameter from xml
		~t_gpara();
		
	protected:
		NAV_REFERENCE_FRAME _nav_frame;             ///< navigation frame used         
		MSF_TYPE _msf_type;							///< type of mutil_slam_fusion
		double _imu_ts=0.0;							///< imu sample interval
		string _result_path;                        ///< file path for outputing results
	};

}
#endif
