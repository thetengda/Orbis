/**
 * @file         gsetsensors.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        control set from XML for MSF based on factor graph optimization
 * @version      1.0
 * @date         2025-11-04
 *
 * @copyright Copyright (c) 2024, Wuhan University. All rights reserved.
 *
 */

#ifndef GSETSENSORS_H
#define GSETSENSORS_H

#include "gset/gsetbase.h"
#include <Eigen/Eigen>
#include "gexport/ExportLibGREAT.h"

#define XMLKEY_GOMSF_SENSORS "fgomsf_sensors"
using namespace gnut;
namespace gfgomsf
{
	/** @brief class enum of MSF type. */
	enum MSF_TYPE 
	{
		INS_MODE,       ///< only INS mode
		GNSS_MODE,		///< only GNSS mode
		GINS_TC_MODE,   ///< GNSS/INS tightly coupled mode
		VIL_SLAM_MODE,	///< only VILO mode
		GVIL_LC_MODE,	///< GNSS/VILO Loosely combination mode
		GVIL_TC_MODE	///< GNSS/VILO tightly combination mode
	};

	/** @brief class enum of GNSS position type. */
	enum GNSS_POS_TYPE 
	{
		PPP,	///< PPP position
		RTK,	///< RTK position
		PPP_RTK	///< PPP-RTK position
	};

	/** @brief class enum of GNSS navigition frame. */
	enum NAV_REFERENCE_FRAME
	{
		E_F,    ///< ECEF frame
		L_F,    ///< local level frame (ENU frame)
		W_F,    ///< world frame  (slam used)
		B_F     ///< body frame   
	};

	/**
	* @brief change string to MSF_TYPE.
	* @param[in]	s			type of MSF in string form
	* @return		MSF_TYPE	type of MSF
	*/
	LibGREAT_LIBRARY_EXPORT  MSF_TYPE str2msf(const string &s);

	/**
	* @brief		class for set fgomsf_sensors xml
	*/
	class LibGREAT_LIBRARY_EXPORT t_gset_gomsf_sensors : public virtual t_gsetbase
	{
	public:	

		/**
		* @brief settings check.
		*/
		void check();

		/**
		* @brief settings help.
		*/
		void help();

		/**
		* @brief  get MSF type.
		* @return	MSF_TYPE	type of MSF
		*/
		MSF_TYPE msf_type();

		/**
		* @brief  get type of GNSS position.
		* @return	GNSS_POS_TYPE	type of GNSS position
		*/		

		string result_path();

	};
}
#endif