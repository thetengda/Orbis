#ifndef GPRECISEBIASFGO_H
#define GPRECISEBIASFGO_H
/**
*
* @verbatim
History
-1.0  2022-03-24  hyChang         create the file

@endverbatim
* Copyright (c) 2020, Wuhan University. All rights reserved.
*
* @file            gprecisebiasFGO.h
* @brief           the class for FGO precise GNSS procrssing
*
* @author       GREAT Wuhan University
* @version      1.0.0
* @date         2022-03-24
*
*/

#include "gexport/ExportLibGREAT.h"
#include "gmodels/gprecisebiasGPP.h"

namespace gfgo
{
	/**
	*@brief       Class for t_gprecisebiasFGO, derive from t_gprecisebiasGPP
	*/
	class LibGREAT_LIBRARY_EXPORT t_gprecisebiasFGO : public t_gprecisebiasGPP
	{
	public:
		/**
		* @brief constructor.
		*
		* @param[in]  data             the data pointer
		* @param[in]  spdlog           the log
		* @param[in]  setting          the set pointer
		*/
		explicit t_gprecisebiasFGO(t_gallproc *data, t_spdlog spdlog, t_gsetbase *setting);

		/** @brief default destructor. */
		~t_gprecisebiasFGO();

		/**
		* @brief combine EQU.
		*
		* @param[in]  isFGO              for FGO or EKF/LSQ
		* @param[in]  calculate_equ      calculate t_gbaseEquation or not
		* @param[in]  epoch              the current time
		* @param[in]  params          the parameters
		* @param[in]  obsdata          the observation data
		* @param[in]  gobs              the observation object
		* @param[in]  result          the EQU result
		* @return      bool              combine EQU mode
		*/
		bool cmb_equ(bool isFGO, bool calculate_equ, t_gtime &epoch, t_gallpar &params, t_gsatdata &obsdata, t_gobs &gobs, t_gbaseEquation &result);

		/**
		* @brief prepare observation for GPP.(FGO)
		*
		* @param[in]  epoch              the current time
		* @param[in]  nav             the navigation data
		* @param[in]  gallobj          the all object
		* @param[in]  pars              the pars
		* @return      bool              prepare mode
		*/
		bool _prepare_obs_GPP_FGO(const t_gtime &crt_epo, t_gallnav *gallnav, t_gallobj *gallobj, t_gallpar &pars);

		/**
		* @brief prt all observation.
		*
		* @param[in]  crt_epo          the current time
		* @param[in]  obsdata          the satellite data
		* @param[in]  pars             the paramter
		* @param[in]  gobs            the obs
		* @param[in]  coeff           the coeff
		* @return      bool             the prt mode
		*/
		bool _prt_obs_ALL(const t_gtime &crt_epo, t_gsatdata &obsdata, t_gallpar &pars, t_gobs &gobs, vector<pair<int, double>> &coeff);

	};
}


#endif