#ifndef ML_IF_ING_FACTOR
#define ML_IF_ING_FACTOR
/**
* @file multi_carrierphase_IF_integration_factor.h
* @details
* @note  construct IF pseudorange factor for GNSS/INS integration (GAL & BDS)
* @verbatim
		History
		-0.1    hyChang        2022-06-10 creat the file.

  @endverbatim

* @author hyChang
* @version 0.1
* @date 2022-06-10
* @license
*/
#include"gfgo/gutility.h"
#include "gfgo/gprecisebiasFGO.h"
using namespace great;
namespace gfgo
{
	class LibGREAT_LIBRARY_EXPORT MultiCarrierphaseIFINGFactor : public ceres::SizedCostFunction<1, 7, 1, 1, 1, 1>//res, CRD, CLK, TRP, ISB, AMB_IF (ISB for GAL & BDS)
	{
	public:
		MultiCarrierphaseIFINGFactor(const t_gtime &cur_time, const string &site, const t_gallpar  &params, const t_gsatdata &IF_sat_data, t_gprecisebiasFGO *bias_model, const pair<FREQ_SEQ, GOBSBAND> &freq_band1, const pair<FREQ_SEQ, GOBSBAND> &freq_band2, const Eigen::Vector3d &lever_arm);
		void updatePara(t_gallpar & params_tmp, const Eigen::Vector3d & Pi, const double& Clk, const double&Trp, const double&ISB, const double&Amb) const;
		virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
	protected:
		t_gtime _cur_time;
		string _site;
		t_gallpar _params;
		t_gprecisebiasFGO *_gprecise_bias_model = nullptr;
		pair<FREQ_SEQ, GOBSBAND> _freq_band1;
		pair<FREQ_SEQ, GOBSBAND> _freq_band2;
		t_gsatdata _IF_sat_data;
		Eigen::Vector3d _lever_arm;
	};
}

#endif