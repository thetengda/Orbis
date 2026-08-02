#ifndef P_IF_FACTOR
#define P_IF_FACTOR
/**
* @file pseudorange_IF_factor.h
* @details
* @note  construct IF pseudorange factor
* @verbatim
		History
		-0.1    hyChang        2022-05-12 creat the file.

  @endverbatim

* @author hyChang
* @version 0.1
* @date 2022-05-12
* @license
*/
#include"gfgo/gutility.h"
#include "gexport/ExportLibGREAT.h"
#include "gfgo/gprecisebiasFGO.h"
using namespace great;
namespace gfgo
{
	class LibGREAT_LIBRARY_EXPORT PseudorangeIFFactor : public ceres::SizedCostFunction<1, 3, 1, 1>//res, CRD, CLK, TRP (for GPS only)
	{
	public:
		PseudorangeIFFactor(const t_gtime &cur_time, const string &site, const t_gallpar  &params, const t_gsatdata &IF_sat_data, t_gprecisebiasFGO *bias_model, const pair<FREQ_SEQ, GOBSBAND> &freq_band1, const pair<FREQ_SEQ, GOBSBAND> &freq_band2);
		void updatePara(t_gallpar & params_tmp, const Eigen::Vector3d & Pi, const double& Clk, const double&Trp) const;
		virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
	protected:
		t_gtime _cur_time;
		string _site;
		t_gallpar _params;
		t_gprecisebiasFGO *_gprecise_bias_model = nullptr;
		pair<FREQ_SEQ, GOBSBAND> _freq_band1;
		pair<FREQ_SEQ, GOBSBAND> _freq_band2;
		t_gsatdata _IF_sat_data;
	};
}

#endif
