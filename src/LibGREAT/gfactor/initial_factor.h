#ifndef INITIAL_FACTOR
#define INITIAL_FACTOR
/**
* @file initial_factor.h
* @details
* @note  construct initial factor
* @verbatim
		History
		-0.1    hyChang        2022-05-23 creat the file.

  @endverbatim

* @author hyChang
* @version 0.1
* @date 2022-05-22
* @license
*/
#include"gfgo/gutility.h"
#include "gexport/ExportLibGREAT.h"
namespace gfgo
{
	class LibGREAT_LIBRARY_EXPORT InitialFactor : public ceres::SizedCostFunction<1, 1>//res, X
	{
	public:
		InitialFactor(const double & value, const double & sqrt_info);
		virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
	protected:
		double _value;
		double _sqrt_info;
	};
}

#endif