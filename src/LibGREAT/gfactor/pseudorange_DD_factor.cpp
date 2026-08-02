/**
 * @file         pseudorange_DD_factor.cpp
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        construct Double-Different pseudo-range factor
 * @version      1.0
 * @date         2025-11-04
 *
 * @copyright Copyright (c) 2025, Wuhan University. All rights reserved.
 *
 */

#include "pseudorange_DD_factor.h"


gfgo::PseudorangeDDFactor::PseudorangeDDFactor(const t_gtime & cur_time, const pair<string, string> &base_rover_site, const t_gallpar & params, const vector<pair<t_gsatdata, t_gsatdata>>& DD_sat_data, t_gprecisebiasFGO * bias_model, const pair<FREQ_SEQ, GOBSBAND> &freq_band):
	_cur_time(cur_time), _base_rover_site(base_rover_site),_params(params), _DD_sat_data(DD_sat_data), _gprecise_bias_model(bias_model), _freq_band(freq_band)
{
	
}

void gfgo::PseudorangeDDFactor::updatePara(t_gallpar & params_tmp, const Eigen::Vector3d & Pi, const Eigen::Vector3d & Vi) const
{
	params_tmp = _params;
	int i = 0;
	i = params_tmp.getParam(_base_rover_site.second, par_type::CRD_X, "");
	if (i >= 0)
	{
		params_tmp[i].value(Pi.x());
	}
	i = params_tmp.getParam(_base_rover_site.second, par_type::CRD_Y, "");
	if (i >= 0)
	{
		params_tmp[i].value(Pi.y());
	}

	i = params_tmp.getParam(_base_rover_site.second, par_type::CRD_Z, "");
	if (i >= 0)
	{
		params_tmp[i].value(Pi.z());
	}
}

void gfgo::PseudorangeDDFactor::trans2Eigen(const vector<vector<pair<int, double>>>& B, const vector<double>& P, const vector<double>& l, Eigen::Matrix<double, 2, 3>& B_new, Eigen::Matrix<double, 2, 2>& P_new, Eigen::Matrix<double, 2, 1> & l_new) const
{
	B_new.setZero();
	P_new.setZero();
	l_new.setZero();
	for (int i = 0; i < B.size(); i++)
	{
		for (int j = 0; j < B[i].size(); j++)
		{
			B_new(i, j) = B[i][j].second;
		}
	}
	for (int i = 0; i < 2; i++)
	{
		P_new(i, i) = P[i];
	}
	
	for (int i = 0; i < 2; i++)
	{
		l_new(i) = l[i];			
	}
}

bool gfgo::PseudorangeDDFactor::Evaluate(double const * const * parameters, double * residuals, double ** jacobians) const
{
	Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
	//Eigen::Vector3d Vi(parameters[1][0], parameters[1][1], parameters[1][2]); 
	double sqrt_info;
	Eigen::Matrix<double, 1, 2> DD_operator(1, 2);
	//construct DD equ	
	unsigned npar_orig = _params.parNumber()-5;
	t_gallpar params_temp;	
	updatePara(params_temp, Pi);	
	t_gtriple xyz;
	params_temp.getCrdParam(_base_rover_site.second, xyz);
	Eigen::Vector3d xyz_before = Eigen::Vector3d(xyz.crd(0), xyz.crd(1), xyz.crd(2));
	//cout << " range YXZ update: " << xyz_before.transpose() << endl;

	vector<vector<pair<int, double>>> B;        ///< coeff of equations
	vector<double> P;                           ///< weight of equations
	vector<double> l;                           ///< res of equations
	Eigen::MatrixXd B_DD;
	double l_DD, P_DD;
	for (auto it : _DD_sat_data)
	{
		t_gbaseEquation tempP;
		pair<t_gsatdata, t_gsatdata> rec_pair = it;
		for (int isite = 0; isite < 2; isite++)
		{
			t_gsatdata *satdata_ptr;
			if (isite == 0) satdata_ptr = &rec_pair.first;
			else satdata_ptr = &rec_pair.second;
			t_gobs  obsP = t_gobs(satdata_ptr->select_range(_freq_band.second));
			t_gtime crt = satdata_ptr->epoch();			
			if (!_gprecise_bias_model->cmb_equ(true,true,crt, params_temp, *satdata_ptr, obsP, tempP))
			{
				cout << "sat " << rec_pair.second.sat() << "  construct pseudorange DD factor error" << endl;
				return false;

			}
		}
		vector<pair<int, double>> B_P;
		double P_P, l_P;
		int ibase = 0;
		int irover = 1;
		for (const auto& b : tempP.B[irover]) {
			if (b.first > npar_orig) continue;
			B_P.push_back(b);
		}
		for (const auto& b : tempP.B[ibase]) {
			if (b.first > npar_orig) continue;
			B_P.emplace_back(b.first, -b.second);
		}
		P_P = 1 / (1 / tempP.P[irover] + 1 / tempP.P[ibase]); l_P = tempP.l[irover] - tempP.l[ibase];	

		B.push_back(B_P);
		P.push_back(P_P);
		l.push_back(l_P);


	}
	int iobs = 1;
	int index_ref = 0;
	int index_sat = 1;
	DD_operator(iobs-1, index_ref) = -1;
	DD_operator(iobs-1, index_sat) = 1;
	Eigen::Matrix<double, 2, 3> B_new;
	Eigen::Matrix<double, 2, 2> P_new;
	Eigen::Matrix<double, 2, 1> l_new;
	trans2Eigen(B,P,l,B_new,P_new,l_new);
	B_DD = DD_operator * B_new;
	l_DD = DD_operator * l_new;
	P_DD = DD_operator * P_new.inverse()*DD_operator.transpose();
	P_DD = 1.0 / P_DD;
	//cout << "residual: " << l_DD << endl;
	//cout << "Jacbian: " << B_DD << endl;
	//cout << "weight: " << P_DD << endl;
	//cout << endl;
	//set ceres value
	sqrt_info = sqrt(P_DD);
	residuals[0]= sqrt_info * l_DD;
	//residuals[0] = l_DD;
	//cout << "range residual: " << residuals[0] << endl;
	if (jacobians)
	{
		if (jacobians[0])
		{
			Eigen::Map<Eigen::Matrix<double, 1, 3, Eigen::RowMajor>> jacobian_XYZ(jacobians[0]);
			jacobian_XYZ = -sqrt_info*B_DD;
			//cout << " range pos jacobian: " << jacobian_XYZ.transpose() << endl;
			//jacobian_XYZ = B_DD;
		}
		/*if (jacobians[1])
		{   TODO：for velocity parameter block;


		}*/

	}
	return true;

	//////////if (!jacobians) return true;

	//////////bool need_recompute = !_jac_cached;
	//////////if (_jac_cached) {
	//////////	const double dp = (Pi - _lin_P).norm();
	//////////	if (dp > 0.02) need_recompute = true;
	//////////}

	//////////// 缓存可用：直接回填
	//////////if (_freeze_jac && _jac_cached && !need_recompute)
	//////////{
	//////////	if (jacobians[0])
	//////////	{
	//////////		Eigen::Map<Eigen::Matrix<double, 1, 3, Eigen::RowMajor>> jacobian_XYZ(jacobians[0]);
	//////////		jacobian_XYZ(0, 0) = _J_p[0];
	//////////		jacobian_XYZ(0, 1) = _J_p[1];
	//////////		jacobian_XYZ(0, 2) = _J_p[2];
	//////////	}
	//////////	return true;
	//////////}
	//////////else
	//////////{
	//////////	// 需要重算：按你原方式算 Jacobian（Map 写法保持）
	//////////	if (jacobians[0])
	//////////	{
	//////////		Eigen::Map<Eigen::Matrix<double, 1, 3, Eigen::RowMajor>> jacobian_XYZ(jacobians[0]);
	//////////		jacobian_XYZ = -sqrt_info * B_DD;  // 你原来就是这句（或等价）
	//////////	}
	//////////}
	//////////// 写缓存（缓存的是写入 ceres 的 whitened Jacobian）
	//////////if (_freeze_jac)
	//////////{
	//////////	if (jacobians[0]) {
	//////////		_J_p[0] = jacobians[0][0];
	//////////		_J_p[1] = jacobians[0][1];
	//////////		_J_p[2] = jacobians[0][2];
	//////////	}
	//////////	else {
	//////////		_J_p[0] = _J_p[1] = _J_p[2] = 0.0;
	//////////	}
	//////////	_lin_P = Pi;
	//////////	_jac_cached = true;
	//////////}

	//////////return true;

}

void gfgo::PseudorangeDDFactor::check(double ** parameters)
{
	double *res = new double[1];
	double **jaco = new double *[2];
	jaco[0] = new double[1 * 3];	
	Evaluate(parameters, res, jaco);
	/*puts("check begins");

	puts("my");*/

	/*std::cout << Eigen::Map<Eigen::Matrix<double, 1, 1>>(res).transpose() << std::endl
		<< std::endl;
	std::cout << Eigen::Map<Eigen::Matrix<double, 1, 3, Eigen::RowMajor>>(jaco[0]) << std::endl
		<< std::endl;*/

	

}
