/**
 * @file         carrierphase_DD_factor.cpp
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        Class for factor graph-based GNSS/SINS integrated solution
 * @version      1.0
 * @date         2025-11-04
 *
 * @copyright Copyright (c) 2025, Wuhan University. All rights reserved.
 *
 */
#include "carrierphase_DD_factor.h"

gfgo::CarrierphaseDDFactor::CarrierphaseDDFactor(const t_gtime & cur_time, const pair<string, string> &base_rover_site, const t_gallpar & params, const vector<pair<t_gsatdata, t_gsatdata>>& DD_sat_data, t_gprecisebiasFGO *bias_model, const pair<FREQ_SEQ, GOBSBAND> &freq_band):
	_cur_time(cur_time),_base_rover_site(base_rover_site), _params(params), _DD_sat_data(DD_sat_data), _gprecise_bias_model(bias_model), _freq_band(freq_band)
{
	
}

void gfgo::CarrierphaseDDFactor::updatePara(t_gallpar & params_tmp, const double &ref_sd_amb, const double &nonref_sd_amb, const Eigen::Vector3d & Pi, const Eigen::Vector3d & Vi) const
{
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
	map<FREQ_SEQ, par_type> ambtype_list = {
				{FREQ_1, par_type::AMB_L1},
				{FREQ_2, par_type::AMB_L2},
				{FREQ_3, par_type::AMB_L3},
				{FREQ_4, par_type::AMB_L4},
				{FREQ_5, par_type::AMB_L5} };
	//for ref sat
	string ref_sat_name = _DD_sat_data[0].first.sat();
	string site = _DD_sat_data[0].second.site(); //for rover

	i = params_tmp.getParam(site, ambtype_list[_freq_band.first], ref_sat_name);
	if (i >= 0)
	{
		params_tmp[i].value(ref_sd_amb);
	}
	//for nonref sat
	string nonref_sat_name = _DD_sat_data[1].first.sat();
	i = params_tmp.getParam(site, ambtype_list[_freq_band.first], nonref_sat_name);
	if (i >= 0)
	{
		params_tmp[i].value(nonref_sd_amb);
	}
}

void gfgo::CarrierphaseDDFactor::trans2Eigen(const vector<vector<pair<int, double>>>& B, const vector<double>& P, const vector<double>& l, Eigen::Matrix<double, 2, 5>& B_new, Eigen::Matrix<double, 2, 2>& P_new, Eigen::Matrix<double, 2, 1>& l_new) const
{
	B_new.setZero();
	P_new.setZero();
	l_new.setZero();
	for (int i = 0; i < B.size(); i++)
	{
		for (int j = 0; j < 3; j++)
		{
			B_new(i, j) = B[i][j].second;
		}
		
		if (i == 0)  B_new(i, 3) = B[i][3].second;
		if (i == 1)  B_new(i, 4) = B[i][3].second;
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

bool gfgo::CarrierphaseDDFactor::Evaluate(double const * const * parameters, double * residuals, double ** jacobians) const
{	
	Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);	
	//Eigen::Vector3d Vi(parameters[1][0], parameters[1][1], parameters[1][2]);
	double ref_SD_ambiguity = parameters[1][0];
	double nonref_SD_ambiguity = parameters[2][0];	
	double sqrt_info;
	Eigen::Matrix<double, 1, 2> DD_operator(1, 2);
	//construct DD equ	
	unsigned npar_orig = _params.parNumber()-5;
	t_gallpar params_temp = _params;
	t_gtriple xyz;
	updatePara(params_temp,ref_SD_ambiguity, nonref_SD_ambiguity, Pi);
	//params_temp.getCrdParam(_base_rover_site.second, xyz);
	//Eigen::Vector3d xyz_before = Eigen::Vector3d(xyz.crd(0), xyz.crd(1), xyz.crd(2));
    //cout << " phase YXZ update: " << xyz_before.transpose() << endl;
	vector<vector<pair<int, double>>> B;        ///< coeff of equations
	vector<double> P;                           ///< weight of equations
	vector<double> l;                           ///< res of equations
	Eigen::MatrixXd B_DD;
	double l_DD, P_DD;
	map<FREQ_SEQ, par_type> ambtype_list = {
				{FREQ_1, par_type::AMB_L1},
				{FREQ_2, par_type::AMB_L2},
				{FREQ_3, par_type::AMB_L3},
				{FREQ_4, par_type::AMB_L4},
				{FREQ_5, par_type::AMB_L5} };

	for (auto it : _DD_sat_data)
	{
		t_gbaseEquation tempL;
		pair<t_gsatdata, t_gsatdata> rec_pair = it;
		for (int isite = 0; isite < 2; isite++)
		{
			t_gsatdata *satdata_ptr;
			if (isite == 0) satdata_ptr = &rec_pair.first;
			else satdata_ptr = &rec_pair.second;
			t_gobs  obsL = t_gobs(satdata_ptr->select_phase(_freq_band.second));
			t_gtime crt = satdata_ptr->epoch();
			if (!_gprecise_bias_model->cmb_equ(true,true,crt, params_temp, *satdata_ptr, obsL, tempL))
			{
				cout << "sat " << rec_pair.second.sat() << "  construct carrierphase DD factor error" << endl;
				return false;
			}
			if (satdata_ptr->site() == _base_rover_site.second)
			{
				int idx = params_temp.getParam(satdata_ptr->site(), ambtype_list[_freq_band.first], satdata_ptr->sat());

				if (idx < 0)
				{
					cout << "sat " << rec_pair.second.sat() << "  construct carrierphase DD factor error" << endl;
					return false;
				}

				tempL.B.back().push_back(make_pair(idx + 1, 1.0));
				tempL.l.back() -= params_temp[idx].value();
			}
		}

		

		vector<pair<int, double>> B_L;
		double P_L, l_L;
		int ibase = 0;
		int irover = 1;
		for (const auto& b : tempL.B[irover]) {
			if (b.first > npar_orig) continue;
			B_L.push_back(b);
		}
		for (const auto& b : tempL.B[ibase]) {
			if (b.first > npar_orig) continue;
			B_L.emplace_back(b.first, -b.second);
		}
		P_L = 1 / (1 / tempL.P[irover] + 1 / tempL.P[ibase]); l_L = tempL.l[irover] - tempL.l[ibase];

		B.push_back(B_L);
		P.push_back(P_L);
		l.push_back(l_L);

		//for (int i = 0; i < B_L.size(); i++)
		//{
		//	cout << B_L[i].second << " ";
		//}
		//cout << endl;
		//cout << endl;
		//cout << endl;
	}

	int iobs = 1;
	int index_ref = 0;
	int index_sat = 1;
	DD_operator(iobs - 1, index_ref) = -1;
	DD_operator(iobs - 1, index_sat) = 1;
	Eigen::Matrix<double, 2, 5> B_new;
	Eigen::Matrix<double, 2, 2> P_new;
	Eigen::Matrix<double, 2, 1> l_new;
	trans2Eigen(B, P, l, B_new, P_new, l_new);
	B_DD = DD_operator * B_new;	
	l_DD = DD_operator * l_new;
	P_DD = DD_operator * P_new.inverse()*DD_operator.transpose();
	P_DD = 1.0 / P_DD;
	cout << "residual: " << l_DD << endl;
	cout << "Jacbian: " << B_DD << endl;
	cout << "weight: " << P_DD << endl;	
	//set ceres value
	sqrt_info = sqrt(P_DD);
	residuals[0] = sqrt_info * l_DD;
	//residuals[0] = l_DD;
	//cout << "carrier residual: " << residuals[0] << endl;
	if (jacobians)
	{
		if (jacobians[0])
		{
			Eigen::Map<Eigen::Matrix<double, 1, 3, Eigen::RowMajor>> jacobian_XYZ(jacobians[0]);
			jacobian_XYZ = -sqrt_info * B_DD.leftCols<3>();
			//jacobian_XYZ = B_DD.leftCols<3>();
			//cout << "carrier Jacbian: " << jacobian_XYZ.transpose() << endl;
		}
		/*if (jacobians[1])
		{   TODO：for velocity parameter block;


		}*/
		if (jacobians[1])
		{
			Eigen::Map<Eigen::Matrix<double, 1, 1, Eigen::RowMajor>> jacobian_ambiguity1(jacobians[1]);
			jacobian_ambiguity1 = -sqrt_info * B_DD.middleCols<1>(3);
			//jacobian_ambiguity1 = B_DD.middleCols<1>(3);
			//cout << jacobian_ambiguity1.transpose() << " ";
		}

		if (jacobians[2])
		{
			Eigen::Map<Eigen::Matrix<double, 1, 1, Eigen::RowMajor>> jacobian_ambiguity2(jacobians[2]);
			jacobian_ambiguity2 = -sqrt_info * B_DD.rightCols<1>();
			//jacobian_ambiguity2 = B_DD.rightCols<1>();
			//cout << jacobian_ambiguity2.transpose() << endl;
		}
	}
	return true;

	////////////// residuals[0] 已经算完了（你原来的 l_DD / P_DD / sqrt_info 那套）

	////////////if (!jacobians) return true;

	//////////////// ====== 保险：如果本轮线性化点变化太大，则强制重算 Jacobian ======
	////////////bool need_recompute = !_jac_cached;
	//////////////bool need_recompute = false;
	////////////// 你说“同一窗口内 Jacobian 基本不变”，通常意味着 Pi 变化很小
	////////////// 给个门限，比如 2cm。你可以按实际调大/调小。
	////////////cout << "_jac_cached:" << _jac_cached << endl;
	////////////cout << endl;
	////////////if (_jac_cached) {
	////////////	const double dp = (Pi - _lin_P).norm();
	////////////	if (dp > 0.02) need_recompute = true;  // 2cm 阈值（建议至少做这个保险）
	////////////}

	////////////if (_freeze_jac && _jac_cached && !need_recompute)
	////////////{
	////////////	if (jacobians[0]) {
	////////////		Eigen::Map<Eigen::Matrix<double, 1, 3, Eigen::RowMajor>> Jxyz(jacobians[0]);
	////////////		Jxyz(0, 0) = _J_p[0];
	////////////		Jxyz(0, 1) = _J_p[1];
	////////////		Jxyz(0, 2) = _J_p[2];
	////////////	}

	////////////	if (jacobians[1]) {
	////////////		Eigen::Map<Eigen::Matrix<double, 1, 1, Eigen::RowMajor>> Jamb1(jacobians[1]);
	////////////		Jamb1(0, 0) = _J_a1;
	////////////	}

	////////////	if (jacobians[2]) {
	////////////		Eigen::Map<Eigen::Matrix<double, 1, 1, Eigen::RowMajor>> Jamb2(jacobians[2]);
	////////////		Jamb2(0, 0) = _J_a2;
	////////////	}
	////////////	return true;
	////////////}
	////////////else 
	////////////{
	////////////	if (jacobians[0])
	////////////	{
	////////////		Eigen::Map<Eigen::Matrix<double, 1, 3, Eigen::RowMajor>> jacobian_XYZ(jacobians[0]);
	////////////		jacobian_XYZ = -sqrt_info * B_DD.leftCols<3>();
	////////////	}

	////////////	if (jacobians[1])
	////////////	{
	////////////		Eigen::Map<Eigen::Matrix<double, 1, 1, Eigen::RowMajor>> jacobian_ambiguity1(jacobians[1]);
	////////////		jacobian_ambiguity1 = -sqrt_info * B_DD.middleCols<1>(3);
	////////////	}

	////////////	if (jacobians[2])
	////////////	{
	////////////		Eigen::Map<Eigen::Matrix<double, 1, 1, Eigen::RowMajor>> jacobian_ambiguity2(jacobians[2]);
	////////////		jacobian_ambiguity2 = -sqrt_info * B_DD.rightCols<1>();
	////////////	}
	////////////}

	////////////// 写缓存
	////////////if (_freeze_jac)
	////////////{
	////////////	if (jacobians[0]) {
	////////////		_J_p[0] = jacobians[0][0];
	////////////		_J_p[1] = jacobians[0][1];
	////////////		_J_p[2] = jacobians[0][2];
	////////////	}
	////////////	else {
	////////////		_J_p[0] = _J_p[1] = _J_p[2] = 0.0;
	////////////	}

	////////////	_J_a1 = (jacobians[1] ? jacobians[1][0] : 0.0);
	////////////	_J_a2 = (jacobians[2] ? jacobians[2][0] : 0.0);

	////////////	_lin_P = Pi;
	////////////	_jac_cached = true;
	////////////}
	////////////return true;
}

void gfgo::CarrierphaseDDFactor::check(double ** parameters)
{
	double *res = new double[1];
	double **jaco = new double *[3];
	jaco[0] = new double[1 * 3];
	jaco[1] = new double[1 * 1];
	jaco[2] = new double[1 * 1];
	Evaluate(parameters, res, jaco);
	/*puts("CarrierphaseDDFactor check begins");
	puts("my: ");
	std::cout << Eigen::Map<Eigen::Matrix<double, 1, 1>>(res).transpose() << std::endl
		<< std::endl;
	std::cout << Eigen::Map<Eigen::Matrix<double, 1, 3, Eigen::RowMajor>>(jaco[0]) << std::endl
		<< std::endl;
	std::cout << Eigen::Map<Eigen::Matrix<double, 1, 1, Eigen::RowMajor>>(jaco[1]) << std::endl
		<< std::endl;
	std::cout << Eigen::Map<Eigen::Matrix<double, 1, 1, Eigen::RowMajor>>(jaco[2]) << std::endl
		<< std::endl;*/	


}
