/**
 * @file         pseudorange_DD_integration_factor.cpp
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        construct  pseudo-range factor for GNSS/INS integration
 * @version      1.0
 * @date         2025-11-04
 *
 * @copyright Copyright (c) 2025, Wuhan University. All rights reserved.
 *
 */

#include "pseudorange_DD_integration_factor.h"

gfgo::PseudorangeDDINGFactor::PseudorangeDDINGFactor(const t_gtime & cur_time, const pair<string, string>& base_rover_site, const t_gallpar & params, const vector<pair<t_gsatdata, t_gsatdata>>& DD_sat_data, t_gprecisebiasFGO * bias_model, const pair<FREQ_SEQ, GOBSBAND>& freq_band, const Eigen::Vector3d &lever_arm):
	_cur_time(cur_time), _base_rover_site(base_rover_site), _params(params), _DD_sat_data(DD_sat_data), _gprecise_bias_model(bias_model), _freq_band(freq_band),  _lever_arm(lever_arm)
{



}
void gfgo::PseudorangeDDINGFactor::updatePara(t_gallpar & params_tmp, const Eigen::Vector3d & Pi, const Eigen::Vector3d & Vi) const
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
void gfgo::PseudorangeDDINGFactor::trans2Eigen(const vector<vector<pair<int, double>>>& B, const vector<double>& P, const vector<double>& l, Eigen::Matrix<double, 2, 3>& B_new, Eigen::Matrix<double, 2, 2>& P_new, Eigen::Matrix<double, 2, 1>& l_new) const
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
//bool gfgo::PseudorangeDDINGFactor::Evaluate(double const * const * parameters, double * residuals, double ** jacobians) const
//{
//	/*if (_cur_time.sow() == 200340 && _DD_sat_data.at(0).first.sat() == "G10" && _DD_sat_data.at(1).first.sat() == "G25")
//	{
//		cout << endl;
//	}*/
//
//	Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
//	Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);
//	//Eigen::Vector3d Vi(parameters[1][0], parameters[1][1], parameters[1][2])
//	Eigen::Matrix3d Reb = Qi.toRotationMatrix();
//	Eigen::Vector3d P_INS = Pi + Reb * _lever_arm;
//	double sqrt_info;
//	Eigen::Matrix<double, 1, 2> DD_operator(1, 2);
//	//construct DD equ	
//	unsigned npar_orig = _params.parNumber() - 5;
//	t_gallpar params_temp;
//	updatePara(params_temp, P_INS);
//	vector<vector<pair<int, double>>> B;        ///< coeff of equations
//	vector<double> P;                           ///< weight of equations
//	vector<double> l;                           ///< res of equations
//	Eigen::MatrixXd B_DD;
//	double l_DD, P_DD;
//	for (auto it : _DD_sat_data)
//	{
//		t_gbaseEquation tempP;
//		pair<t_gsatdata, t_gsatdata> rec_pair = it;
//		for (int isite = 0; isite < 2; isite++)
//		{
//			t_gsatdata *satdata_ptr;
//			if (isite == 0) satdata_ptr = &rec_pair.first;
//			else satdata_ptr = &rec_pair.second;
//			t_gobs  obsP = t_gobs(satdata_ptr->select_range(_freq_band.second));
//			t_gtime crt = satdata_ptr->epoch();
//			if (!_gprecise_bias_model->cmb_equ(crt, params_temp, *satdata_ptr, obsP, tempP))
//			{
//				cout << "sat " << rec_pair.second.sat() << "  construct pseudorange DD factor error" << endl;
//				return false;
//
//			}
//		}
//		vector<pair<int, double>> B_P;
//		double P_P, l_P;
//		int ibase = 0;
//		int irover = 1;
//		for (const auto& b : tempP.B[irover]) {
//			if (b.first > npar_orig) continue;
//			B_P.push_back(b);
//		}
//		for (const auto& b : tempP.B[ibase]) {
//			if (b.first > npar_orig) continue;
//			B_P.emplace_back(b.first, -b.second);
//		}
//		P_P = 1 / (1 / tempP.P[irover] + 1 / tempP.P[ibase]); l_P = tempP.l[irover] - tempP.l[ibase];
//
//		B.push_back(B_P);
//		P.push_back(P_P);
//		l.push_back(l_P);
//
//#ifndef DEBUG_FACTOR
//		//cout << "sat: " << it.first.sat() << " " << "band: " << _freq_band.first << " " << endl;
//		//cout << "code: " << "weight: " << P_P << " " << "residual: " << l_P << " " << "Jacbian: ";
//		/*for (int i = 0; i < B_P.size(); i++)
//		{
//			cout << B_P[i].second << " ";
//		}
//		cout << endl;*/
//#endif // !1
//	}
//	int iobs = 1;
//	int index_ref = 0;
//	int index_sat = 1;
//	DD_operator(iobs - 1, index_ref) = -1;
//	DD_operator(iobs - 1, index_sat) = 1;
//	Eigen::Matrix<double, 2, 3> B_new;
//	Eigen::Matrix<double, 2, 2> P_new;
//	Eigen::Matrix<double, 2, 1> l_new;
//	trans2Eigen(B, P, l, B_new, P_new, l_new);
//	B_DD = DD_operator * B_new;
//	l_DD = DD_operator * l_new;
//	P_DD = DD_operator * P_new.inverse()*DD_operator.transpose();
//	P_DD = 1.0 / P_DD;
//	//cout << "residual: " << l_DD << endl;
//	//cout << "Jacbian: " << B_DD << endl;
//	//cout << "weight: " << P_DD << endl;
//	//cout << endl;
//	//set ceres value
//	sqrt_info = sqrt(P_DD);
//	residuals[0] = sqrt_info * l_DD;
//
//	if (jacobians)
//	{
//
//		if (jacobians[0])
//		{
//			Eigen::Map<Eigen::Matrix<double, 1, 7, Eigen::RowMajor>> jacobian_pose(jacobians[0]);
//			Eigen::Matrix<double, 1, 6> jaco_pose;
//			jaco_pose.leftCols<3>() = B_DD.leftCols<3>();
//			jaco_pose.rightCols<3>() = B_DD.leftCols<3>() * -Reb * t_gfgo_utility::skewSymmetric(_lever_arm);
//			jacobian_pose.leftCols<6>() = -sqrt_info * jaco_pose;
//			jacobian_pose.rightCols<1>().setZero();
//
//			/*if (_cur_time.sow() == 200340 && _DD_sat_data.at(0).first.sat() == "G10" && _DD_sat_data.at(1).first.sat() == "G25")
//				cout << jacobian_pose << endl;*/
//
//			//std::cout << "range Jacbian: " << jacobian_pose.transpose() << endl;
//		}
//		/*if (jacobians[1])
//		{   TODO：for velocity parameter block;
//
//
//		}*/
//
//	}
//	return true;
//}

bool gfgo::PseudorangeDDINGFactor::Evaluate(double const * const * parameters,
	double * residuals,
	double ** jacobians) const
{
	// -------------------------
	// 0) pose + 杆臂
	// -------------------------
	const Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
	const Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);
	const Eigen::Matrix3d Reb = Qi.toRotationMatrix();
	const Eigen::Vector3d P_INS = Pi + Reb * _lever_arm;

	// -------------------------
	// 1) params_temp（保持你原逻辑：cmb_equ 需要它）
	// -------------------------
	const unsigned npar_orig = _params.parNumber() - 5;

	t_gallpar params_temp;
	updatePara(params_temp, P_INS);

	// -------------------------
	// 2) 直接构造两条 SD（ref / nonref）的 B_row(1x3), l, P
	//    不再构造 vector<vector<pair>> B/P/l，也不再 trans2Eigen()
	// -------------------------
	if (_DD_sat_data.size() < 2) return false;

	Eigen::Matrix<double, 2, 3> B_new; B_new.setZero(); // 两条 SD 的 xyz 系数
	Eigen::Matrix<double, 2, 1> l_new; l_new.setZero(); // 两条 SD 残差
	Eigen::Matrix<double, 2, 2> P_new; P_new.setZero(); // 两条 SD 权（对角）

	for (int isat = 0; isat < 2; ++isat)
	{
		const auto& rec_pair = _DD_sat_data[isat];
		t_gbaseEquation tempP;

		// base/rover 两站
		for (int isite = 0; isite < 2; ++isite)
		{
			t_gsatdata* satdata_ptr = (isite == 0)
				? const_cast<t_gsatdata*>(&rec_pair.first)
				: const_cast<t_gsatdata*>(&rec_pair.second);

			t_gobs  obsP = t_gobs(satdata_ptr->select_range(_freq_band.second));
			t_gtime crt = satdata_ptr->epoch(); 
			if (!_gprecise_bias_model->cmb_equ(true, true, crt, params_temp, *satdata_ptr, obsP, tempP))
			{
				std::cout << "sat " << rec_pair.second.sat()
					<< " construct pseudorange DD factor error" << std::endl;
				return false;
			}
		}

		// SD：rover - base（与你原来一致）
		const int ibase = 0;
		const int irover = 1;

		// 你原来是把两个站的项拼到 B_P，再 trans2Eigen 只取前三项(second)
		// 为保持与你 trans2Eigen 一致，这里仍按“顺序 + 过滤”生成 B_P，
		// 然后直接取 B_P[0..2].second 填入 B_new。
		std::vector<std::pair<int, double>> B_P;
		B_P.reserve(3); // 伪距只需要 xyz 三项

		for (const auto& b : tempP.B[irover]) {
			if (b.first > (int)npar_orig) continue;
			B_P.push_back(b);
		}
		for (const auto& b : tempP.B[ibase]) {
			if (b.first > (int)npar_orig) continue;
			B_P.emplace_back(b.first, -b.second);
		}

		const double l_P = tempP.l[irover] - tempP.l[ibase];
		const double P_P = 1.0 / (1.0 / tempP.P[irover] + 1.0 / tempP.P[ibase]);

		// 按你 trans2Eigen 的假设：B_P[0..2] 就是 xyz 系数
		// 为安全起见，做一下长度检查（避免越界）
		if ((int)B_P.size() < 3) {
			std::cout << "sat " << rec_pair.second.sat()
				<< " pseudorange SD equation has insufficient terms" << std::endl;
			return false;
		}

		B_new(isat, 0) = B_P[0].second;
		B_new(isat, 1) = B_P[1].second;
		B_new(isat, 2) = B_P[2].second;

		l_new(isat, 0) = l_P;
		P_new(isat, isat) = P_P;
	}

	// -------------------------
	// 3) 直接做 DD（sat - ref），替代 DD_operator * xxx
	// -------------------------
	const Eigen::Matrix<double, 1, 3> B_DD = B_new.row(1) - B_new.row(0);
	const double l_DD = l_new(1, 0) - l_new(0, 0);

	// 你原来：P_DD = 1 / (D * P^{-1} * D^T)
	// P_new 是对角权阵 => var_DD = 1/P_ref + 1/P_sat
	const double var_ref = 1.0 / P_new(0, 0);
	const double var_sat = 1.0 / P_new(1, 1);
	const double var_DD = var_ref + var_sat;
	const double P_DD = 1.0 / var_DD;

	const double sqrt_info = std::sqrt(P_DD);
	residuals[0] = sqrt_info * l_DD;

	// -------------------------
	// 4) Jacobians（保持你原来的形式）
	// -------------------------
	if (jacobians)
	{
		if (jacobians[0])
		{
			Eigen::Map<Eigen::Matrix<double, 1, 7, Eigen::RowMajor>> jacobian_pose(jacobians[0]);
			Eigen::Matrix<double, 1, 6> jaco_pose;

			jaco_pose.leftCols<3>() = B_DD.leftCols<3>();
			jaco_pose.rightCols<3>() = B_DD.leftCols<3>() * (-Reb) * t_gfgo_utility::skewSymmetric(_lever_arm);

			jacobian_pose.leftCols<6>() = -sqrt_info * jaco_pose;
			jacobian_pose.rightCols<1>().setZero();
		}
	}

	return true;
}



void gfgo::PseudorangeDDINGFactor::check(double ** parameters)
{
	double *res = new double[1];
	double **jaco = new double *[2];
	jaco[0] = new double[1 * 7];
	Evaluate(parameters, res, jaco);
	/*puts("check begins");

	puts("my");*/

	/*std::cout << Eigen::Map<Eigen::Matrix<double, 1, 1>>(res).transpose() << std::endl
		<< std::endl;
	std::cout << Eigen::Map<Eigen::Matrix<double, 1, 3, Eigen::RowMajor>>(jaco[0]) << std::endl
		<< std::endl;*/
		
}
