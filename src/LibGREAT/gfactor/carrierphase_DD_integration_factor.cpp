/**
 * @file         carrierphase_DD_integration_factor.cpp
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        construct  carrier-phase factor for GNSS/INS integration
 * @version      1.0
 * @date         2025-11-04
 *
 * @copyright Copyright (c) 2025, Wuhan University. All rights reserved.
 *
 */

#include "carrierphase_DD_integration_factor.h"

gfgo::CarrierphaseDDINGFactor::CarrierphaseDDINGFactor(const t_gtime & cur_time, const pair<string, string>& base_rover_site, const t_gallpar & params, const vector<pair<t_gsatdata, t_gsatdata>>& DD_sat_data, t_gprecisebiasFGO * bias_model, const pair<FREQ_SEQ, GOBSBAND>& freq_band, const Eigen::Vector3d &lever_arm):
	_cur_time(cur_time), _base_rover_site(base_rover_site), _params(params), _DD_sat_data(DD_sat_data), _gprecise_bias_model(bias_model), _freq_band(freq_band),_lever_arm(lever_arm)
{

}

void gfgo::CarrierphaseDDINGFactor::updatePara(t_gallpar & params_tmp, const double & ref_sd_amb, const double & nonref_sd_amb, const Eigen::Vector3d & Pi, const Eigen::Vector3d & Vi) const
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

void gfgo::CarrierphaseDDINGFactor::trans2Eigen(const vector<vector<pair<int, double>>>& B, const vector<double>& P, const vector<double>& l, Eigen::Matrix<double, 2, 5>& B_new, Eigen::Matrix<double, 2, 2>& P_new, Eigen::Matrix<double, 2, 1>& l_new) const
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

bool gfgo::CarrierphaseDDINGFactor::Evaluate(double const * const * parameters, double * residuals, double ** jacobians) const
{

	Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
	Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);
	Eigen::Matrix3d Reb = Qi.toRotationMatrix();
	Eigen::Vector3d P_INS = Pi + Reb * _lever_arm;

	//Eigen::Vector3d Vi(parameters[1][0], parameters[1][1], parameters[1][2]); 
	double ref_SD_ambiguity = parameters[1][0];
	double nonref_SD_ambiguity = parameters[2][0];
	double sqrt_info;
	Eigen::Matrix<double, 1, 2> DD_operator(1, 2);
	//construct DD equ	
	unsigned npar_orig = _params.parNumber() - 5;
	t_gallpar params_temp = _params;
	t_gtriple xyz;


	updatePara(params_temp, ref_SD_ambiguity, nonref_SD_ambiguity, P_INS);

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


		//cout << "sat: " << it.first.sat() << " " << "band: " << _freq_band.first << " " << endl;
		//cout << "phase: " << "weight: " << P_L << " " << "residual: " << l_L << " " << "Jacbian: ";

		//for (int i = 0; i < B_L.size(); i++)
		//{
		//	cout << B_L[i].second << " ";
		//}
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
	//set ceres value
	sqrt_info = sqrt(P_DD);
	residuals[0] = sqrt_info * l_DD;
	if (jacobians)
	{		
		if (jacobians[0])
		{
			Eigen::Map<Eigen::Matrix<double, 1, 7, Eigen::RowMajor>> jacobian_pose(jacobians[0]);
			Eigen::Matrix<double, 1, 6> jaco_pose;
			jaco_pose.leftCols<3>() = B_DD.leftCols<3>();
			jaco_pose.rightCols<3>() = B_DD.leftCols<3>() * -Reb * t_gfgo_utility::skewSymmetric(_lever_arm);
			jacobian_pose.leftCols<6>() = -sqrt_info * jaco_pose;
			jacobian_pose.rightCols<1>().setZero();			

			//jacobian_XYZ = B_DD.leftCols<3>();
			//std::cout << "carrier Jacbian: " << jacobian_pose.transpose() << endl;
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
}



//bool gfgo::CarrierphaseDDINGFactor::Evaluate(double const * const * parameters,
//	double * residuals,
//	double ** jacobians) const
//{
//	// -------------------------
//	// 0) pose + 杆臂
//	// -------------------------
//	const Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
//	const Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);
//	const Eigen::Matrix3d Reb = Qi.toRotationMatrix();
//	const Eigen::Vector3d P_INS = Pi + Reb * _lever_arm;
//
//	const double ref_SD_ambiguity = parameters[1][0];
//	const double nonref_SD_ambiguity = parameters[2][0];
//
//	// -------------------------
//	// 1) 频点 -> ambiguity par_type（替换 map，避免每次构造）
//	// -------------------------
//	par_type amb_type;
//	switch (_freq_band.first) {
//	case FREQ_1: amb_type = par_type::AMB_L1; break;
//	case FREQ_2: amb_type = par_type::AMB_L2; break;
//	case FREQ_3: amb_type = par_type::AMB_L3; break;
//	case FREQ_4: amb_type = par_type::AMB_L4; break;
//	case FREQ_5: amb_type = par_type::AMB_L5; break;
//	default:     amb_type = par_type::AMB_L1; break;
//	}
//
//	// -------------------------
//	// 2) 临时参数（保持你原逻辑：cmb_equ 需要 params_temp）
//	// -------------------------
//	unsigned npar_orig = _params.parNumber() - 5;
//	t_gallpar params_temp = _params;
//
//	// 你原来这里更新了坐标和两颗星的 SD ambiguity
//	updatePara(params_temp, ref_SD_ambiguity, nonref_SD_ambiguity, P_INS);
//
//	// -------------------------
//	// 3) 我们直接生成 trans2Eigen 所需的 (B_new, P_new, l_new)
//	//    避免构造 vector<vector<pair>> B/P/l 和 trans2Eigen() 的二次搬运
//	// -------------------------
//	Eigen::Matrix<double, 2, 5> B_new; B_new.setZero();
//	Eigen::Matrix<double, 2, 2> P_new; P_new.setZero();
//	Eigen::Matrix<double, 2, 1> l_new; l_new.setZero();
//
//	// _DD_sat_data 你原代码既 for(auto it: ...) 又 _DD_sat_data[0/1]，
//	// 合理推断它 size=2（ref, nonref），每个元素是 pair<base_obs, rover_obs> 或同类
//	if (_DD_sat_data.size() < 2) return false;
//
//	for (int isat = 0; isat < 2; ++isat)
//	{
//		const auto& rec_pair = _DD_sat_data[isat];
//		t_gbaseEquation tempL;
//
//		// base/rover 两站
//		for (int isite = 0; isite < 2; ++isite)
//		{
//			t_gsatdata* satdata_ptr = (isite == 0)
//				? const_cast<t_gsatdata*>(&rec_pair.first)
//				: const_cast<t_gsatdata*>(&rec_pair.second);
//
//			t_gobs  obsL = t_gobs(satdata_ptr->select_phase(_freq_band.second));
//			t_gtime crt = satdata_ptr->epoch();
//
//			if (!_gprecise_bias_model->cmb_equ(crt, params_temp, *satdata_ptr, obsL, tempL))
//			{
//				std::cout << "sat " << rec_pair.second.sat()
//					<< " construct carrierphase DD factor error" << std::endl;
//				return false;
//			}
//
//			// rover 站：添加该星该频点 ambiguity 项（你原逻辑）
//			if (satdata_ptr->site() == _base_rover_site.second)
//			{
//				const int idx = params_temp.getParam(satdata_ptr->site(), amb_type, satdata_ptr->sat());
//				if (idx < 0)
//				{
//					std::cout << "sat " << rec_pair.second.sat()
//						<< " construct carrierphase DD factor error" << std::endl;
//					return false;
//				}
//				tempL.B.back().push_back(std::make_pair(idx + 1, 1.0));
//				tempL.l.back() -= params_temp[idx].value();
//			}
//		}
//
//		// 组合单差：rover - base（与你原来一致）
//		const int ibase = 0;
//		const int irover = 1;
//
//		// 只保留 <= npar_orig 的项（你原来就做了这个过滤）
//		// 并且保持顺序（因为 trans2Eigen 完全依赖顺序）
//		std::vector<std::pair<int, double>> B_L;
//		B_L.reserve(4); // 通常只需要 XYZ(3) + amb(1)
//
//		for (const auto& b : tempL.B[irover]) {
//			if (b.first > (int)npar_orig) continue;
//			B_L.push_back(b);
//		}
//		for (const auto& b : tempL.B[ibase]) {
//			if (b.first > (int)npar_orig) continue;
//			B_L.emplace_back(b.first, -b.second);
//		}
//
//		const double l_L = tempL.l[irover] - tempL.l[ibase];
//		const double P_L = 1.0 / (1.0 / tempL.P[irover] + 1.0 / tempL.P[ibase]);
//
//		// --------- 关键：直接按 trans2Eigen 的规则填矩阵 ---------
//		// trans2Eigen 假设：B_L[0..2] 是 xyz，B_L[3] 是 ambiguity
//		// 这里保持完全一致
//		B_new(isat, 0) = B_L[0].second;
//		B_new(isat, 1) = B_L[1].second;
//		B_new(isat, 2) = B_L[2].second;
//
//		if (isat == 0) B_new(isat, 3) = B_L[3].second; // ref 行 -> col3
//		else           B_new(isat, 4) = B_L[3].second; // nonref 行 -> col4
//
//		P_new(isat, isat) = P_L;   // 权（与你 trans2Eigen 一致）
//		l_new(isat) = l_L;
//	}
//
//	// -------------------------
//	// 4) 直接做 DD（sat - ref），替代 DD_operator * xxx
//	// -------------------------
//	const Eigen::Matrix<double, 1, 5> B_DD = B_new.row(1) - B_new.row(0);
//	const double l_DD = l_new(1) - l_new(0);
//
//	// 你原实现：P_DD = 1 / (D * P^{-1} * D^T)
//	// 由于 P_new 是对角权阵 => P^{-1} 是对角方差阵
//	// 所以 var_DD = var_sat + var_ref = 1/P_sat + 1/P_ref
//	const double var_ref = 1.0 / P_new(0, 0);
//	const double var_sat = 1.0 / P_new(1, 1);
//	const double var_DD = var_ref + var_sat;
//	const double P_DD = 1.0 / var_DD;
//
//	const double sqrt_info = std::sqrt(P_DD);
//	residuals[0] = sqrt_info * l_DD;
//
//	// -------------------------
//	// 5) Jacobians（保持你原来的形式）
//	// -------------------------
//	if (jacobians)
//	{
//		if (jacobians[0])
//		{
//			Eigen::Map<Eigen::Matrix<double, 1, 7, Eigen::RowMajor>> jacobian_pose(jacobians[0]);
//			Eigen::Matrix<double, 1, 6> jaco_pose;
//
//			jaco_pose.leftCols<3>() = B_DD.leftCols<3>();
//			jaco_pose.rightCols<3>() = B_DD.leftCols<3>() * (-Reb) * t_gfgo_utility::skewSymmetric(_lever_arm);
//
//			jacobian_pose.leftCols<6>() = -sqrt_info * jaco_pose;
//			jacobian_pose.rightCols<1>().setZero();
//		}
//
//		if (jacobians[1]) {
//			jacobians[1][0] = -sqrt_info * B_DD(0, 3);
//		}
//		if (jacobians[2]) {
//			jacobians[2][0] = -sqrt_info * B_DD(0, 4);
//		}
//	}
//
//	return true;
//}




void gfgo::CarrierphaseDDINGFactor::check(double ** parameters)
{
	double *res = new double[1];
	double **jaco = new double *[3];
	jaco[0] = new double[1 * 7];
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
