#include"multi_pseudorange_IF_integration_factor.h"

gfgo::MultiPseudorangeIFINGFactor::MultiPseudorangeIFINGFactor(const t_gtime &cur_time, const string &site, const t_gallpar  &params, const t_gsatdata &IF_sat_data, t_gprecisebiasFGO *bias_model, const pair<FREQ_SEQ, GOBSBAND> &freq_band1, const pair<FREQ_SEQ, GOBSBAND> &freq_band2, const Eigen::Vector3d &lever_arm) :
	_cur_time(cur_time), _site(site), _params(params), _IF_sat_data(IF_sat_data), _gprecise_bias_model(bias_model), _freq_band1(freq_band1), _freq_band2(freq_band2), _lever_arm(lever_arm)
{}

void gfgo::MultiPseudorangeIFINGFactor::updatePara(t_gallpar & params_tmp, const Eigen::Vector3d & Pi, const double& Clk, const double&Trp, const double&ISB) const
{
	params_tmp = _params;
	string sat = _IF_sat_data.sat();
	GSYS sys = _IF_sat_data.gsys();
	int i = 0;
	i = params_tmp.getParam(_site, par_type::CRD_X, "");
	if (i >= 0)
	{
		//cout << "before : PX : " << params_tmp[i].value() << " -> " << Pi.x() << " delta = " << Pi.x() - params_tmp[i].value() << endl;
		params_tmp[i].value(Pi.x());
}
	else
	{
		cerr << "can not update PX for sat : " << sat << endl;
	}

	i = params_tmp.getParam(_site, par_type::CRD_Y, "");
	if (i >= 0)
	{
		//cout << "before : PY : " << params_tmp[i].value() << " -> " << Pi.y() << " delta = " << Pi.y() - params_tmp[i].value() << endl;
		params_tmp[i].value(Pi.y());
	}
	else
	{
		cerr << "can not update PY for sat : " << sat << endl;
	}

	i = params_tmp.getParam(_site, par_type::CRD_Z, "");
	if (i >= 0)
	{
		//cout << "before : PZ : " << params_tmp[i].value() << " -> " << Pi.z() << " delta = " << Pi.z() - params_tmp[i].value() << endl;
		params_tmp[i].value(Pi.z());
	}
	else
	{
		cerr << "can not update PZ for sat : " << sat << endl;
	}

	i = params_tmp.getParam(_site, par_type::CLK, "");
	if (i >= 0)
	{
		//cout << "before : CLK : " << params_tmp[i].value() << " -> " << Clk << " delta = " << Clk - params_tmp[i].value() << endl;
		params_tmp[i].value(Clk);
	}
	else
	{
		cerr << "can not update CLK for sat : " << sat << endl;
	}

	i = params_tmp.getParam(_site, par_type::TRP, "");
	if (i >= 0)
	{
		//cout << "before : TRP : " << params_tmp[i].value() << " -> " << Trp << " delta = " << Trp - params_tmp[i].value() << endl;
		params_tmp[i].value(Trp);
	}
	else
	{
		cerr << "can not update TRP for sat : " << sat << endl;
	}

	if (sys == GSYS::GPS)
	{
		cerr << "MultiPseudorangeIFFactor is not for GPS" << endl;
	}
	else if (sys == GSYS::GAL)
	{
		i = params_tmp.getParam(_site, par_type::GAL_ISB, "");
		if (i >= 0)
		{
			//cout << "before : ISB_GAL : " << params_tmp[i].value() << " -> " << ISB << " delta = " << ISB - params_tmp[i].value() << endl;
			params_tmp[i].value(ISB);
		}
		else
		{
			cerr << "can not update ISB for sat : " << sat << endl;
		}
	}
	else if (sys == GSYS::BDS)
	{
		i = params_tmp.getParam(_site, par_type::BDS_ISB, "");
		if (i >= 0)
		{
			//cout << "before : ISB_BDS : " << params_tmp[i].value() << " -> " << ISB << " delta = " << ISB - params_tmp[i].value() << endl;
			params_tmp[i].value(ISB);
		}
		else
		{
			cerr << "can not update ISB for sat : " << sat << endl;
		}
	}
	else if (sys == GSYS::GLO)
	{
		i = params_tmp.getParam(_site, par_type::GLO_ISB, "");
		if (i >= 0)
		{
			params_tmp[i].value(ISB);
		}
		else
		{
			cerr << "can not update ISB for sat : " << sat << endl;
		}
	}
	else if (sys == GSYS::QZS)
	{
		i = params_tmp.getParam(_site, par_type::QZS_ISB, "");
		if (i >= 0)
		{
			params_tmp[i].value(ISB);
		}
		else
		{
			cerr << "can not update ISB for sat : " << sat << endl;
		}
	}
	else
	{
		cerr << "not support the system now : " << sys << endl;
	}

	}

bool gfgo::MultiPseudorangeIFINGFactor::Evaluate(double const * const * parameters, double * residuals, double ** jacobians) const
{
	Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
	Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);
	double Clk = parameters[1][0];
	double Trp = parameters[2][0];
	double ISB = parameters[3][0];
	Eigen::Matrix3d Reb = Qi.toRotationMatrix();
	Eigen::Vector3d Pi_ARP = Pi + Reb * _lever_arm;

	t_gallpar params_temp;
	updatePara(params_temp, Pi_ARP, Clk, Trp, ISB);
	t_gsatdata obsdata = _IF_sat_data;

	GOBSBAND b1 = _freq_band1.second;
	GOBSBAND b2 = _freq_band2.second;

	if (b1 == BAND || b2 == BAND)
	{
		return false;
	}

	t_gobs gobs1(obsdata.select_range(b1));
	t_gobs gobs2(obsdata.select_range(b2));

	auto gsys = obsdata.gsys();

	// IF coef
	double coef1, coef2;
	obsdata.coef_ionofree(b1, coef1, b2, coef2);

	//t_glsqEquationMatrix equ_IF;
	//combine f1 and f2
	t_gbaseEquation temp_equ;
	t_gtime epoch = obsdata.epoch();
	if (!_gprecise_bias_model->cmb_equ(false, true, epoch, params_temp, obsdata, gobs1, temp_equ))
	{
		cerr << "epoch : " << epoch.sow() << "sat " << obsdata.sat() << "  construct pseudorange IF factor error ( P1 )" << endl;
		return false;
	}
	if (!_gprecise_bias_model->cmb_equ(false, true, epoch, params_temp, obsdata, gobs2, temp_equ))
	{
		cerr << "epoch : " << epoch.sow() << "sat " << obsdata.sat() << "  construct pseudorange IF factor error ( P2 )" << endl;
		return false;
	}

	if (temp_equ.B[0].size() != temp_equ.B[1].size())
	{
		cerr << "coeff size is not equal in f1 and f2" << endl;
		return false;
	}

	vector<pair<int, double>> coef_IF;
	double P_IF = 0.0, l_IF = 0.0;

	for (int i = 0; i < temp_equ.B[0].size(); i++)
	{
		if (temp_equ.B[0][i].first != temp_equ.B[1][i].first)
		{
			cout << "coeff par is not the same in f1 and f2" << endl;
			return false;
		}
		// combine coeff
		coef_IF.emplace_back(temp_equ.B[0][i].first, coef1 * temp_equ.B[0][i].second + coef2 * temp_equ.B[1][i].second);
	}
	// combine P and l
	P_IF = 1.0 / (pow(coef1, 2) / temp_equ.P[0] + pow(coef2, 2) / temp_equ.P[1]);
	l_IF = coef1 * temp_equ.l[0] + coef2 * temp_equ.l[1];

#if DEBUG_GNSS_FACTOR
	cout << "sat = " << obsdata.sat() << endl;
	for (int i = 0; i < coef_IF.size(); i++)
	{
		cout << "corf_IF[" << i << "] = " << coef_IF[i].first << ", " << coef_IF[i].second << endl;
	}
	cout << fixed << setprecision(8) << "P_IF = " << P_IF << endl;
	cout << fixed << setprecision(8) << "l_IF = " << l_IF << endl;
#endif

	//set ceres value
	double sqrt_info = sqrt(P_IF);
	residuals[0] = -sqrt_info * l_IF;//omc -> cmo

	Eigen::Matrix<double, 1, 6> B_IF;
	if (gsys == GSYS::GPS)
	{
		cerr << "MultiCarrierphaseIFFactor is not for GPS" << endl;
	}
	else
	{
		for (int i = 0; i < 6; i++)
		{
			B_IF(0, i) = coef_IF[i].second;
		}
	}

	if (jacobians)
	{
		//sqrt_info = -sqrt_info;
		if (jacobians[0])
		{
			Eigen::Map<Eigen::Matrix<double, 1, 7, Eigen::RowMajor>> jacobian_POSE(jacobians[0]);
			Eigen::Matrix<double, 1, 6> jacobian;
			jacobian.leftCols<3>() = B_IF.leftCols<3>();
			jacobian.rightCols<3>() = B_IF.leftCols<3>() * -Reb * t_gfgo_utility::skewSymmetric(_lever_arm);
			jacobian_POSE.leftCols<6>() = sqrt_info * jacobian;
			jacobian_POSE.rightCols<1>().setZero();
			//cout << "POSE jacobian = " << jacobian_POSE << endl;
		}
		if (jacobians[1])
		{
			Eigen::Map<Eigen::Matrix<double, 1, 1, Eigen::RowMajor>> jacobian_CLK(jacobians[1]);
			jacobian_CLK = sqrt_info * B_IF.middleCols<1>(3);
			//cout << "CLK jacobian = " << jacobian_CLK << endl;
		}
		if (jacobians[2])
		{
			Eigen::Map<Eigen::Matrix<double, 1, 1, Eigen::RowMajor>> jacobian_TRP(jacobians[2]);
			jacobian_TRP = sqrt_info * B_IF.middleCols<1>(4);
			//cout << "TRP jacobian = " << jacobian_TRP << endl;
		}
		if (jacobians[3])
		{
			Eigen::Map<Eigen::Matrix<double, 1, 1, Eigen::RowMajor>> jacobian_ISB(jacobians[3]);
			jacobian_ISB = sqrt_info * B_IF.rightCols<1>();
			//cout << "ISB jacobian = " << jacobian_ISB << endl;
		}
	}

	return true;
}
