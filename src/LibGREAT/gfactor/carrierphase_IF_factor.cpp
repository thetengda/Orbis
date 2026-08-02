#include"carrierphase_IF_factor.h"
//#include"gproc/glsqmatrix.h"
gfgo::CarrierphaseIFFactor::CarrierphaseIFFactor(const t_gtime &cur_time, const string &site, const t_gallpar  &params, const t_gsatdata &IF_sat_data, t_gprecisebiasFGO *bias_model, const pair<FREQ_SEQ, GOBSBAND> &freq_band1, const pair<FREQ_SEQ, GOBSBAND> &freq_band2) :
	_cur_time(cur_time), _site(site), _params(params), _IF_sat_data(IF_sat_data), _gprecise_bias_model(bias_model), _freq_band1(freq_band1), _freq_band2(freq_band2)
{}

void gfgo::CarrierphaseIFFactor::updatePara(t_gallpar & params_tmp, const Eigen::Vector3d & Pi, const double& Clk, const double&Trp, const double&Amb) const
{
	params_tmp = _params;
	string sat = _IF_sat_data.sat();
	int i = 0;
	i = params_tmp.getParam(_site, par_type::CRD_X, "");
	if (i >= 0)
	{
#if DEBUG_GNSS_FACTOR
		cout << "before : PX : " << params_tmp[i].value() << " -> " << Pi.x() << " delta = " << Pi.x() - params_tmp[i].value() << endl;
#endif
		params_tmp[i].value(Pi.x());
	}
	else
	{
		cerr << "can not update PX for sat : " << sat << endl;
	}

	i = params_tmp.getParam(_site, par_type::CRD_Y, "");
	if (i >= 0)
	{
#if DEBUG_GNSS_FACTOR
		cout << "before : PY : " << params_tmp[i].value() << " -> " << Pi.y() << " delta = " << Pi.y() - params_tmp[i].value() << endl;
#endif
		params_tmp[i].value(Pi.y());
	}
	else
	{
		cerr << "can not update PY for sat : " << sat << endl;
	}

	i = params_tmp.getParam(_site, par_type::CRD_Z, "");
	if (i >= 0)
	{
#if DEBUG_GNSS_FACTOR
		cout << "before : PZ : " << params_tmp[i].value() << " -> " << Pi.z() << " delta = " << Pi.z() - params_tmp[i].value() << endl;
#endif
		params_tmp[i].value(Pi.z());
	}
	else
	{
		cerr << "can not update PZ for sat : " << sat << endl;
	}

	i = params_tmp.getParam(_site, par_type::CLK, "");
	if (i >= 0)
	{
#if DEBUG_GNSS_FACTOR
		cout << "before : CLK : " << params_tmp[i].value() << " -> " << Clk << " delta = " << Clk - params_tmp[i].value() << endl;
#endif
		params_tmp[i].value(Clk);
	}
	else
	{
		cerr << "can not update CLK for sat : " << sat << endl;
	}

	i = params_tmp.getParam(_site, par_type::TRP, "");
	if (i >= 0)
	{
#if DEBUG_GNSS_FACTOR
		cout << "before : TRP : " << params_tmp[i].value() << " -> " << Trp << " delta = " << Trp - params_tmp[i].value() << endl;
#endif
		params_tmp[i].value(Trp);
	}
	else
	{
		cerr << "can not update TRP for sat : " << sat << endl;
	}

	string site = _site;
	i = params_tmp.getParam(site, par_type::AMB_IF, sat);
	if (i >= 0)
	{
#if DEBUG_GNSS_FACTOR
		cout << "before : AMB : " << params_tmp[i].value() << " -> " << Amb << " delta = " << Amb - params_tmp[i].value() << endl;
#endif
		params_tmp[i].value(Amb);
	}
	else
	{
		cerr << "can not update AMB for sat : " << sat << endl;
	}
}

bool gfgo::CarrierphaseIFFactor::Evaluate(double const * const * parameters, double * residuals, double ** jacobians) const
{
	Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
	double Clk = parameters[1][0];
	double Trp = parameters[2][0];
	double Amb = parameters[3][0];//IONO FREE

	//cout << "\n<< " << _IF_sat_data.sat() << " CarrierphaseIFFactor " << endl;

	t_gallpar params_temp;
	updatePara(params_temp, Pi, Clk, Trp, Amb);
	t_gsatdata obsdata = _IF_sat_data;

	GOBSBAND b1 = _freq_band1.second;
	GOBSBAND b2 = _freq_band2.second;

	if (b1 == BAND || b2 == BAND)
	{
		return false;
	}

	t_gobs gobs1(obsdata.select_phase(b1));
	t_gobs gobs2(obsdata.select_phase(b2));

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
		cerr << "epoch : " << epoch.sow() << "sat " << obsdata.sat() << "  construct carrierphase IF factor error ( L1 )" << endl;
		return false;
	}
	if (!_gprecise_bias_model->cmb_equ(false, true, epoch, params_temp, obsdata, gobs2, temp_equ))
	{
		cerr << "epoch : " << epoch.sow() << "sat " << obsdata.sat() << "  construct carrierphase IF factor error ( L2 )" << endl;
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
			cerr << "coeff par is not the same in f1 and f2" << endl;
			return false;
		}
		// combine coeff
		coef_IF.emplace_back(temp_equ.B[0][i].first, coef1 * temp_equ.B[0][i].second + coef2 * temp_equ.B[1][i].second);
	}
	// combine P and l
	P_IF = 1.0 / (pow(coef1, 2) / temp_equ.P[0] + pow(coef2, 2) / temp_equ.P[1]);
	l_IF = coef1 * temp_equ.l[0] + coef2 * temp_equ.l[1];

	// CORRECT AMB
	if (gobs1.is_phase() && (!obsdata.is_carrier_range(b1) || !obsdata.is_carrier_range(b2)))
	{
		int idx = params_temp.getParam(obsdata.site(), par_type::AMB_IF, obsdata.sat());
		if (idx < 0)
		{
			return false;
		}
		// update B l
		coef_IF.emplace_back(idx + 1, 1.0);
		l_IF -= params_temp[idx].value();
	}

	/*if (obsdata.sat() == "G30")
	{
		P_IF = P_IF * 0.99;
	}*/
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
	for (int i = 0; i < 6; i++)
	{
		B_IF(0, i) = coef_IF[i].second;
	}
	
	if (jacobians)
	{
		//sqrt_info = -sqrt_info;
		if (jacobians[0])
		{
			Eigen::Map<Eigen::Matrix<double, 1, 3, Eigen::RowMajor>> jacobian_CRD(jacobians[0]);
			jacobian_CRD = sqrt_info * B_IF.leftCols<3>();
			//cout << "CRD jacobian = " << jacobian_CRD << endl;
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
			Eigen::Map<Eigen::Matrix<double, 1, 1, Eigen::RowMajor>> jacobian_AMB_IF(jacobians[3]);
			jacobian_AMB_IF = sqrt_info * B_IF.rightCols<1>();
			//cout << "AMB_IF jacobian = " << jacobian_AMB_IF << endl;
		}
	}

	return true;
}