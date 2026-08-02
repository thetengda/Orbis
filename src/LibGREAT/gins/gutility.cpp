/**
 * @file         gutility.cpp
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        implement common algorithm and define some enums
 * @version      1.0
 * @date         2025-01-01
 *
 * @copyright Copyright (c) 2025, Wuhan University. All rights reserved.
 *
 */

#include "gutility.h"
#include "gset/gsetins.h"
#include "gset/gsetign.h"
using namespace great;

int great::sign(double d)
{
	if (d > t_gglv::EPS)
		return 1;
	else if (d < -t_gglv::EPS)
		return -1;
	else
		return 0;
}

Eigen::Vector3d great::Cart2Geod(const Eigen::Vector3d& CartPos, bool b)
{
	double XYZ[3] = { CartPos(0), CartPos(1), CartPos(2) }, ell[3];
	gnut::xyz2ell(XYZ, ell, b);
	return Eigen::Vector3d(ell);
}

Eigen::Vector3d great::Geod2Cart(const Eigen::Vector3d& GeodPos, bool b)
{
	double ell[3] = { GeodPos(0), GeodPos(1), GeodPos(2) }, XYZ[3];
	gnut::ell2xyz(ell, XYZ, b);
	return Eigen::Vector3d(XYZ);
}

Eigen::Vector3d great::XYZ2ENU(const Eigen::Vector3d& xyz, const Eigen::Vector3d& xyz_ref)
{
	double XYZ[3] = { xyz(0), xyz(1), xyz(2) }, XYZ_REF[3] = { xyz_ref(0), xyz_ref(1), xyz_ref(2) }, neu[3];
	gnut::xyz2neu(XYZ, XYZ_REF, neu);
	// east north up
	return Eigen::Vector3d(neu[1], neu[0], neu[2]);
}

Eigen::MatrixXd great::NewMat2Eigen(const Matrix& newmat)
{
	int m = newmat.Nrows(), n = newmat.Ncols();
	Eigen::MatrixXd res(m, n);
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			res(i, j) = newmat(i + 1, j + 1);
		}
	}
	return res;
}

Matrix great::Eigen2newMat(const Eigen::MatrixXd& eigen)
{
	int m = eigen.rows(), n = eigen.cols();
	Matrix res;
	res.ReSize(m, n);
	res = 0.0;
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			res(i + 1, j + 1) = eigen(i, j);
		}
	}
	return res;
}

SymmetricMatrix great::Eigen2BaseMatrix(const Eigen::MatrixXd& eigen)
{
	int m = eigen.rows(), n = eigen.cols();
	SymmetricMatrix res(m);
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			res(i + 1, j + 1) = eigen(i, j);
		}
	}
	return res;
}

Eigen::MatrixXd great::BaseMatrix2Eigen(const SymmetricMatrix& newmat)
{
	int m = newmat.Nrows(), n = newmat.Ncols();
	Eigen::MatrixXd res(m, n);
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			res(i, j) = newmat(i + 1, j + 1);
		}
	}
	return res;
}

Eigen::VectorXd great::Columns2VectorXd(const ColumnVector& newmat)
{
	int m = newmat.Nrows();
	Eigen::VectorXd res(m);
	for (int i = 0; i < m; i++)
	{
		res(i) = newmat(i + 1);
	}
	return res;
}

ColumnVector great::VectorXd2Columns(const Eigen::VectorXd& eigen)
{
	int m = eigen.size();
	ColumnVector res;
	res.ReSize(m);
	res = 0.0;
	for (int i = 0; i < m; i++)
	{
		res(i + 1) = eigen(i);
	}
	return res;
}

namespace great {
	
	const double t_gglv::Re = 6378137.0;
	const double t_gglv::f = (1.0 / 298.257);
	const double t_gglv::wie = 7.2921151467e-5;
	const double t_gglv::g0 = 9.7803267714;
	const double t_gglv::mg = 1.0e-3 * g0;
	const double t_gglv::ug = 1.0e-6 * g0;
	const double t_gglv::PI = 3.14159265358979;
	const double t_gglv::EPS = 1e-10;
	const double t_gglv::INF = 1e100;
	const double t_gglv::deg = PI / 180.0;
	const double t_gglv::min = deg / 60.0;
	const double t_gglv::sec = min / 60.0;
	const double t_gglv::ppm = 1.0e-6;
	const double t_gglv::hur = 3600.0;
	const double t_gglv::dps = deg / 1.0;
	const double t_gglv::dpss = deg / sqrt(1.0);
	const double t_gglv::dph = deg / hur;
	const double t_gglv::dpsh = deg / sqrt(hur);
	const double t_gglv::dphpsh = dph / sqrt(hur);
	const double t_gglv::mpsh = 1 / sqrt(hur);
	const double t_gglv::mpspsh = 1 / 1 / sqrt(hur);
	const double t_gglv::ppmpsh = ppm / sqrt(hur);
	const double t_gglv::secpsh = sec / sqrt(hur);
	const double t_gglv::mgpsHz = mg / sqrt(1.0);
	const double t_gglv::ugpsHz = ug / sqrt(1.0);
	const double t_gglv::mgpsh = mg / sqrt(hur);
	const double t_gglv::ugpsh = ug / sqrt(hur);
	
}


great::t_scheme::t_scheme()
{
}

great::t_scheme::t_scheme(t_gsetbase* set)
{
	t = 0.0;
	ts = dynamic_cast<t_gsetins*>(set)->ts();
	freq = dynamic_cast<t_gsetins*>(set)->freq();
	nSamples = dynamic_cast<t_gsetins*>(set)->subsample();
	Cps = dynamic_cast<t_gsetins*>(set)->cps();
	align = dynamic_cast<t_gsetins*>(set)->align();
	start = dynamic_cast<t_gsetins*>(set)->start();
	end = dynamic_cast<t_gsetins*>(set)->end();
	align_time = dynamic_cast<t_gsetins*>(set)->align_time();


	_imu_inst_rot = dynamic_cast<t_gsetign*>(set)->imu_inst_rot();
	_imu_inst_trans = dynamic_cast<t_gsetign*>(set)->imu_inst_trans();
	_imu_scale = dynamic_cast<t_gsetign*>(set)->imu_scale();
	delay = dynamic_cast<t_gsetign*>(set)->delay_t();
	delay_odo = dynamic_cast<t_gsetign*>(set)->delay_odo();
	max_pdop = dynamic_cast<t_gsetign*>(set)->max_pdop();
	min_sat = dynamic_cast<t_gsetign*>(set)->min_sat();
	_NHC = dynamic_cast<t_gsetign*>(set)->NHC();
	_ZUPT = dynamic_cast<t_gsetign*>(set)->ZUPT();
	_Hgt = dynamic_cast<t_gsetign*>(set)->Hgt();
}
