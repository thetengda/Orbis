/**
 * @file         gbase.cpp
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        base function
 * @version      1.0
 * @date         2025-01-01
 *
 * @copyright Copyright (c) 2025, Wuhan University. All rights reserved.
 *
 */

#include "gbase.h"
using namespace great;
using namespace Eigen;

Matrix3d t_gbase::a2mat(const Vector3d& att)
{
    double sp = sin(att(0)), cp = cos(att(0));
    double sr = sin(att(1)), cr = cos(att(1));
    double sy = sin(att(2)), cy = cos(att(2));

    Matrix3d m;
    m << cy*cr - sy*sp*sr, -sy*cp, cy*sr + sy*sp*cr,
        sy*cr + cy*sp*sr, cy*cp, sy*sr - cy*sp*cr,
        -cp*sr, sp, cp*cr;
    return m;
}

Vector3d t_gbase::m2att(const Matrix3d& m)
{
    Vector3d att;
    att(0) = asin(m(2, 1));
    att(1) = atan2(-m(2, 0), m(2, 2));
    att(2) = atan2(-m(0, 1), m(1, 1));
    return att;
}

t_gquat t_gbase::a2qua(const Vector3d& att)
{
    //double yaw = att(0) / 2.0, pitch = att(1) / 2.0, roll = att(2) / 2.0;
    double pitch = att(0) / 2.0, roll = att(1) / 2.0, yaw = att(2) / 2.0;
    double    sp = sin(pitch), sr = sin(roll), sy = sin(yaw),
        cp = cos(pitch), cr = cos(roll), cy = cos(yaw);
    t_gquat qnb;
    qnb.q0 = cp*cr*cy - sp*sr*sy;
    qnb.q1 = sp*cr*cy - cp*sr*sy;
    qnb.q2 = cp*sr*cy + sp*cr*sy;
    qnb.q3 = cp*cr*sy + sp*sr*cy;
    return qnb;
}

Vector3d t_gbase::q2att(const t_gquat& qnb)
{
    double    q11 = qnb.q0*qnb.q0, q12 = qnb.q0*qnb.q1, q13 = qnb.q0*qnb.q2, q14 = qnb.q0*qnb.q3,
        q22 = qnb.q1*qnb.q1, q23 = qnb.q1*qnb.q2, q24 = qnb.q1*qnb.q3,
        q33 = qnb.q2*qnb.q2, q34 = qnb.q2*qnb.q3,
        q44 = qnb.q3*qnb.q3;
    Vector3d att;
    //att(0) = atan2(-2 * (q23 - q14), q11 - q22 + q33 - q44);
    //att(1) = asin(2 * (q34 + q12));
    //att(2) = atan2(-2 * (q24 - q13), q11 - q22 - q33 + q44);

    att(0) = asin(2 * (q34 + q12));
    att(1) = atan2(-2 * (q24 - q13), q11 - q22 - q33 + q44);
    att(2) = atan2(-2 * (q23 - q14), q11 - q22 + q33 - q44);
    //if (att(2) > 0 && att(2) < 180)att(2) = 2 * t_gglv::PI - att(2);
    //if (att(2) > -180 && att(2) < 0)att(2) = -att(2);
    return att;
}

t_gquat t_gbase::rv2q(const Vector3d& rv)
{
#define F1    (   2 * 1)    
#define F2    (F1*2 * 2)
#define F3    (F2*2 * 3)
#define F4    (F3*2 * 4)
#define F5    (F4*2 * 5)
    double c, f, n2 = rv.norm()*rv.norm();
    if (n2 < (t_gglv::PI / 180.0*t_gglv::PI / 180.0))
    {
        double n4 = n2*n2;
        c = 1.0 - n2*(1.0 / F2) + n4*(1.0 / F4);
        f = 0.5 - n2*(1.0 / F3) + n4*(1.0 / F5);
    }
    else
    {
        double n_2 = sqrt(n2) / 2.0;
        c = cos(n_2);
        f = sin(n_2) / n_2*0.5;
    }
    return t_gquat(c, f*rv(0), f*rv(1), f*rv(2));
}

Vector3d t_gbase::q2rv(const t_gquat& q)
{
    t_gquat dq;
    dq = q;
    if (dq.q0 < 0)  { dq.q0 = -dq.q0, dq.q1 = -dq.q1, dq.q2 = -dq.q2, dq.q3 = -dq.q3; }
    if (dq.q0 > 1.0) dq.q0 = 1.0;
    double n2 = acos(dq.q0), f;
    if (n2 > 1.0e-20)
    {
        f = 2.0 / (sin(n2) / n2);
    }
    else
    {
        f = 2.0;
    }
    return Vector3d(dq.q1, dq.q2, dq.q3)*f;
}

t_gquat t_gbase::m2qua(const Matrix3d& Cnb)
{
    double q0, q1, q2, q3, qq4;
    if (Cnb(0, 0) >= Cnb(1, 1) + Cnb(2, 2))
    {
        q1 = 0.5*sqrt(1 + Cnb(0, 0) - Cnb(1, 1) - Cnb(2, 2));  qq4 = 4 * q1;
        q0 = (Cnb(2, 1) - Cnb(1, 2)) / qq4; q2 = (Cnb(0, 1) + Cnb(1, 0)) / qq4; q3 = (Cnb(0, 2) + Cnb(2, 0)) / qq4;
    }
    else if (Cnb(1, 1) >= Cnb(0, 0) + Cnb(2, 2))
    {
        q2 = 0.5*sqrt(1 - Cnb(0, 0) + Cnb(1, 1) - Cnb(2, 2));  qq4 = 4 * q2;
        q0 = (Cnb(0, 2) - Cnb(2, 0)) / qq4; q1 = (Cnb(0, 1) + Cnb(1, 0)) / qq4; q3 = (Cnb(1, 2) + Cnb(2, 1)) / qq4;
    }
    else if (Cnb(2, 2) >= Cnb(0, 0) + Cnb(1, 1))
    {
        q3 = 0.5*sqrt(1 - Cnb(0, 0) - Cnb(1, 1) + Cnb(2, 2));  qq4 = 4 * q3;
        q0 = (Cnb(1, 0) - Cnb(0, 1)) / qq4; q1 = (Cnb(0, 2) + Cnb(2, 0)) / qq4; q2 = (Cnb(1, 2) + Cnb(2, 1)) / qq4;
    }
    else
    {
        q0 = 0.5*sqrt(1 + Cnb(0, 0) + Cnb(1, 1) + Cnb(2, 2));  qq4 = 4 * q0;
        q1 = (Cnb(2, 1) - Cnb(1, 2)) / qq4; q2 = (Cnb(0, 2) - Cnb(2, 0)) / qq4; q3 = (Cnb(1, 0) - Cnb(0, 1)) / qq4;
    }
    double nq = sqrt(q0*q0 + q1*q1 + q2*q2 + q3*q3);
    q0 /= nq; q1 /= nq; q2 /= nq; q3 /= nq;
    return t_gquat(q0, q1, q2, q3);
}

Matrix3d t_gbase::q2mat(const t_gquat& qnb)
{
    double    q11 = qnb.q0*qnb.q0, q12 = qnb.q0*qnb.q1, q13 = qnb.q0*qnb.q2, q14 = qnb.q0*qnb.q3,
        q22 = qnb.q1*qnb.q1, q23 = qnb.q1*qnb.q2, q24 = qnb.q1*qnb.q3,
        q33 = qnb.q2*qnb.q2, q34 = qnb.q2*qnb.q3,
        q44 = qnb.q3*qnb.q3;
    Matrix3d Cnb;
    Cnb(0, 0) = q11 + q22 - q33 - q44, Cnb(0, 1) = 2 * (q23 - q14), Cnb(0, 2) = 2 * (q24 + q13),
        Cnb(1, 0) = 2 * (q23 + q14), Cnb(1, 1) = q11 - q22 + q33 - q44, Cnb(1, 2) = 2 * (q34 - q12),
        Cnb(2, 0) = 2 * (q24 - q13), Cnb(2, 1) = 2 * (q34 + q12), Cnb(2, 2) = q11 - q22 - q33 + q44;
    return Cnb;
}
/************************************************************
** Function: Calculate DCM according to rotation vector (Rodrigues formula)
** Input: rotation vector
** m=I+a*(rx)+b*(rx)^2
***********************************************************/
Matrix3d t_gbase::rv2m(const Vector3d& rv)
{
    double n2 = rv.norm(), a, b, n;
    if (n2 < (t_gglv::PI / 180.0*t_gglv::PI / 180.0))
    {
        a = 1 - n2*(1 / 6 - n2 / 120);
        b = 0.5 - n2*(1 / 24 - n2 / 720);
    }
    else
    {
        n = sqrt(n2);
        a = sin(n) / n;
        b = (1 - cos(n)) / n2;
    }
    Matrix3d rx = askew(rv);
    Matrix3d DCM = Matrix3d::Identity() + a*rx + b*rx*rx;
    return DCM;

}

Eigen::Matrix3d great::t_gbase::dv2mat(const Eigen::Vector3d & vb1, const Eigen::Vector3d & vb2, const Eigen::Vector3d & vn1, const Eigen::Vector3d & vn2)
{
    Eigen::Vector3d vb = vb1.cross(vb2), vn = vn1.cross(vn2);
    Eigen::Vector3d vbb = vb.cross(vb1), vnn = vn.cross(vn1);
    Eigen::Matrix3d Mb,Mn;
    Mb.block(0, 0, 1, 3) = vb1 / vb1.norm(); Mn.block(0, 0, 1, 3) = vn1 / vn1.norm();
    Mb.block(0, 0, 1, 3) = vb / vb.norm();  Mn.block(0, 0, 1, 3) = vn / vn.norm();
    Mb.block(0, 0, 1, 3) = vbb / vbb.norm(); Mn.block(0, 0, 1, 3) = vnn / vnn.norm();
    return Mb.transpose()*Mn;
}

/************************************************************
** Function: Calculate Anti-symmetric matrix
** Input: rotation vector
** [ 0, -z,  y
z,  0, -x
-y,  x,  0 ]
***********************************************************/
Matrix3d t_gbase::askew(const Vector3d& v)
{
    Matrix3d vnx;
    vnx << 0, -v(2), v(1),
        v(2), 0, -v(0),
        -v(1), v(0), 0;
    return vnx;
}

Matrix4d t_gbase::m2m4(const Vector3d& v)
{
    Matrix4d M;
    M << 0.0, -v(0), -v(1), -v(2),
        v(0), 0.0, -v(2), v(1),
        v(1), v(2), 0.0, -v(0),
        v(2), -v(1), v(0), 0.0;
    return M;
}

Matrix4d t_gbase::m2m4_(const Vector3d& v)
{
    Matrix4d M;
    M << 0.0, -v(0), -v(1), -v(2),
        v(0), 0.0, v(2), -v(1),
        v(1), -v(2), 0.0, v(0),
        v(2), v(1), -v(0), 0.0;
    return M;
}

Matrix3d t_gbase::Cen(const Vector3d& pos)
{
    double sb = sin(pos(0)), cb = cos(pos(0)), sl = sin(pos(1)), cl = cos(pos(1));
    Matrix3d cen;
    cen << -sl, -sb*cl, cb*cl, cl, -sb*sl, cb*sl, 0, cb, sb;
    return cen;
}

void t_gbase::symmetry(MatrixXd& m)
{
    m = (m + m.transpose()) / 2.0;
}

Eigen::Vector3d great::t_gbase::product(const Eigen::Vector3d& vec, const Eigen::Matrix3d& mat)
{
    Eigen::Vector3d res;
    for (int i = 0; i < 3; i++)
    {
        res(i) = vec(0)*mat(0, i) + vec(1)*mat(1, i) + vec(2)*mat(2, i);
    }
    return res;
}

void great::t_gbase::delrowcol(Eigen::MatrixXd & M, int i)
{
    Eigen::MatrixXd TMP = M;
    int m = TMP.rows(), n = TMP.cols();
    M = Eigen::MatrixXd::Zero(m - 1, n - 1);
    M.block(0, 0, i, i) = TMP.block(0, 0, i, i);
    M.block(0, i, i, n - i - 1) = TMP.block(0, i + 1, i, n - i - 1);
    M.block(i, 0, m - i - 1, i) = TMP.block(i + 1, 0, m - i - 1, i);
    M.block(i, i, m - i - 1, n - i - 1) = TMP.block(i + 1, i + 1, m - i - 1, n - i - 1);
}

void great::t_gbase::delrow(Eigen::MatrixXd & M, int i)
{
    Eigen::MatrixXd TMP = M;
    int m = TMP.rows(), n = TMP.cols();
    M = Eigen::MatrixXd::Zero(m - 1, n);
    M.block(0, 0, i, n) = TMP.block(0, 0, i, n);
    M.block(i, 0, m - i - 1, n) = TMP.block(i + 1, 0, m - i - 1, n);
}

void great::t_gbase::delcol(Eigen::MatrixXd & M, int i)
{
    Eigen::MatrixXd TMP = M;
    int m = TMP.rows(), n = TMP.cols();
    M = Eigen::MatrixXd::Zero(m, n - 1);
    M.block(0, 0, m, i) = TMP.block(0, 0, m, i);
    M.block(0, i, m, n - i - 1) = TMP.block(0, i + 1, m, n - i - 1);
}

void great::t_gbase::delrow(Eigen::VectorXd & V, int i)
{
    Eigen::VectorXd TMP = V;
    int m = TMP.rows();
    V = Eigen::VectorXd::Zero(m - 1, 1);
    V.block(0, 0, i, 1) = TMP.block(0, 0, i, 1);
    V.block(i, 0, m - i - 1, 1) = TMP.block(i + 1, 0, m - i - 1, 1);
}
