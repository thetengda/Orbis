/**
 * @file         gquat.cpp
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        Class for Quaternion
 * @version      1.0
 * @date         2025-01-01
 *
 * @copyright Copyright (c) 2025, Wuhan University. All rights reserved.
 *
 */

#include "gquat.h"
#include "gbase.h"
using namespace great;
using namespace Eigen;

/******************************************* t_gquat ***********************************************/
t_gquat::t_gquat(double q0, double q1, double q2, double q3)
{
    this->q0 = q0; this->q1 = q1;
    this->q2 = q2; this->q3 = q3;
}

t_gquat::t_gquat(const Vector4d& m)
{
    q0 = m(0), q1 = m(1), q2 = m(2), q3 = m(3);
}

t_gquat t_gquat::operator+(const t_gquat &q) const
{
    return t_gquat(q0 + q.q0, q1 + q.q1, q2 + q.q2, q3 + q.q3);
}

t_gquat t_gquat::operator+(const Vector3d& phi) const
{
    t_gquat qtmp = t_gbase::rv2q(-phi);
    return qtmp*(*this);
}

t_gquat t_gquat::operator-(const Vector3d &phi) const
{
    t_gquat qtmp = t_gbase::rv2q(phi);
    return qtmp*(*this);
}

t_gquat t_gquat::operator*(const t_gquat &q) const
{
    t_gquat qtmp;
    qtmp.q0 = q0*q.q0 - q1*q.q1 - q2*q.q2 - q3*q.q3;
    qtmp.q1 = q0*q.q1 + q1*q.q0 + q2*q.q3 - q3*q.q2;
    qtmp.q2 = q0*q.q2 + q2*q.q0 + q3*q.q1 - q1*q.q3;
    qtmp.q3 = q0*q.q3 + q3*q.q0 + q1*q.q2 - q2*q.q1;
    return qtmp;
}

Vector3d t_gquat::operator-(const t_gquat &quat) const
{
    t_gquat dq;

    dq = quat*(t_gquat::conj(*this));
    if (dq.q0 < 0)
    {
        dq.q0 = -dq.q0, dq.q1 = -dq.q1, dq.q2 = -dq.q2, dq.q3 = -dq.q3;
    }
    double n2 = acos(dq.q0), f;
    if (sign(n2) != 0)
    {
        f = 2.0 / (sin(n2) / n2);
    }
    else
    {
        f = 2.0;
    }
    return Vector3d(dq.q1, dq.q2, dq.q3)*f;
}

Vector3d t_gquat::operator*(const Vector3d &v) const
{
    t_gquat qtmp;
    Vector3d vtmp;
    qtmp.q0 = -q1*v(0) - q2*v(1) - q3*v(2);
    qtmp.q1 = q0*v(0) + q2*v(2) - q3*v(1);
    qtmp.q2 = q0*v(1) + q3*v(0) - q1*v(2);
    qtmp.q3 = q0*v(2) + q1*v(1) - q2*v(0);
    vtmp(0) = -qtmp.q0*q1 + qtmp.q1*q0 - qtmp.q2*q3 + qtmp.q3*q2;
    vtmp(1) = -qtmp.q0*q2 + qtmp.q2*q0 - qtmp.q3*q1 + qtmp.q1*q3;
    vtmp(2) = -qtmp.q0*q3 + qtmp.q3*q0 - qtmp.q1*q2 + qtmp.q2*q1;
    return vtmp;
}

t_gquat& t_gquat::operator*=(const t_gquat &q)
{
    return (*this = *this*q);
}

void t_gquat::normlize(t_gquat& q)
{
    double nq = sqrt(q.q0*q.q0 + q.q1*q.q1 + q.q2*q.q2 + q.q3*q.q3);
    q.q0 /= nq, q.q1 /= nq, q.q2 /= nq, q.q3 /= nq;
}

t_gquat t_gquat::conj(const t_gquat& q)
{
    return t_gquat(q.q0, -q.q1, -q.q2, -q.q3);
}

Eigen::Matrix4d great::t_gquat::left()
{
    Eigen::Matrix4d mat;
    mat << q0, -q1, -q2, -q3,
        q1, q0, -q3, q2,
        q2, q3, q0, -q1,
        q3, -q2, q1, q0;
    return mat;
}

Eigen::Matrix4d great::t_gquat::right()
{
    Eigen::Matrix4d mat;
    mat << q0, -q1, -q2, -q3,
        q1, q0, q3, -q2,
        q2, -q3, q0, q1,
        q3, q2, -q1, q0;
    return mat;
}

Eigen::Vector3d great::t_gquat::vec()
{
    return Eigen::Vector3d(q1, q2, q3);
}
