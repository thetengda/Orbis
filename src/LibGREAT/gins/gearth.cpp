/**
 * @file         gearth.cpp
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        Class for Earth parameter updating
 * @version      1.0
 * @date         2025-01-01
 *
 * @copyright Copyright (c) 2025, Wuhan University. All rights reserved.
 *
 */

#include "gearth.h"

using namespace great;
using namespace Eigen;

/******************************************* t_gearth ***********************************************/
t_gearth::t_gearth(double a0, double f0, double g0)
{
    a = a0; f = f0;
    gn = Vector3d(0, 0, -g0);
    b = (1 - f)*a;
    wie = t_gglv::wie;
    e = sqrt(a*a - b*b) / a;
    e2 = e*e;
}

/************************************************************
** Function: Update Earth Parameter including Cen, wnie etc.
** Input: Geodestic position pos and local coordination velocity vn
***********************************************************/
void t_gearth::Update(const Vector3d &pos, const Vector3d &vn)
{
    this->pos = pos;  this->vn = vn;
    sb = sin(pos(0)), cb = cos(pos(0)), tb = sb / cb;
    sb2 = sb*sb, sb4 = sb2*sb2;
    sl = sin(pos(1)), cl = cos(pos(1));
    Cen << -sl, -sb*cl, cb*cl, cl, -sb*sl, cb*sl, 0, cb, sb;
    Cne = Cen.transpose();
    double sq = 1 - e2*sb*sb, sq2 = sqrt(sq);
    RMh = a*(1 - e2) / sq / sq2 + pos(2);    f_RMh = 1.0 / RMh;
    RNh = a / sq2 + pos(2);    cbRNh = cb*RNh;  f_RNh = 1.0 / RNh; f_cbRNh = 1.0 / cbRNh;
    wnie << 0, wie*cb, wie*sb;
    weie = Cen*wnie;
    wnen << -vn(1)*f_RMh, vn(0)*f_RNh, wnen(1)*tb;
    wnen(2) = wnen(1)*tb;
    wnin = wnie + wnen;
    gn(2) = -(t_gglv::g0*(1 + 5.27094e-3*sb2 + 2.32718e-5*sb4) - 3.086e-6*pos(2));
    gcc = gn - (wnie + wnin).cross(vn);
}

Vector3d t_gearth::v2dp(const Vector3d& vn, double ts)
{
    return Vector3d(vn(1)*f_RMh, vn(0)*f_cbRNh, vn(2))*ts;
}

