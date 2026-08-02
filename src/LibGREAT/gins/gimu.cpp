/**
 * @file         gimu.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        IMU core algorithm
 * @version      1.0
 * @date         2025-01-01
 *
 * @copyright Copyright (c) 2025, Wuhan University. All rights reserved.
 *
 */

#include "gimu.h"
#include <iostream>
#include <iomanip>
#include "gset/gsetout.h"
#include "gset/gsetins.h"
#include <Eigen/src/Core/util/DisableStupidWarnings.h>
using namespace great;
using namespace Eigen;

t_gimu::t_gimu()
{
    phim = dvbm = wm_1 = vm_1 = Vector3d::Zero();
    pf1 = pf2 = 1;
}

Vector3d t_gimu::phi_cone_crt(const vector<Vector3d>& wm)
{
    Vector3d dphim;
    int n = wm.size();
    if (n == 1)
        dphim = Vector3d::Zero();
    else if (n == 2)
        dphim = 2. / 3 * wm[0].cross(wm[1]);
    else if (n == 3)
        dphim = 27. / 40 * wm[1].cross(wm[2]) + 9. / 20 * wm[0].cross(wm[2]) + 27. / 40 * wm[0].cross(wm[1]);
    else if (n == 4)
        dphim = 232. / 315 * wm[2].cross(wm[3]) + 46. / 105 * wm[1].cross(wm[3]) +
        18. / 35 * wm[0].cross(wm[3]) + 178. / 315 * wm[1].cross(wm[2]) +
        46. / 105 * wm[0].cross(wm[2]) + 232. / 315 * wm[0].cross(wm[1]);
    else if (n == 5)
        dphim = 18575. / 24192 * wm[3].cross(wm[4]) + 2675. / 6048 * wm[2].cross(wm[4]) +
        11225. / 24192 * wm[1].cross(wm[4]) + 125. / 252 * wm[0].cross(wm[4]) +
        2575. / 6048 * wm[2].cross(wm[3]) + 425. / 672 * wm[1].cross(wm[3]) +
        13975. / 24192 * wm[0].cross(wm[3]) + 1975. / 3024 * wm[1].cross(wm[2]) +
        325. / 1512 * wm[0].cross(wm[2]) + 21325. / 24192 * wm[0].cross(wm[1]);
    else
        dphim = Vector3d::Zero();
    return dphim;
}

Vector3d t_gimu::phi_poly_crt(const vector<Vector3d>& wm)
{
    Vector3d dphim;
    int n = wm.size();
    if (n == 1)
    {
        // One Plus Previous
        if (pf1 == 1){ wm_1 = wm[0]; pf1 = 0; }
        dphim = 1. / 12 * wm_1.cross(wm[0]);
    }
    else if (n == 2)
        dphim = 2. / 3 * wm[0].cross(wm[1]);
    else if (n == 3)
        dphim = 33. / 80 * wm[0].cross(wm[2]) + 57. / 80 * wm[1].cross(wm[2] - wm[0]);
    else if (n == 4)
        dphim = 736. / 945 * (wm[0].cross(wm[1]) + wm[2].cross(wm[3])) +
        334. / 945 * (wm[0].cross(wm[2]) + wm[1].cross(wm[3])) +
        526. / 945 * wm[0].cross(wm[3]) + 654. / 945 * wm[1].cross(wm[2]);
    else if (n == 5)
        dphim = 123425. / 145152 * (wm[0].cross(wm[1]) + wm[3].cross(wm[4])) +
        34875. / 145152 * (wm[0].cross(wm[2]) + wm[2].cross(wm[4])) +
        90075. / 145152 * (wm[0].cross(wm[3] + wm[1].cross(wm[4]))) +
        66625. / 145152 * wm[0].cross(wm[4]) +
        103950. / 145152 * (wm[1].cross(wm[2]) + wm[2].cross(wm[3])) +
        55400. / 145152 * (wm[1].cross(wm[3]));
    else
        dphim = Vector3d::Zero();  
    return dphim;
}

Vector3d t_gimu::vm_cone_crt(const vector<Vector3d>& wm, const vector<Vector3d>& vm)
{
    Vector3d scullm;
    int n = wm.size();
    if (n == 1)
        scullm = Vector3d::Zero();
    else if (n == 2)
        scullm = 2. / 3 * (wm[0].cross(vm[1]) + vm[0].cross(wm[1]));
    else if (n == 3)
        scullm = 27. / 40 * (wm[1].cross(vm[2]) + vm[1].cross(wm[2])) +
        9. / 20 * (wm[0].cross(vm[2]) + vm[0].cross(wm[2])) +
        27. / 40 * (wm[0].cross(vm[1]) + vm[0].cross(wm[1]));
    else if (n == 4)
        scullm = 232. / 315 * (wm[2].cross(vm[3]) + vm[2].cross(wm[3])) +
        46. / 105 * (wm[1].cross(vm[3]) + vm[1].cross(wm[3])) +
        18. / 35 * (wm[0].cross(vm[3]) + vm[0].cross(wm[3])) +
        178. / 315 * (wm[1].cross(vm[2]) + vm[1].cross(wm[2])) +
        46. / 105 * (wm[0].cross(vm[2]) + vm[0].cross(wm[2])) +
        232. / 315 * (wm[0].cross(vm[1]) + vm[0].cross(wm[1]));
    return scullm;
}

Vector3d t_gimu::vm_poly_crt(const vector<Vector3d>& wm, const vector<Vector3d>& vm)
{
    Vector3d scullm;
    int n = wm.size();
    if (n == 1)
    {
        // One Plus Previous
        if (pf2 == 1){ wm_1 = wm[0]; vm_1 = vm[0]; pf2 = 0; }
        scullm = 1. / 12 * (wm_1.cross(vm[0]) + vm_1.cross(wm[0]));
        wm_1 = wm[0]; vm_1 = vm[0];
    }
    else if (n == 2)
        scullm = 2. / 3 * (wm[0].cross(vm[1]) + vm[0].cross(wm[1]));
    else if (n == 3)
        scullm = 33. / 80 * (wm[0].cross(vm[2]) + vm[0].cross(wm[2])) +
        57. / 80 * (wm[0].cross(vm[1]) + vm[0].cross(wm[1]) + wm[1].cross(vm[2] + vm[1].cross(wm[2])));
    else if (n == 4)
        scullm = 736. / 945 * (wm[0].cross(vm[1]) + wm[2].cross(vm[3]) + vm[0].cross(wm[1]) + vm[2].cross(wm[3])) +
        334. / 945 * (wm[0].cross(vm[2]) + wm[1].cross(vm[3]) + vm[0].cross(wm[2]) + vm[1].cross(wm[3])) +
        526. / 945 * (wm[0].cross(vm[3]) + vm[0].cross(wm[3])) +
        654. / 945 * (wm[1].cross(vm[2]) + vm[1].cross(wm[2]));
    else
        scullm = Vector3d::Zero();
    return scullm;

}

void t_gimu::Update(const vector<Vector3d>& wm, const vector<Vector3d>& vm, const t_scheme& scm)
{
    Vector3d dphim, scullm, wmm = Vector3d::Zero(), vmm = Vector3d::Zero();
    for (size_t i = 0; i < wm.size(); i++)
    {
        wmm += wm[i];
        vmm += vm[i];
    }
    if (scm.Cps == 0)
    {
        dphim = phi_cone_crt(wm);
        scullm = vm_cone_crt(wm, vm);
    }
    else if (scm.Cps == 2)
    {
        dphim = phi_poly_crt(wm);
        scullm = vm_poly_crt(wm, vm);
    }
    Vector3d rotm = 1. / 2 * wmm.cross(vmm);
    phim = wmm + dphim;
    dvbm = vmm + rotm + scullm;
}



