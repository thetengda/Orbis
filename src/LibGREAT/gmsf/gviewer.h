/**
 * @file         gviewer.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        visualization
 * @version      1.0
 * @date         2025-01-01
 *
 * @copyright Copyright (c) 2025, Wuhan University. All rights reserved.
 *
 */

#ifndef GVIEWER_H
#define GVIEWER_H

#include <Eigen/Core>
#include <GLFW/glfw3.h>
#include <vector>
#include <utility>
#include <thread>
#include <mutex>
#include <math.h>
#include <map>
#ifdef _WIN32
#include <Windows.h>
#else
#include <unistd.h>
#endif
#include "gexport/ExportLibGREAT.h"

using namespace std;
typedef vector<pair<Eigen::Matrix3d, Eigen::Vector3d>> Frames;
typedef vector<Eigen::Vector3d> VPointCloud;
typedef vector<Eigen::Vector3d> Trajectory;

namespace great
{
    class LibGREAT_LIBRARY_EXPORT gviewer
    {
    public:
        gviewer();
        ~gviewer();
    public:
        void Show();
        void ClearView();
        void SetFrames(const Frames &t);
        void SetFrames(const vector<Frames> &vt);
        void AddNewFrame(const pair<Eigen::Matrix3d, Eigen::Vector3d> &f);

        void SetPointCloud(const VPointCloud &pc);
        void SetPointCloud(const vector<VPointCloud> &vpc);
        void SetPlaneCloud(const vector<Eigen::Vector3d> &_pcs, const vector<vector<Eigen::Vector3d>> & _near_points);
        void AddNewPoint(const Eigen::Vector3d &p);
        void AddNewPoint(const vector<Eigen::Vector3d> &p);

        void SetTrajectory(const Trajectory &t);
        void SetTrajectory(const vector<Trajectory> &vt);
        void AddNewPos(const Eigen::Vector3d &p);
        void Hide();
    private:
        void Run();

    private:
        static vector<gviewer*> vptr;
        GLFWwindow* window;
        thread *t;
        mutex m_mutex;
        vector<Trajectory> mv_trajectory;
        vector<VPointCloud> mv_pointCloud;
        vector<Eigen::Vector3d> pcs;
        vector<vector<Eigen::Vector3d>> near_points;
        vector<Frames> mv_frames;
    };
}
#endif