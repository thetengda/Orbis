/**
 * @file         gpublish.cpp
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        visualization
 * @version      1.0
 * @date         2025-01-01
 *
 * @copyright Copyright (c) 2025, Wuhan University. All rights reserved.
 *
 */
#include "gpublish.h"


great::t_gpublish::t_gpublish()
{
    
}

great::t_gpublish::~t_gpublish()
{
    if (viewer != nullptr) delete viewer;
}

void great::t_gpublish::Initialize()
{
    viewer = new gviewer();
    viewer->Show();
}

void great::t_gpublish::UpdateNewState(const IMUState& imu_state)
{
    if (viewer != nullptr)
    {
        if (firstFlag_imu)
        {
            if (imu_state.position.norm() >= 0.0)
            {
                _init_imupos = imu_state.position;
                Eigen::Vector3d BLH = Cart2Geod(_init_imupos, false);
                R_e_n = t_gbase::Cen(BLH).transpose();
                firstFlag_imu = false;
            }
        }

        viewer->SetFrames(vector<pair<Eigen::Matrix3d, Eigen::Vector3d>>(0));
        vector<pair<Eigen::Matrix3d, Eigen::Vector3d>> frames;
        Eigen::Vector3d Position;
        Eigen::Matrix3d atti;

        Position = R_e_n * (imu_state.position - _init_imupos);
        atti = imu_state.orientation.toRotationMatrix();
        frames.push_back(make_pair(atti, Position));

        //att
        viewer->SetFrames(frames);
        //trajectory
        viewer->AddNewPos(Position);
    }

}

void great::t_gpublish::AddMapPoints(const vector<Eigen::Vector3d> & map_points)
{
    vector<Eigen::Vector3d> map_points_e;
    for (const auto& pt : map_points)
    {
        map_points_e.push_back(R_e_n*pt);
    }
    viewer->AddNewPoint(map_points_e);
}


void great::t_gpublish::UpdateVisualPoints(const vector<Eigen::Vector3d> & map_points)
{
    if (viewer != nullptr)
    {
        vector<Eigen::Vector3d> visual_points;
        for (auto point : map_points)
        {
            Eigen::Vector3d vis_point = R_e_n * (point + _init_campos - _init_imupos);
            visual_points.push_back(vis_point);
        }
        viewer->SetPointCloud(visual_points);
    }
}

void great::t_gpublish::UpdatePlanePoints(const vector<Eigen::Vector3d> & pcs, const vector<vector<Eigen::Vector3d>> &_near_points, Eigen::Matrix3d R_l_e, Eigen::Vector3d t_l_e)
{
    if (viewer != nullptr)
    {
        vector<vector<Eigen::Vector3d>> surroundings;
        vector<Eigen::Vector3d> centers;
        for (int i = 0; i < pcs.size(); i++)
        {
            vector<Eigen::Vector3d> one;
            centers.push_back(R_e_n*(R_l_e*pcs.at(i) + t_l_e - _init_imupos));
            for (int j = 0; j < _near_points.at(i).size(); j++)
            {
                one.push_back(R_e_n * (R_l_e*_near_points.at(i).at(j) + t_l_e - _init_imupos));
            }
            surroundings.push_back(one);
        }
        
        viewer->SetPlaneCloud(centers, surroundings);
    }

}
