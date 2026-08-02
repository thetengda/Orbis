#ifndef GFGOESTIMATOR_H
#define GFGOESTIMATOR_H
/**
 * @file         gfgo_estimator.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        state of estimator
 * @version      1.0
 * @date         2025-01-01
 *
 * @copyright Copyright (c) 2025, Wuhan University. All rights reserved.
 *
 */
#include"gall/gallpar.h"
#include "gfgo_para.h"
#include "gfactor/gimu_factor.h"
#include "gfactor/gmarginalization_factor.h"
#include "gfactor/gnss_info_factor.h"
#include "gfactor/ginitial_gnss_amb_factor.h"

using namespace gnut;
using namespace great;

namespace gfgo
{   
	/**
	*@brief  Class for sliding window based factor graph optimization
	*/
	class LibGREAT_LIBRARY_EXPORT  t_gfgo : virtual public t_gfgo_para
	{
	public:
		t_gfgo(t_gsetbase* set);
		~t_gfgo();
		void clear_state();			
				
	protected:
		virtual void _double_to_vector(void) = 0;					/// data storage type conversion(double-->vector)
		virtual void _vector_to_double(void) = 0;					/// data storage type conversion(vector-->double)
	
	//for SLAM
	protected:			
		MarginalizationInfo *_last_marginalization_info = nullptr;			/// info of last marginalized frame
		vector<double *> _last_marginalization_parameter_blocks;	/// block to store info of last marginalized parameters
		/// parameters  blocks for Ceres solver 
		double _para_pose[WINDOW_SIZE + 1][SIZE_POSE];				///< pose 
		double _para_speed_bias[WINDOW_SIZE + 1][SIZE_SPEEDBIAS];	///< speed and bias
		double _para_feature[NUM_OF_F][SIZE_FEATURE];				///< feature inverse depth
		double _para_ex_pose[2][SIZE_POSE];							///< external pose 
		double _para_retrive_pose[SIZE_POSE];						///< unused
		double _para_td[1][1];										///< time delay
		double _para_tr[1][1];										///< unused    
	protected:
		vector<double*> t_array; ///<for global pose graph optimization,the translation variables are stored in it
		vector<double*> q_array; ///<for global pose graph optimization,the rotation variables are stored in it
    
	//for GNSS
	protected:
		double _para_CRD[GWINDOW_SIZE][SIZE_CRD];                     ///<ECEF coordinate (XYZ) in sliding window
		double _para_SPEED[GWINDOW_SIZE][SIZE_SPEED];                 ///<ECEF Velocity (Vx,Vy,Vz) in sliding window
		double _para_amb[NUM_OF_ARC][SIZE_AMB];                    ///carrier-phase ambiguity arc tracked in sliding window 



		//for PPP IF:
		double _para_CLK[GWINDOW_SIZE + 1][SIZE_CLK];                    ///Receiver clock error in sliding window
		double _para_TRP[GWINDOW_SIZE + 1][SIZE_TRP];                    ///Tropospheric delay in sliding window (random walk)
		double _para_ISB_GAL[GWINDOW_SIZE + 1][SIZE_ISB];				      ///GAL ISB in sliding window (random walk)
		double _para_ISB_BDS[GWINDOW_SIZE + 1][SIZE_ISB];				      ///GAL ISB in sliding window (random walk)
		double _para_ISB_GLO[GWINDOW_SIZE + 1][SIZE_ISB];				      ///GLO ISB in sliding window (random walk)
		double _para_AMB_IF[NUM_OF_ARC][SIZE_AMB];					  ///IF carrier-phase ambiguity arc tracked in sliding window
	};
}

#endif
