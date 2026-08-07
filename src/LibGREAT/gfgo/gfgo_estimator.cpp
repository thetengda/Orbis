/**
 * @file         gfgo_estimator.cpp
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        state of estimator
 * @version      1.0
 * @date         2025-01-01
 *
 * @copyright Copyright (c) 2025, Wuhan University. All rights reserved.
 *
 */

#include"gfgo_estimator.h"

#include <cstring>

namespace gfgo
{
	t_gfgo::t_gfgo(t_gsetbase * set):
	t_gfgo_para(set)
	{
		clear_state();


	}

	t_gfgo::~t_gfgo()
	{		
		 delete _last_marginalization_info;
	}

	void t_gfgo::clear_state()
	{	
		if (_last_marginalization_info != nullptr) delete _last_marginalization_info;		
		_last_marginalization_info = nullptr;		
		_last_marginalization_parameter_blocks.clear();
		std::memset(_para_ISB_QZS, 0, sizeof(_para_ISB_QZS));
		std::memset(_para_AMB_RAW, 0, sizeof(_para_AMB_RAW));
		std::memset(_para_SION, 0, sizeof(_para_SION));
	}
	
}
