/**
 * @file         gamb_manager.cpp
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        manager of ambiguity parameters
 * @version      1.0
 * @date         2025-01-01
 *
 * @copyright Copyright (c) 2025, Wuhan University. All rights reserved.
 *
 */

#include "gamb_manager.h"

int gfgomsf::t_gamb_per_ID::endRover()
{
	return _start_rover_count + _rover_list.size() - 1;
}

int gfgomsf::t_gamb_per_ID::startRover()
{
	return _start_rover_count;
}

void gfgomsf::t_gamb_per_ID::setBeg(t_gtime beg)
{
	_beg = beg;
}

void gfgomsf::t_gamb_per_ID::setEnd(t_gtime end)
{
	_end = end;
}

void gfgomsf::t_gamb_per_ID::setSlip()
{
	_is_slip = true;
}

void gfgomsf::t_gamb_per_ID::addRover(int rover_count)
{
	assert(rover_count > _start_rover_count);
	_rover_list.push_back(rover_count);
}

double gfgomsf::t_gamb_per_ID::get_initial_value()
{
	return _initial_amb;
}

double gfgomsf::t_gamb_per_ID::get_est_value()
{
	return _estimated_amb;
}

double gfgomsf::t_gamb_per_ID::is_estimated()
{
	return _is_estimated;
}

void gfgomsf::t_gamb_per_ID::set_est_value(double est_value)
{
	_estimated_amb = est_value;
	_is_estimated = true;
}

void gfgomsf::t_gamb_per_ID::set_fixed_value(double fixed_value)
{
	_fixed_amb = fixed_value;
}

void gfgomsf::t_gamb_per_ID::setFix()
{
	_is_fix = true;
}

string gfgomsf::t_gamb_per_ID::getPRN()
{
	return _sat_prn;
}

pair<FREQ_SEQ, GOBSBAND> gfgomsf::t_gamb_per_ID::getFB()
{
	return _freq_band;
}

gfgomsf::t_gamb_manager::t_gamb_manager()
{
}

gfgomsf::t_gamb_manager::t_gamb_manager(map<GSYS, map<FREQ_SEQ, GOBSBAND>> band_index):
	_band_index(band_index)
{
}

gfgomsf::t_gamb_manager::~t_gamb_manager()
{
}

void gfgomsf::t_gamb_manager::clearState()
{
	_ambiguity.clear();
	_sat_map.clear();
	cur_sats.clear();
	search_index.clear();
	//new add
	//_band_index.clear();
	ambiguity_ids.clear();	
}





void gfgomsf::t_gamb_manager::removeSat(const int &sat_global_id, const int &rover_index)
{

	//delet satellite and the corresponding ambiguities
	auto it = _sat_map.find(sat_global_id);
	if (it != _sat_map.end())
	{
		bool is_new_sat = false;
		shared_ptr<t_gsat_map> sat_map = _sat_map[sat_global_id];
		for (int i = 0; i < sat_map->_amb_ids.size(); i++)
		{
			int id = sat_map->_amb_ids[i];

			auto it_id = _ambiguity.find(id);
			assert(it_id != _ambiguity.end());

			//case 1: new amb-> delet amb and sat
			if (it_id->second->_rover_list.size() == 1 && it_id->second->endRover() == rover_index)
			{
				is_new_sat = true;
				//removed_ambs.push_back(id);
				_ambiguity.erase(id);
				auto it_find = find_if(ambiguity_ids.begin(), ambiguity_ids.end(), [id](const int &amb_id)
				{
					return id == amb_id;
				});
				if (it_find != ambiguity_ids.end())
					ambiguity_ids.erase(it_find);

			}
			else
			{
				//case 2: old amb -> only update the roverlist(as lost tracking)
				if (it_id->second->endRover() == rover_index)
				{
					auto iter_rover = it_id->second->_rover_list.begin();
					for (; iter_rover != it_id->second->_rover_list.end(); iter_rover++)
					{
						if (*iter_rover == rover_index)
						{
							it_id->second->_rover_list.erase(iter_rover);
							break;
						}

					}
				}

			}
		}
		if (is_new_sat)
			_sat_map.erase(sat_global_id);
		else
		{
			it->second->time_span.pop_back();
		}
	}
}



void gfgomsf::t_gamb_manager::slidingWindow()
{
	// Move the sliding window to remove the first frame from the ambiguity set.
	// update amb vector
	// For each rover corresponding to all ambiguities in the ambiguity set, move forward one step.
	for (auto it = _ambiguity.begin(), it_next= _ambiguity.begin(); it != _ambiguity.end(); it = it_next)
	{		
		it_next++;
		if ( it->second->_rover_list[0] == 0)
		{
			if (it->second->_rover_list.size() == 1)
			{
				int id = it->first;
				auto iter_find = find_if(ambiguity_ids.begin(), ambiguity_ids.end(), [id](const int &a)
				{
					return a == id;
				});
				if (iter_find != ambiguity_ids.end())
					ambiguity_ids.erase(iter_find);

				_ambiguity.erase(it);
				continue;
			}
			else
			{
				it->second->_rover_list.erase(it->second->_rover_list.begin());
			}
		}

		for (auto itor = it->second->_rover_list.begin(); itor != it->second->_rover_list.end(); itor++)
		{
			*itor = *itor - 1;
		}

	
		it->second->_start_rover_count = it->second->_rover_list[0];			
	}
	
	//update sat
	for (auto it = _sat_map.begin(), it_next = _sat_map.begin(); it != _sat_map.end(); it = it_next)
	{
		it_next++;
		int amb_count = 0;
		
		// Determine whether the current satellite's ambiguity is still within the ambiguity set of the window.
		for (int i = 0; i < it->second->_amb_ids.size(); i++)
		{
			auto find_it = _ambiguity.find(it->second->_amb_ids[i]);
			if (find_it == _ambiguity.end())
				amb_count++;		

		}
		
		if (amb_count>0)
		{			
			// If one is missing, delete the current satellite.
			cout << "sat: " << it->second->_sat_name << " slide out!!!"<< endl;
			_sat_map.erase(it);
		}
	}
	
}

int gfgomsf::t_gamb_manager::getAmbCount()
{
	int cnt = 0;
	for (auto it :_ambiguity)
	{	
		cnt++;
	}
	return cnt;
}

int gfgomsf::t_gamb_manager::getAmbStartRoverID(const int & amb_search_id)
{	
	for (auto it : _ambiguity)
	{
		if (it.first == amb_search_id)
		{
			return it.second->startRover();
		}
		
	}

	return -1;
}

int gfgomsf::t_gamb_manager::getAmbEndRoverID(const int & amb_search_id)
{
	for (auto it : _ambiguity)
	{
		if (it.first == amb_search_id)
		{
			return it.second->endRover();
		}

	}

	return -1;
}

Eigen::VectorXd gfgomsf::t_gamb_manager::getAmbVector()
{	
	Eigen::VectorXd amb_vec(getAmbCount());
	int amb_index = -1;
	for (auto it : _ambiguity)
	{
		amb_index++;		
		amb_vec(amb_index) = it.second->get_initial_value();
	}
	return amb_vec;	
}

/*
int gfgomsf::t_gamb_manager::get_last_epoch_sats(double time)
{
	cur_sat_prn.clear();
	auto it = _sat_map.begin();
	for (; it != _sat_map.end(); it++)
	{
		if (it->second->time_span.back() == time)
		{
			cur_sats.push_back(make_pair(it->second->_sat_name, it->first));
		}
	}
	return 0;
}
*/

double gfgomsf::t_gamb_manager::getAmb(int amb_id)
{	
	return _ambiguity[amb_id]->get_est_value();
}

double gfgomsf::t_gamb_manager::getInitialAmb(const int & amb_id)
{
	return _ambiguity[amb_id]->get_initial_value();
}

vector<int> gfgomsf::t_gamb_manager::getMarginAmb()
{
	vector<int> margin_amb;
	for (auto it : _ambiguity)
	{
		if(it.second->endRover()==0)
			margin_amb.push_back(it.first);
	}

	for (int i = 0; i < removed_ambs.size(); i++)
	{
		int id = removed_ambs[i];
		auto it = find_if(margin_amb.begin(), margin_amb.end(), [id](const int &a)
		{
			return a == id;
		});
	
		if (it == margin_amb.end())
			margin_amb.push_back(id);
	}
	removed_ambs.clear();
	return margin_amb;
}

vector<int> gfgomsf::t_gamb_manager::getCurWinAmb()
{
	vector<int> cur_amb;
	for (auto it : _ambiguity)
	{
		if(it.second->endRover() == 0)
			continue;
		cur_amb.push_back(it.first);
	}
	return cur_amb;
}


void gfgomsf::t_gamb_manager::updateAmb(int amb_id, double value)
{

	_ambiguity[amb_id]->set_est_value(value);

}


void gfgomsf::t_gamb_manager::generateAmbSearchIndex()
{
	search_index.clear();	
	for (auto it : _ambiguity)
	{		
		pair<int, FREQ_SEQ> sat_freq = make_pair(it.second->_sat_global_id, it.second->getFB().first);
		auto iter = search_index.find(sat_freq);
		assert(iter == search_index.end());
		search_index[sat_freq] = it.first;
	}	
}


void gfgomsf::t_gamb_manager::setEstAmb(const Eigen::VectorXd & x)
{
	int amb_index = -1;
	for (auto &it : _ambiguity)
	{
		it.second->set_est_value(x(++amb_index));
	}
}

int gfgomsf::t_gamb_manager::getAmbSearchIndex(const pair<int, FREQ_SEQ>& sat_freq)
{	
	auto iter = search_index.find(sat_freq);
	if (iter != search_index.end())
		return search_index[sat_freq];
	else
		return -1;
}

void gfgomsf::t_gamb_manager::addNewSat(const t_gtime & cur_time, const int & rover_index, const int & sat_index, int & amb_index, const t_gsatdata & sat_data, t_gallpar params)
{
	string sat_name = sat_data.sat();
	GSYS gnss_system = sat_data.gsys();
	t_gtime obs_time = sat_data.epoch();
	vector<t_gpar> amb_para_per_sat;
	for (unsigned int i = 0; i < params.parNumber(); i++)
	{
		if (params[i].parType == par_type::AMB_IF ||
			params[i].parType == par_type::AMB_L1 ||
			params[i].parType == par_type::AMB_L2 ||
			params[i].parType == par_type::AMB_L3 ||
			params[i].parType == par_type::AMB_L4 ||
			params[i].parType == par_type::AMB_L5)
		{
			if (params[i].prn == sat_name)
			{
				amb_para_per_sat.push_back(params[i]);

			}
		}
	}
	shared_ptr<t_gsat_map> cur_sat(new t_gsat_map(gnss_system, sat_name, obs_time.sow(), sat_index));
	_sat_map[sat_index] = cur_sat;
	addAmb(cur_time, amb_para_per_sat, gnss_system, rover_index, sat_index, amb_index);
}




void gfgomsf::t_gamb_manager::addRover(double time, string sat_name, const int & rover_index)
{
	auto iter = find_if(_sat_map.rbegin(), _sat_map.rend(), [sat_name](pair<int, shared_ptr<t_gsat_map>> it)
	{
		return it.second->_sat_name == sat_name;

	});
	assert(iter != _sat_map.end());
	iter->second->time_span.push_back(time);
	int old_sat_index = iter->second->_global_id;
	vector<int> amb_ids = _sat_map[old_sat_index]->_amb_ids;
	for (int i = 0; i < amb_ids.size(); i++)
	{
		_ambiguity[amb_ids[i]]->addRover(rover_index);
	}
}

void gfgomsf::t_gamb_manager::addRover( string sat_name, const int & rover_index)
{
	map<int, shared_ptr<t_gsat_map>>::const_reverse_iterator iter = _sat_map.rbegin();
	while (iter != _sat_map.rend())
	{
		if (iter->second->_sat_name == sat_name)
			break;
		++iter;
	}
	/*
	map<int, shared_ptr<t_gsat_map>>::const_reverse_iterator iter = _sat_map.rbegin();
	while (iter != _sat_map.rend())
	{
		if (iter->second->_sat_name == sat_name)
			break;
		++iter;
	}
	*/
	assert(iter != _sat_map.end());
	
	int old_sat_index = iter->second->_global_id;	
	vector<int> amb_ids = _sat_map[old_sat_index]->_amb_ids;
	for (int i = 0; i < amb_ids.size(); i++)
	{
		_ambiguity[amb_ids[i]]->addRover(rover_index);
	}
}

bool gfgomsf::t_gamb_manager::is_sat_tracking(string prn)
{
	string sat_name = prn;
	auto iter_find = find_if(cur_sats.begin(), cur_sats.end(), [sat_name](const pair<string, int> &a)
	{
		return a.first == sat_name;
	});
	if (iter_find != cur_sats.end())
		return true;

	return false;
}

int gfgomsf::t_gamb_manager::get_last_epoch_sats(double time)
{
	cur_sats.clear();
	auto it = _sat_map.begin();
	for (; it != _sat_map.end(); it++)
	{
		if (it->second->time_span.back() == time)
		{
			cur_sats.push_back(make_pair(it->second->_sat_name, it->first));
		}
	}
	return 0;
}

int gfgomsf::t_gamb_manager::get_sat_id(string prn)
{
	string sat_name = prn;
	auto iter_find = find_if(cur_sats.begin(), cur_sats.end(), [sat_name](const pair<string, int> &a)
	{
		return a.first == sat_name;
	});
	if (iter_find != cur_sats.end())
	{
		return iter_find->second;
	}

	return -1;
}

bool gfgomsf::t_gamb_manager::addAmb(const t_gtime & cur_time, vector<t_gpar> amb_para, const GSYS & gnss_system, int rover_index, const int & sat_index, int & amb_index)
{
	map<FREQ_SEQ, GOBSBAND> crt_bands = _band_index[gnss_system];
	string sat_name = amb_para[0].prn;

	if (amb_para.size() > 0)
	{
		for (int i = 0; i < amb_para.size(); i++)
		{
			amb_index++;
			FREQ_SEQ freq = ambtype_list[amb_para[i].parType];
			GOBSBAND band = crt_bands[freq];
			pair<FREQ_SEQ, GOBSBAND> freq_band = make_pair(freq, band);
			double sd_amb_value = amb_para[i].value();
			shared_ptr<t_gamb_per_ID> amb_i(new t_gamb_per_ID(freq_band, gnss_system, sat_name, sat_index, amb_index, rover_index, sd_amb_value));
			amb_i->setBeg(cur_time);
			_ambiguity[amb_index] = amb_i;
			_sat_map[sat_index]->_amb_ids.push_back(amb_index);
			ambiguity_ids.push_back(amb_index);			
		}
	}
	else
	{
		cout << "sat %s has no amb para!!!" << sat_name << endl;
		return false;
	}

	return true;
}

void gfgomsf::t_gamb_manager::addGpara(t_gallpar & params, const int & idx, double value)
{
	for (auto it : _ambiguity)
	{
		if (it.first == idx)
		{
			par_type amb_type = freq_ambtype_list[it.second->getFB().first];
			string site = params[0].site;
			string sat = it.second->getPRN();
			t_gpar newPar(site, amb_type, params.parNumber() + 1, sat);
			newPar.value(value);
			params.addParam(newPar);
		}
	}	
}


void gfgomsf::t_gamb_manager::checkDDSat(const int & sat_id)
{
	auto it = _sat_map.find(sat_id);
	if (it != _sat_map.end())
		it->second->DD_used = true;

}
void gfgomsf::t_gamb_manager::updateAmb()
{
	removed_sat.clear();
	for (auto it : _sat_map)
	{
		if (it.second->DD_used)
		{			
			continue;
		}
			
		removed_sat[it.second->_sat_name] = it.first;
		for (int i = 0; i < it.second->_amb_ids.size(); i++)
		{
			int id = it.second->_amb_ids[i];
			_ambiguity.erase(id);
			auto it_find = find_if(ambiguity_ids.begin(), ambiguity_ids.end(), [id](const int &amb_id)
			{
				return id == amb_id;
			});
			if (it_find != ambiguity_ids.end())
				ambiguity_ids.erase(it_find);
		}
		_sat_map.erase(it.first);
	}

	for (auto &it : _sat_map)
	{
		if (it.second->DD_used = true)
			it.second->DD_used = false;
	}
}

