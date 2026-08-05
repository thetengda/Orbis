#include "gambIF_manager.h"
gfgomsf::t_gambIF_per_ID::t_gambIF_per_ID(pair<FREQ_SEQ, GOBSBAND> freq_band1, pair<FREQ_SEQ, GOBSBAND> freq_band2, GSYS sys, string prn, int sat_id, int amb_id, int start_rover, double initial_amb) :
	t_gamb_per_ID(sys, prn, sat_id, amb_id, start_rover, initial_amb)
{
	_freq_band1 = freq_band1;
	_freq_band2 = freq_band2;
}

gfgomsf::t_gambIF_manager::t_gambIF_manager()
{
}

gfgomsf::t_gambIF_manager::t_gambIF_manager(map<GSYS, map<FREQ_SEQ, GOBSBAND>> band_index) :
	t_gamb_manager(band_index)
{
}

gfgomsf::t_gambIF_manager::~t_gambIF_manager()
{
}

void gfgomsf::t_gambIF_manager::clearState()
{
	t_gamb_manager::clearState();
	_ambiguityIF.clear();
	search_index_IF.clear();
}

bool gfgomsf::t_gambIF_manager::addNewSat(const t_gtime & cur_time, const int &rover_index, const int &sat_index, int &amb_index, const t_gsatdata &sat_data, t_gallpar params)
{
	string sat_name = sat_data.sat();
	GSYS gnss_system = sat_data.gsys();
	t_gtime obs_time = sat_data.epoch();//add for MultiWindow
	vector<t_gpar> amb_para_per_sat;
	for (unsigned int i = 0; i < params.parNumber(); i++)
	{
		if (params[i].parType == par_type::AMB_IF)
		{
			if (params[i].prn == sat_name)
			{
				amb_para_per_sat.push_back(params[i]);
			}
		}
	}
	//shared_ptr<t_gsat_map> cur_sat(new t_gsat_map(gnss_system, sat_name, sat_index));//delete for MultiWindow
	shared_ptr<t_gsat_map> cur_sat(new t_gsat_map(gnss_system, sat_name, obs_time.sow(), sat_index));//add for MultiWindow
	_sat_map[sat_index] = cur_sat;
	const bool added = addAmb(cur_time, amb_para_per_sat, gnss_system,
		rover_index, sat_index, amb_index);
	if (!added)
		_sat_map.erase(sat_index);
	return added;
}

// bool gfgomsf::t_gambIF_manager::addAmb(const t_gtime & cur_time, vector<t_gpar> amb_para, const GSYS & gnss_system, int rover_index, const int & sat_index, int & amb_index)
// {
// 	map<FREQ_SEQ, GOBSBAND> crt_bands = _band_index[gnss_system];
// 	string sat_name = amb_para[0].prn;
// 	if (amb_para.size() > 0)
// 	{
// 		for (int i = 0; i < amb_para.size(); i++)
// 		{
// 			amb_index++;
// 			pair<FREQ_SEQ, GOBSBAND> freq_band1 = make_pair(FREQ_1, crt_bands[FREQ_1]);
// 			pair<FREQ_SEQ, GOBSBAND> freq_band2 = make_pair(FREQ_2, crt_bands[FREQ_2]);
// 			double if_amb_value = amb_para[i].value();
// 			shared_ptr<t_gambIF_per_ID> amb_i(new t_gambIF_per_ID(freq_band1, freq_band2, gnss_system, sat_name, sat_index, amb_index, rover_index, if_amb_value));
// 			amb_i->setBeg(cur_time);
// 			_ambiguityIF[amb_index] = amb_i;
// 			_sat_map[sat_index]->_amb_ids.push_back(amb_index);
// 			ambiguity_ids.push_back(amb_index);
// 		}
// 	}
// 	else
// 	{
// 		cout << "sat %s has no amb para!!!" << sat_name << endl;
// 		return false;
// 	}
//
// 	return true;
// }

bool gfgomsf::t_gambIF_manager::addAmb(const t_gtime& cur_time,
									   vector<t_gpar> amb_para,
									   const GSYS& gnss_system,
									   int rover_index,
									   const int& sat_index,
									   int& amb_index)
{
	if (amb_para.empty())
	{
		cout << "sat id " << sat_index << " has no amb para!!!" << endl;
		return false;
	}
	auto sat_it = _sat_map.find(sat_index);
	if (sat_it == _sat_map.end() || !sat_it->second)
		return false;
	auto bands_it = _band_index.find(gnss_system);
	if (bands_it == _band_index.end())
		return false;
	const map<FREQ_SEQ, GOBSBAND> &crt_bands = bands_it->second;
	auto f1 = crt_bands.find(FREQ_1);
	auto f2 = crt_bands.find(FREQ_2);
	if (f1 == crt_bands.end() || f2 == crt_bands.end() ||
		f1->second == BAND || f2->second == BAND)
		return false;

	string sat_name = amb_para[0].prn;

	for (int i = 0; i < amb_para.size(); i++)
	{
		amb_index++;
		pair<FREQ_SEQ, GOBSBAND> freq_band1 = make_pair(FREQ_1, f1->second);
		pair<FREQ_SEQ, GOBSBAND> freq_band2 = make_pair(FREQ_2, f2->second);
		double if_amb_value = amb_para[i].value();

		shared_ptr<t_gambIF_per_ID> amb_i(
			new t_gambIF_per_ID(freq_band1, freq_band2, gnss_system,
								sat_name, sat_index, amb_index,
								rover_index, if_amb_value));

		amb_i->setBeg(cur_time);
		amb_i->setEnd(LAST_TIME);
		_ambiguityIF[amb_index] = amb_i;
		sat_it->second->_amb_ids.push_back(amb_index);
		ambiguity_ids.push_back(amb_index);
	}

	return true;
}

void gfgomsf::t_gambIF_manager::addRover(string sat_name, const int &rover_index)
{
	auto iter = find_if(_sat_map.begin(), _sat_map.end(), [sat_name](pair<int, shared_ptr<t_gsat_map>> it)
	{
		return it.second->_sat_name == sat_name;

	});
	assert(iter != _sat_map.end());

	int old_sat_index = iter->second->_global_id;
	vector<int> amb_ids = _sat_map[old_sat_index]->_amb_ids;
	for (int i = 0; i < amb_ids.size(); i++)
	{
		_ambiguityIF[amb_ids[i]]->addRover(rover_index);
	}
}

//add for MultiWindow
void gfgomsf::t_gambIF_manager::addRover(double time, string sat_name, const int &rover_index)
{
	/*auto iter = find_if(_sat_map.begin(), _sat_map.end(), [sat_name](pair<int, shared_ptr<t_gsat_map>> it)
	{
		return it.second->_sat_name == sat_name;

	});*/
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
		_ambiguityIF[amb_ids[i]]->addRover(rover_index);
	}
}

void gfgomsf::t_gambIF_manager::updateAmb(int amb_id, double value)
{
	_ambiguityIF[amb_id]->set_est_value(value);
}

double gfgomsf::t_gambIF_manager::getAmb(int amb_id)
{
	return _ambiguityIF[amb_id]->get_est_value();
}

bool gfgomsf::t_gambIF_manager::isEstimated(const int & amb_id)
{
	return _ambiguityIF[amb_id]->is_estimated();
}

double gfgomsf::t_gambIF_manager::getInitialAmb(const int & amb_id)
{
	return _ambiguityIF[amb_id]->get_initial_value();
}

void gfgomsf::t_gambIF_manager::generateAmbSearchIndex()
{
	search_index_IF.clear();
	for (auto it : _ambiguityIF)
	{
		auto iter = search_index_IF.find(it.second->_sat_global_id);
		assert(iter == search_index_IF.end());
		search_index_IF[it.second->_sat_global_id] = it.first;
	}
}

int gfgomsf::t_gambIF_manager::getAmbSearchIndex(const int& sat_id)
{
	auto iter = search_index_IF.find(sat_id);
	if (iter != search_index_IF.end())
		return search_index_IF[sat_id];
	else
		return -1;
}



void gfgomsf::t_gambIF_manager::removeSat(const int &sat_global_id, const int &rover_index)
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

			auto it_id = _ambiguityIF.find(id);
			assert(it_id != _ambiguityIF.end());

			//case 1: new amb-> delet amb and sat
			if (it_id->second->_rover_list.size() == 1 && it_id->second->endRover() == rover_index)
			{
				is_new_sat = true;
				_ambiguityIF.erase(id);
				auto it_find = find_if(ambiguity_ids.begin(), ambiguity_ids.end(), [id](const int &amb_id)
				{
					return id == amb_id;
				});
				if (it_find != ambiguity_ids.end())
					ambiguity_ids.erase(it_find);
			}
			else//case 2: old amb -> only update the roverlist(as lost tracking)
			{
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
				}//if (it_id->second->endRover() == rover_index)
			}//else
		}//for (int i = 0; i < sat_map->_amb_ids.size(); i++)
		if (is_new_sat)
			_sat_map.erase(sat_global_id);
		else if (!it->second->time_span.empty())//add for MultiWindow
		{
			it->second->time_span.pop_back();
		}
	}//if (it != _sat_map.end())
}

void gfgomsf::t_gambIF_manager::addGpara(t_gallpar & params, const int & idx, double value)
{
	if (params.parNumber() <= 0)
		return;

	for (auto it : _ambiguityIF)
	{
		if (it.first == idx)
		{
			par_type amb_type = par_type::AMB_IF;
			string site = params[0].site;
			string sat = it.second->getPRN();
			t_gpar newPar(site, amb_type, params.parNumber() + 1, sat);
			newPar.value(value);
			params.addParam(newPar);
		}
	}
}

vector<int> gfgomsf::t_gambIF_manager::getMarginAmb()
{
	vector<int> margin_amb;
	for (auto it : _ambiguityIF)
	{
		if (it.second->endRover() == 0)
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

int gfgomsf::t_gambIF_manager::getAmbStartRoverID(const int & amb_search_id)
{
	for (auto it : _ambiguityIF)
	{
		if (it.first == amb_search_id)
		{
			return it.second->startRover();
		}
	}
	return -1;
}

int gfgomsf::t_gambIF_manager::getAmbEndRoverID(const int & amb_search_id)
{
	for (auto it : _ambiguityIF)
	{
		if (it.first == amb_search_id)
		{
			return it.second->endRover();
		}

	}
	return -1;
}

vector<int> gfgomsf::t_gambIF_manager::getCurWinAmb()
{
	vector<int> cur_amb;
	for (auto it : _ambiguityIF)
	{
		if (it.second->endRover() == 0)
			continue;

		cur_amb.push_back(it.first);
	}
	return cur_amb;
}

vector<int> gfgomsf::t_gambIF_manager::getCurWinAmb(int window_type)
{
	assert(window_type == 0 || window_type == 1);
	vector<int> cur_amb;

	//case 0: normal GNSS  window
	if (window_type == 0)
	{
		for (auto it : _ambiguityIF)
		{
			if (it.second->endRover() == 0)
				continue;

			cur_amb.push_back(it.first);
		}
	}
	//case 1: gnss/odo fusion case
	else
	{
		for (auto it : _ambiguityIF)
		{
			cur_amb.push_back(it.first);
		}
	}
	return cur_amb;
}

void gfgomsf::t_gambIF_manager::slidingWindow()
{
	//update amb vector
	for (auto it = _ambiguityIF.begin(), it_next = _ambiguityIF.begin(); it != _ambiguityIF.end(); it = it_next)
	{
		it_next++;
		if (it->second->_rover_list[0] == 0)
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

				_ambiguityIF.erase(it);
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

	}//_ambiguityIF

	//update sat
	for (auto it = _sat_map.begin(), it_next = _sat_map.begin(); it != _sat_map.end(); it = it_next)
	{
		it_next++;
		int amb_count = 0;

		for (int i = 0; i < it->second->_amb_ids.size(); i++)
		{
			auto find_it = _ambiguityIF.find(it->second->_amb_ids[i]);
			if (find_it == _ambiguityIF.end())
				amb_count++;

		}

		if (amb_count > 0)
		{
			//cout << "sat: " << it->second->_sat_name << " slide out!!!" << endl;
			_sat_map.erase(it);
		}
	}//_sat_map
}

void gfgomsf::t_gambIF_manager::slideNew(const int & rover_number)
{
	//update amb vector
	for (auto it = _ambiguityIF.begin(), it_next = _ambiguityIF.begin(); it != _ambiguityIF.end(); it = it_next)
	{
		it_next++;

		auto iter_rover = it->second->_rover_list.begin();
		for (; iter_rover != it->second->_rover_list.end(); iter_rover++)
		{
			if (*iter_rover == rover_number)
			{
				it->second->_rover_list.erase(iter_rover);
				break;
			}

		}

		if (it->second->_rover_list.empty())
		{
			int id = it->first;
			auto iter_find = find_if(ambiguity_ids.begin(), ambiguity_ids.end(), [id](const int &a)
			{
				return a == id;
			});
			if (iter_find != ambiguity_ids.end())
				ambiguity_ids.erase(iter_find);

			_ambiguityIF.erase(it);
			continue;
		}

	}

	//update sat
	for (auto it = _sat_map.begin(), it_next = _sat_map.begin(); it != _sat_map.end(); it = it_next)
	{
		it_next++;
		int amb_count = 0;

		for (int i = 0; i < it->second->_amb_ids.size(); i++)
		{
			auto find_it = _ambiguityIF.find(it->second->_amb_ids[i]);
			if (find_it == _ambiguityIF.end())
				amb_count++;

		}

		if (amb_count > 0)
		{
			cout << "sat: " << it->second->_sat_name << " slide out!!!" << endl;
			_sat_map.erase(it);
		}

	}

}

void gfgomsf::t_gambIF_manager::stopAmbPropagating(const int& rover_number)
{
	for (auto it = _ambiguityIF.begin(), it_next = _ambiguityIF.begin(); it != _ambiguityIF.end(); it = it_next)
	{
		it_next++;

		auto iter_rover = it->second->_rover_list.begin();
		for (; iter_rover != it->second->_rover_list.end(); iter_rover++)
		{
			if (*iter_rover == rover_number)
			{
				it->second->_rover_list.erase(iter_rover);
				break;
			}
		}
	}
}
