#include "gambRAW_manager.h"
#include "gfgo/gfgo_para.h"

#include <algorithm>
#include <cmath>
#include <sstream>

using namespace std;
using namespace gnut;

namespace gfgomsf
{
    namespace
    {
        bool rawCarrierAvailable(
            const t_gsatdata &sat_data,
            const map<GSYS, map<FREQ_SEQ, GOBSBAND>> &band_index,
            FREQ_SEQ freq)
        {
            auto system_it = band_index.find(sat_data.gsys());
            if (system_it == band_index.end())
                return false;
            auto band_it = system_it->second.find(freq);
            if (band_it == system_it->second.end() || band_it->second == BAND)
                return false;

            const t_gobs phase(sat_data.select_phase(band_it->second, true));
            const double value = sat_data.obs_L(phase);
            return phase.gobs() != GOBS::X &&
                   std::isfinite(value) && std::fabs(value) > 0.0;
        }

        bool rawAmbiguityFrequency(par_type type, FREQ_SEQ &freq)
        {
            switch (type)
            {
            case par_type::AMB_L1: freq = FREQ_1; return true;
            case par_type::AMB_L2: freq = FREQ_2; return true;
            case par_type::AMB_L3: freq = FREQ_3; return true;
            case par_type::AMB_L4: freq = FREQ_4; return true;
            case par_type::AMB_L5: freq = FREQ_5; return true;
            default: return false;
            }
        }
    }

    t_gambRAW_manager::t_gambRAW_manager() : t_gamb_manager()
    {
    }

    t_gambRAW_manager::t_gambRAW_manager(
        const map<GSYS, map<FREQ_SEQ, GOBSBAND>> &band_index)
        : t_gamb_manager(band_index)
    {
    }

    t_gambRAW_manager::~t_gambRAW_manager() = default;

    bool t_gambRAW_manager::addNewSat(const t_gtime &cur_time,
                                       const int &rover_index,
                                       const int &sat_index,
                                       int &amb_index,
                                       const t_gsatdata &sat_data,
                                       t_gallpar params)
    {
        const string sat_name = sat_data.sat();
        const GSYS gnss_system = sat_data.gsys();
        auto sat_it = _sat_map.find(sat_index);
        if (sat_it == _sat_map.end())
        {
            _sat_map[sat_index] = shared_ptr<t_gsat_map>(
                new t_gsat_map(gnss_system, sat_name,
                               sat_data.epoch().sow(), sat_index));
            sat_it = _sat_map.find(sat_index);
        }

        vector<t_gpar> ambiguity_parameters;
        for (unsigned int i = 0; i < params.parNumber(); ++i)
        {
            FREQ_SEQ freq = FREQ_X;
            if (params[i].prn == sat_name &&
                rawAmbiguityFrequency(params[i].parType, freq) &&
                rawCarrierAvailable(sat_data, _band_index, freq))
                ambiguity_parameters.push_back(params[i]);
        }

        // A code-only satellite remains in _sat_map so that its SION and
        // code equations can be used. It simply has no ambiguity entries.
        return addAmb(cur_time, ambiguity_parameters, gnss_system,
                      rover_index, sat_index, amb_index);
    }

    bool t_gambRAW_manager::addAmb(const t_gtime &cur_time,
                                    vector<t_gpar> amb_para,
                                    const GSYS &gnss_system,
                                    int rover_index,
                                    const int &sat_index,
                                    int &amb_index)
    {
        if (amb_para.empty())
            return false;

        auto sat_it = _sat_map.find(sat_index);
        if (sat_it == _sat_map.end() || !sat_it->second)
            return false;

        auto bands_it = _band_index.find(gnss_system);
        if (bands_it == _band_index.end())
            return false;

        bool added = false;
        for (const auto &parameter : amb_para)
        {
            FREQ_SEQ freq = FREQ_X;
            if (!rawAmbiguityFrequency(parameter.parType, freq))
                continue;

            auto band_it = bands_it->second.find(freq);
            if (band_it == bands_it->second.end() || band_it->second == BAND)
                continue;

            // Do not create two live arcs for the same satellite/frequency.
            if (activeAmbiguity(sat_index, freq) >= 0)
                continue;

			// Parameter storage in t_gfgo is currently a fixed-size array.  Never
			// publish an arc whose ID cannot be represented by that storage; the
			// caller can then reject the observation explicitly instead of silently
			// building a code-only graph with an out-of-range ambiguity address.
			if (amb_index + 1 >= NUM_OF_ARC)
				return added;

            ++amb_index;
            shared_ptr<t_gamb_per_ID> ambiguity(new t_gamb_per_ID(
                make_pair(freq, band_it->second), gnss_system,
                parameter.prn, sat_index, amb_index, rover_index,
                parameter.value()));
            ambiguity->setBeg(cur_time);
            ambiguity->setEnd(LAST_TIME);
            _ambiguity[amb_index] = ambiguity;
            sat_it->second->_amb_ids.push_back(amb_index);
            ambiguity_ids.push_back(amb_index);
            added = true;
        }
        return added;
    }

    int t_gambRAW_manager::activeAmbiguity(int sat_index, FREQ_SEQ freq) const
    {
        auto sat_it = _sat_map.find(sat_index);
        if (sat_it == _sat_map.end() || !sat_it->second)
            return -1;

        for (auto id_it = sat_it->second->_amb_ids.rbegin();
             id_it != sat_it->second->_amb_ids.rend(); ++id_it)
        {
            auto ambiguity_it = _ambiguity.find(*id_it);
            if (ambiguity_it != _ambiguity.end() &&
                ambiguity_it->second->getFB().first == freq &&
                ambiguity_it->second->_end == LAST_TIME)
                return *id_it;
        }
        return -1;
    }

    void t_gambRAW_manager::eraseAmbiguity(int amb_id)
    {
        _ambiguity.erase(amb_id);
        auto id_it = find(ambiguity_ids.begin(), ambiguity_ids.end(), amb_id);
        if (id_it != ambiguity_ids.end())
            ambiguity_ids.erase(id_it);
    }

    bool t_gambRAW_manager::resetFrequencyArc(const t_gtime &cur_time,
                                               int rover_index,
                                               int sat_index,
                                               FREQ_SEQ freq,
                                               int &amb_index,
                                               t_gallpar params)
    {
        auto sat_it = _sat_map.find(sat_index);
        if (sat_it == _sat_map.end() || !sat_it->second)
            return false;

        auto bands_it = _band_index.find(sat_it->second->_gnss_system);
        if (bands_it == _band_index.end())
            return false;
        auto band_it = bands_it->second.find(freq);
        if (band_it == bands_it->second.end() || band_it->second == BAND)
            return false;

        t_gpar new_parameter;
        bool has_new_parameter = false;
        for (unsigned int i = 0; i < params.parNumber(); ++i)
        {
            FREQ_SEQ parameter_freq = FREQ_X;
            if (params[i].prn == sat_it->second->_sat_name &&
                rawAmbiguityFrequency(params[i].parType, parameter_freq) &&
                parameter_freq == freq)
            {
                new_parameter = params[i];
                has_new_parameter = true;
                break;
            }
        }

        // Validate the replacement arc before retiring the old one.  This
        // keeps a cycle-slip reset atomic when the parameter list or band map
        // is incomplete.
        if (!has_new_parameter)
            return false;

        const int old_id = activeAmbiguity(sat_index, freq);

		if (amb_index + 1 >= NUM_OF_ARC)
			return false;

        ++amb_index;
        const int new_id = amb_index;
        shared_ptr<t_gamb_per_ID> new_ambiguity(new t_gamb_per_ID(
            make_pair(freq, band_it->second), sat_it->second->_gnss_system,
            new_parameter.prn, sat_index, new_id, rover_index,
            new_parameter.value()));
        new_ambiguity->setBeg(cur_time);
        new_ambiguity->setEnd(LAST_TIME);
        _ambiguity[new_id] = new_ambiguity;
        sat_it->second->_amb_ids.push_back(new_id);
        ambiguity_ids.push_back(new_id);

        // The replacement exists, so now retire only the old arc for this
        // frequency.  Other frequencies and the stable satellite/SION ID are
        // left untouched.
        if (old_id >= 0)
        {
            auto old_it = _ambiguity.find(old_id);
            if (old_it != _ambiguity.end())
            {
                auto &rovers = old_it->second->_rover_list;
                for (auto rover_it = rovers.begin(); rover_it != rovers.end();)
                {
                    if (*rover_it == rover_index)
                        rover_it = rovers.erase(rover_it);
                    else
                        ++rover_it;
                }
                if (rovers.empty())
                {
                    auto amb_id_it = find(sat_it->second->_amb_ids.begin(),
                                          sat_it->second->_amb_ids.end(), old_id);
                    if (amb_id_it != sat_it->second->_amb_ids.end())
                        sat_it->second->_amb_ids.erase(amb_id_it);
                    eraseAmbiguity(old_id);
                }
                else
                {
                    old_it->second->setEnd(cur_time);
                }
            }
        }
        return true;
    }

    void t_gambRAW_manager::addRoverToLatest(const string &sat_name,
                                             int rover_index, double time,
                                             bool append_time)
    {
        auto sat_it = _sat_map.rbegin();
        for (; sat_it != _sat_map.rend(); ++sat_it)
        {
            if (sat_it->second && sat_it->second->_sat_name == sat_name)
                break;
        }
        if (sat_it == _sat_map.rend() || !sat_it->second)
            return;

        if (append_time &&
            (sat_it->second->time_span.empty() ||
             sat_it->second->time_span.back() != time))
            sat_it->second->time_span.push_back(time);

        for (const int amb_id : sat_it->second->_amb_ids)
        {
            auto ambiguity_it = _ambiguity.find(amb_id);
            if (ambiguity_it == _ambiguity.end() ||
                ambiguity_it->second->_end != LAST_TIME)
                continue;
            if (ambiguity_it->second->_rover_list.empty() ||
                ambiguity_it->second->endRover() < rover_index)
                ambiguity_it->second->addRover(rover_index);
        }
    }

    void t_gambRAW_manager::addRover(string sat_name,
                                     const int &rover_index)
    {
        addRoverToLatest(sat_name, rover_index, 0.0, false);
    }

    void t_gambRAW_manager::addRover(double time, string sat_name,
                                     const int &rover_index)
    {
        addRoverToLatest(sat_name, rover_index, time, true);
    }

    void t_gambRAW_manager::removeSat(const int &sat_global_id,
                                      const int &rover_index)
    {
        auto sat_it = _sat_map.find(sat_global_id);
        if (sat_it == _sat_map.end() || !sat_it->second)
            return;

        vector<int> remaining_ambiguities;
        for (const int amb_id : sat_it->second->_amb_ids)
        {
            auto ambiguity_it = _ambiguity.find(amb_id);
            if (ambiguity_it == _ambiguity.end())
                continue;

            auto &rovers = ambiguity_it->second->_rover_list;
            bool removed_current = false;
            for (auto rover_it = rovers.begin(); rover_it != rovers.end();)
            {
                if (*rover_it == rover_index)
                {
                    rover_it = rovers.erase(rover_it);
                    removed_current = true;
                }
                else
                    ++rover_it;
            }

            if (removed_current && rovers.empty())
                eraseAmbiguity(amb_id);
            else
                remaining_ambiguities.push_back(amb_id);
        }
        sat_it->second->_amb_ids.swap(remaining_ambiguities);

        if (!sat_it->second->time_span.empty())
            sat_it->second->time_span.pop_back();
        if (sat_it->second->time_span.empty() &&
            sat_it->second->_amb_ids.empty())
            _sat_map.erase(sat_it);
    }

    void t_gambRAW_manager::slidingWindow(double outgoing_time)
    {
        for (auto it = _ambiguity.begin(); it != _ambiguity.end();)
        {
            auto current = it++;
            auto &rovers = current->second->_rover_list;
            if (rovers.empty())
            {
                eraseAmbiguity(current->first);
                continue;
            }
            if (rovers.front() == 0)
            {
                if (rovers.size() == 1)
                {
                    eraseAmbiguity(current->first);
                    continue;
                }
                rovers.erase(rovers.begin());
            }
            for (int &rover : rovers)
                --rover;
            current->second->_start_rover_count = rovers.front();
        }

        for (auto sat_it = _sat_map.begin(); sat_it != _sat_map.end();)
        {
            auto outgoing_it = std::find(sat_it->second->time_span.begin(),
                                         sat_it->second->time_span.end(),
                                         outgoing_time);
            if (outgoing_it != sat_it->second->time_span.end())
                sat_it->second->time_span.erase(outgoing_it);

            vector<int> remaining;
            for (const int amb_id : sat_it->second->_amb_ids)
            {
                if (_ambiguity.find(amb_id) != _ambiguity.end())
                    remaining.push_back(amb_id);
            }
            sat_it->second->_amb_ids.swap(remaining);
            if (sat_it->second->_amb_ids.empty() &&
                sat_it->second->time_span.empty())
                sat_it = _sat_map.erase(sat_it);
            else
                ++sat_it;
        }
    }

    void t_gambRAW_manager::generateAmbSearchIndex()
    {
        search_index.clear();
        for (const auto &ambiguity_pair : _ambiguity)
        {
            const int amb_id = ambiguity_pair.first;
            const auto &ambiguity = ambiguity_pair.second;
            // The search index is used for the current epoch.  Historical
            // arcs remain in _ambiguity for window factors and marginalization
            // but must never be selected for a new phase equation.
            if (ambiguity->_end != LAST_TIME)
                continue;
            const pair<int, FREQ_SEQ> key(ambiguity->_sat_global_id,
                                          ambiguity->getFB().first);
            auto current = search_index.find(key);
            if (current == search_index.end())
            {
                search_index[key] = amb_id;
                continue;
            }

            auto current_ambiguity = _ambiguity.find(current->second);
            if (current_ambiguity == _ambiguity.end() ||
                amb_id > current->second)
                current->second = amb_id;
        }
    }

    int t_gambRAW_manager::getLatestSatId(const string &sat_name) const
    {
        for (auto sat_it = _sat_map.rbegin(); sat_it != _sat_map.rend(); ++sat_it)
        {
            if (sat_it->second && sat_it->second->_sat_name == sat_name)
                return sat_it->first;
        }
        return -1;
    }

    bool t_gambRAW_manager::hasActiveArc(int sat_index, FREQ_SEQ freq) const
    {
        return activeAmbiguity(sat_index, freq) >= 0;
    }

    bool t_gambRAW_manager::getArcTime(int amb_id, t_gtime &beg,
                                       t_gtime &end) const
    {
        auto ambiguity_it = _ambiguity.find(amb_id);
        if (ambiguity_it == _ambiguity.end() || !ambiguity_it->second)
            return false;
        beg = ambiguity_it->second->_beg;
        end = ambiguity_it->second->_end;
        return true;
    }

    bool t_gambRAW_manager::containsArc(int amb_id) const noexcept
    {
        const auto ambiguity_it = _ambiguity.find(amb_id);
        return ambiguity_it != _ambiguity.end() && ambiguity_it->second != nullptr;
    }

    bool t_gambRAW_manager::getArcInfo(int amb_id, RawArcInfo &info) const
    {
        const auto ambiguity_it = _ambiguity.find(amb_id);
        if (ambiguity_it == _ambiguity.end() || !ambiguity_it->second)
            return false;

        const auto &ambiguity = ambiguity_it->second;
        const FREQ_SEQ frequency = ambiguity->getFB().first;
        info.arc_id = amb_id;
        info.sat_id = ambiguity->_sat_global_id;
        info.start_node = ambiguity->_start_rover_count;
        info.end_node = ambiguity->_rover_list.empty() ? -1 : ambiguity->endRover();
        info.sat = ambiguity->_sat_prn;
        info.stable_key = ambiguityKey(info.sat, frequency, amb_id);
        info.freq = frequency;
        info.beg = ambiguity->_beg;
        info.end = ambiguity->_end;
        const auto search = search_index.find(make_pair(info.sat_id, frequency));
        info.current_search_arc = search != search_index.end() && search->second == amb_id;
        return true;
    }

    string t_gambRAW_manager::ambiguityKey(const string &sat, FREQ_SEQ freq, int arc)
    {
        ostringstream os;
        os << sat << ":F" << static_cast<int>(freq) << ":A" << arc;
        return os.str();
    }

    string t_gambRAW_manager::ionosphereKey(const string &sat, int node)
    {
        ostringstream os;
        os << sat << ":N" << node;
        return os.str();
    }
}
