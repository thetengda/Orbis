#ifndef GAMBRAW_MANAGER_H
#define GAMBRAW_MANAGER_H

#include "gfgognss/gamb_manager.h"

namespace gfgomsf
{
    /**
     * RAW PPP ambiguity manager.
     *
     * RAW keeps one stable satellite ID for the lifetime of a tracked
     * satellite. Each carrier frequency owns an independent ambiguity arc;
     * a slip therefore retires only that frequency's arc and does not change
     * the satellite/SION address used by the factor graph.
     */
    class LibGREAT_LIBRARY_EXPORT t_gambRAW_manager : public t_gamb_manager
    {
    public:
        t_gambRAW_manager();
        explicit t_gambRAW_manager(const map<GSYS, map<FREQ_SEQ, GOBSBAND>> &band_index);
        ~t_gambRAW_manager();

        bool addNewSat(const t_gtime &cur_time, const int &rover_index,
                       const int &sat_index, int &amb_index,
                       const t_gsatdata &sat_data, t_gallpar params);
        bool addAmb(const t_gtime &cur_time, vector<t_gpar> amb_para,
                    const GSYS &gnss_system, int rover_index,
                    const int &sat_index, int &amb_index);
        bool resetFrequencyArc(const t_gtime &cur_time, int rover_index,
                               int sat_index, FREQ_SEQ freq, int &amb_index,
                               t_gallpar params);

        void addRover(string sat_name, const int &rover_index);
        void addRover(double time, string sat_name, const int &rover_index);
        void removeSat(const int &sat_global_id, const int &rover_index);
        void slidingWindow(double outgoing_time);
        void generateAmbSearchIndex();

        int getLatestSatId(const string &sat_name) const;
        bool hasActiveArc(int sat_index, FREQ_SEQ freq) const;
        bool getArcTime(int amb_id, t_gtime &beg, t_gtime &end) const;

        static string ambiguityKey(const string &sat, FREQ_SEQ freq, int arc);
        static string ionosphereKey(const string &sat, int node);

    private:
        int activeAmbiguity(int sat_index, FREQ_SEQ freq) const;
        void eraseAmbiguity(int amb_id);
        void addRoverToLatest(const string &sat_name, int rover_index,
                              double time, bool append_time);
    };
}

#endif
