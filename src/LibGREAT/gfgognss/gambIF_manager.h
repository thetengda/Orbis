#ifndef AMBMANAGERIF_H
#define AMBMANAGERIF_H
/**
* @file gamb_manager_LC.h
* @details
* @note  Class for ambiguity management (IF)
* @verbatim
        History
        -0.1    hyChang        2022-05-17 creat the file.

  @endverbatim

* @author hyChang
* @version 0.1
* @date 2022-05-17
* @license
*/
#include "gfgognss/gamb_manager.h"
#include "gproc/gpvtflt.h"
#include "gexport/ExportLibGREAT.h"
namespace gfgomsf
{
	class t_gambIF_per_ID : public t_gamb_per_ID
	{
	public:
		t_gambIF_per_ID(pair<FREQ_SEQ, GOBSBAND> freq_band1, pair<FREQ_SEQ, GOBSBAND> freq_band2, GSYS sys, string prn, int sat_id, int amb_id, int start_rover, double initial_amb);
	protected:
		pair<FREQ_SEQ, GOBSBAND> _freq_band1;
		pair<FREQ_SEQ, GOBSBAND> _freq_band2;
	};

	class LibGREAT_LIBRARY_EXPORT t_gambIF_manager : public t_gamb_manager
	{
	public:
		t_gambIF_manager();
		t_gambIF_manager(map<GSYS, map<FREQ_SEQ, GOBSBAND>> band_index);
		~t_gambIF_manager();
		void clearState();
		void addNewSat(const t_gtime & cur_time, const int &rover_index, const int &sat_index, int &amb_index, const t_gsatdata &sat_data, t_gallpar params);
		bool addAmb(const t_gtime & cur_time, vector<t_gpar> amb_para, const  GSYS & gnss_system, int rover_index, const int &sat_index, int &amb_index);
		void addRover(string sat_name, const int &rover_index);
		void addRover(double time, string sat_name, const int &rover_index);//add for MultiWindow
		void updateAmb(int amb_id, double value);
		double getAmb(int amb_id);
		double getInitialAmb(const int & amb_id);
		void generateAmbSearchIndex();
		int getAmbSearchIndex(const int& sat_id);
		void removeSat(const int &sat_global_id, const int &rover_index);
		void addGpara(t_gallpar & params, const int & idx, double value);
		vector<int> getMarginAmb();
		int getAmbStartRoverID(const int & amb_search_id);
		int getAmbEndRoverID(const int & amb_search_id);
		vector<int> getCurWinAmb();
		vector<int> getCurWinAmb(int window_type);
		void slidingWindow();
		void slideNew(const int & rover_number);
		void stopAmbPropagating(const int& rover_number);
		bool isEstimated(const int & amb_id);

	protected:
		map<int, shared_ptr<t_gambIF_per_ID>> _ambiguityIF;//amb_id -> t_gambIF_per_ID
		map<int, int> search_index_IF;//sat_id -> amb_id
	};
}
#endif