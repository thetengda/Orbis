/**
 * @file         GREAT_GINSFGO.cpp
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        Main Function for tightly combination RTK/INS basecd on factor graph optimization
 * @version      1.0
 * @date         2025-11-04
 *
 * @copyright Copyright (c) 2025, Wuhan University. All rights reserved.
 *
 */

#include "gcfg_fgomsf.h"
#include <chrono>
#include <thread>

using namespace std;
using namespace gnut;
using namespace std::chrono;
using namespace gfgo;
using namespace gfgomsf;

void catch_signal(int) { std::cout << "Program interrupted by Ctrl-C [SIGINT,2]\n"; }

// MAIN
// ----------
int main(int argc, char** argv)
{
	// Only to cout the Reminder here
	signal(SIGINT, catch_signal);

	t_gcfg_fgomsf gset;
	gset.app("G-Nut/PVT", "0.9.0", "$Rev: 2448 $", "(gnss@pecny.cz)", __DATE__, __TIME__);
	// Get the arguments from the command line
	gset.arg(argc, argv, true, false);

	auto log_type = dynamic_cast<t_gsetout*>(&gset)->log_type();
	auto log_level = dynamic_cast<t_gsetout*>(&gset)->log_level();
	auto log_name = dynamic_cast<t_gsetout*>(&gset)->log_name();
	auto log_pattern = dynamic_cast<t_gsetout*>(&gset)->log_pattern();
	spdlog::set_level(log_level);
	spdlog::set_pattern(log_pattern);
	spdlog::flush_on(spdlog::level::err);
	t_grtlog great_log = t_grtlog(log_type, log_level, log_name);
	auto my_logger = great_log.spdlog();

	bool isBase = false;
	if (dynamic_cast<t_gsetgen*>(&gset)->list_base().size()) isBase = true;
	
	// Prepare site list from gset
	set<string> sites = dynamic_cast<t_gsetgen*>(&gset)->recs();
	// Prepare input files list form gset
	multimap<IFMT, string> inp = gset.inputs_all();
	// Get sample intval from gset. if not, init with the default value
	int sample = int(dynamic_cast<t_gsetgen*>(&gset)->sampling());
	if (!sample) sample = int(dynamic_cast<t_gsetgen*>(&gset)->sampling_default());

	t_gtime beg = dynamic_cast<t_gsetgen*>(&gset)->beg();
	t_gtime end = dynamic_cast<t_gsetgen*>(&gset)->end();

	// DECLARATIONS/INITIALIZATIONS
	// gobs for the obs data
	t_gallobs* gobs = new t_gallobs();  gobs->spdlog(my_logger); gobs->gset(&gset);
	// gallnav for all the navigation data, gorb = gorbit data
	t_gallprec* gorb = new t_gallprec(); gorb->spdlog(my_logger);

	t_gdata* gdata = 0;
	// gpcv for the atx data read from the atx file
	t_gallpcv* gpcv = 0; if (gset.input_size("atx") > 0) { gpcv = new t_gallpcv;  gpcv->spdlog(my_logger); }
	// gotl for the blq data read from the blq file, which will be used for Ocean tidal corrections
	t_gallotl* gotl = 0; if (gset.input_size("blq") > 0) { gotl = new t_gallotl;  gotl->spdlog(my_logger); }
	// gbia for the DCB data read from the biasinex and bianern files
	t_gallbias* gbia = 0; if (gset.input_size("biasinex") > 0 ||
		gset.input_size("bias") > 0) {
		gbia = new t_gallbias; gbia->spdlog(my_logger);
	}
	t_gupd* gupd = nullptr;
	if (gset.input_size("upd") > 0) {
		gupd = new t_gupd;  gupd->spdlog(my_logger);
	}
	// gobj for the gpcv and gotl, which means the model can be used by all satellites and stations
	//t_gallobj* gobj = new t_gallobj(gpcv, gotl); gobj->glog(&glog);

	t_gallobj* gobj = new t_gallobj(my_logger, gpcv, gotl); gobj->spdlog(my_logger);
	t_gnavde* gde = new t_gnavde;
	t_gpoleut1* gerp = new t_gpoleut1;
	t_gifcb* gifcb = nullptr; if (gset.input_size("ifcb") > 0) { gifcb = new t_gifcb;  gifcb->spdlog(my_logger); }
	// vgppp for the process of ppp with fliter

	t_gimudata* gimu = new t_gimudata(); gimu->spdlog(my_logger);
	
	vector<t_gfgo_gins*> vgins;

	// runepoch for the time costed each epoch (i guess)
	t_gtime runepoch(t_gtime::GPS);
	// lstepoch for the time of all epoches (i guess)
	t_gtime lstepoch(t_gtime::GPS);

	if (gset.input_size("sp3") == 0 && gset.input_size("rinexc") == 0)
	{
		gorb->use_clknav(true);
		gorb->use_posnav(true);
	}
	else if (gset.input_size("sp3") > 0 && gset.input_size("rinexc") == 0)
	{
		gorb->use_clksp3(true);
	}

	// SET OBJECTS
	set<string>::const_iterator itOBJ;
	set<string> obj = dynamic_cast<t_gsetrec*>(&gset)->objects();
	for (itOBJ = obj.begin(); itOBJ != obj.end(); ++itOBJ) {
		string name = *itOBJ;
		shared_ptr<t_grec> rec = dynamic_cast<t_gsetrec*>(&gset)->grec(name, my_logger);
		gobj->add(rec);
	}
	// Multi gcoder for multi-thread decoding data
	vector<t_gcoder*> gcoder;
	// Multi gior for multi-thread receiving data
	vector<t_gio*> gio;
	// multi-thread 
	//vector<thread> gthread;
	t_gio* tgio = 0;
	t_gcoder* tgcoder = 0;

	// DATA READING
	multimap<IFMT, string>::const_iterator itINP = inp.begin();
	for (size_t i = 0; i < inp.size() && itINP != inp.end(); ++i, ++itINP)
	{
		// Get the file format/path, which will be used in decoder
		IFMT   ifmt(itINP->first);
		string path(itINP->second);
		string id("ID" + int2str(i));

		// For different file format, we prepare different data container and decoder for them.
		// Decode bnc obs data
		//if (ifmt == BNCOBS_INP) { gdata = gobs; tgcoder = new t_bncobs(&gset, "", 40960); }
		if (ifmt == IMU_INP) { gdata = gimu; tgcoder = new t_imufile(&gset, "", 40960); }
		else if (ifmt == SP3_INP) { gdata = gorb; tgcoder = new t_sp3(&gset, "", 8172); }
		else if (ifmt == RINEXO_INP) { gdata = gobs; tgcoder = new t_rinexo(&gset, "", 4096); }
		else if (ifmt == RINEXC_INP) { gdata = gorb; tgcoder = new t_rinexc(&gset, "", 4096); }
		else if (ifmt == RINEXN_INP) { gdata = gorb; tgcoder = new t_rinexn(&gset, "", 4096); }
		else if (ifmt == BIASINEX_INP) { gdata = gbia; tgcoder = new t_biasinex(&gset, "", 20480); }
		else if (ifmt == ATX_INP) { gdata = gpcv; tgcoder = new t_atx(&gset, "", 4096); }
		else if (ifmt == BLQ_INP) { gdata = gotl; tgcoder = new t_blq(&gset, "", 4096); }
		else if (ifmt == UPD_INP) { gdata = gupd; tgcoder = new  t_upd(&gset, "", 4096); }
		else if (ifmt == DE_INP) { gdata = gde; tgcoder = new t_dvpteph405(&gset, "", 4096); }
		else if (ifmt == EOP_INP) { gdata = gerp; tgcoder = new t_poleut1(&gset, "", 4096); }
		else if (ifmt == IFCB_INP) { gdata = gifcb; tgcoder = new t_ifcb(&gset, "", 4096); }
		else {
			SPDLOG_LOGGER_ERROR(my_logger, "unrecognized format " + t_gsetinp::ifmt2str(ifmt));
			gdata = 0;
		}

		// Check the file path
		if (path.substr(0, 7) == "file://") {
			SPDLOG_LOGGER_ERROR(my_logger, "path is not file (skipped)!");
			tgio = new t_gfile(my_logger);
			tgio->spdlog(my_logger);
			tgio->path(path);
		}

		// READ DATA FROM FILE
		if (tgcoder) {
			// Put the file into gcoder
			tgcoder->clear();
			tgcoder->path(path);
			tgcoder->spdlog(my_logger);
			// Put the data container into gcoder
			tgcoder->add_data(id, gdata);
			if (ifmt != IMU_INP )tgcoder->add_data("OBJ", gobj);
			
			// Put the gcoder into the gio
			// Note, gcoder contain the gdata and gio contain the gcoder
			tgio->coder(tgcoder);

			
			
				runepoch = t_gtime::current_time(t_gtime::GPS);
				// Read the data from file here
				tgio->run_read();
				lstepoch = t_gtime::current_time(t_gtime::GPS);
				// Write the information of reading process to log file
				SPDLOG_LOGGER_INFO(my_logger, "main", "READ: " + path + " time: "
					+ dbl2str(lstepoch.diff(runepoch)) + " sec");
				// Delete 
				delete tgio;
				delete tgcoder;
			
		}
	}
	cout << "gimu:" << gimu->size() << endl;
	// set antennas for satllites (must be before PCV assigning)
	gobj->read_satinfo(beg);
	// assigning PCV pointers to objects
	gobj->sync_pcvs();
	//if (gbia == NULL) gbia = new t_gallbias();

	// add all data
	t_gallproc* data = new t_gallproc();
	data->Add_Data("gobs", gobs);
	data->Add_Data("gorb", gorb);
	data->Add_Data("gobj", gobj);
	if (gbia)data->Add_Data(t_gdata::type2str(gbia->id_type()), gbia);
	if (gotl)data->Add_Data(t_gdata::type2str(gotl->id_type()), gotl);
	if (gde)data->Add_Data(t_gdata::type2str(gde->id_type()), gde);
	if (gerp)data->Add_Data(t_gdata::type2str(gerp->id_type()), gerp);
	if (gupd && dynamic_cast<t_gsetamb*>(&gset)->fix_mode() == FIX_MODE::SEARCH && !isBase)
	{
		data->Add_Data("gupd", gupd);
	}

	int frequency = dynamic_cast<t_gsetproc*>(&gset)->frequency();
	set<string> system = dynamic_cast<t_gsetgen*>(&gset)->sys();
	if (frequency == 3 && system.find("GPS") != system.end() && !isBase)
	{
		data->Add_Data("gifcb", gifcb);
	}

	t_gtime::gmsleep(sample * 1000);
	auto tic_start = system_clock::now();

	// PVT PROCESSING - loop over sites from settings
	int i = 0, nsite = sites.size();
	if (isBase ) nsite = gset.list_rover().size();
	set<string>::iterator it = sites.begin();
	while (i < nsite) {
		string site_base = "";
		string site = *it;
		if (isBase) {
			site_base = (gset.list_base())[i];
			site = (gset.list_rover())[i];
			if ((gobs->beg_obs(site_base) == LAST_TIME || gobs->end_obs(site_base) == FIRST_TIME ||
				site_base.empty() || gobs->isSite(site_base) == false)) {
				SPDLOG_LOGGER_INFO(my_logger, "No two site/data for processing!");
				i++;
				continue;
			}
		}
		if (gobs->beg_obs(site) == LAST_TIME ||
			gobs->end_obs(site) == FIRST_TIME ||
			site.empty() || gobs->isSite(site) == false) {
			SPDLOG_LOGGER_INFO(my_logger, "No site/data for processing!");
			if (!isBase) it++;
			i++;
			continue;
		}

		size_t fgo_size = 1;
		bool _beg_end = true;
		int idx; 

		runepoch = t_gtime::current_time(t_gtime::GPS);


		for (size_t i = 0; i < fgo_size; i++)
		{
			vgins.push_back(0); idx = vgins.size() - 1;
			vgins[idx] = new t_gfgo_gins(site, site_base, &gset,my_logger, data);
			vgins[idx]->get_imu(gimu);
			SPDLOG_LOGGER_INFO(my_logger, "main", "Multi-thread fgo-gvil processing started");
			SPDLOG_LOGGER_INFO(my_logger, "main", beg.str_ymdhms("  beg: ") + end.str_ymdhms("  end: "));
			vgins[idx]->ProcessBatch(beg, end, _beg_end);
			_beg_end = false;
		}
		// Write the log file
		lstepoch = t_gtime::current_time(t_gtime::GPS);
		SPDLOG_LOGGER_INFO(my_logger, "main", site_base + site + "PVT/INS processing finished : duration  "
			+ dbl2str(lstepoch.diff(runepoch)) + " sec");
		if (!isBase) it++;
		i++;
	}

	for (size_t i = 0; i < gio.size(); ++i) { delete gio[i]; }; gio.clear();
	for (size_t i = 0; i < gcoder.size(); ++i) { delete gcoder[i]; }; gcoder.clear();

	for (unsigned int i = 0; i < vgins.size(); ++i) {
		if (vgins[i])
		{
			delete vgins[i];
		}
	}

	if (gobs) delete gobs;
	if (gpcv) delete gpcv;
	if (gotl) delete gotl;
	if (gobj) delete gobj;
	if (gorb) delete gorb;
	if (gbia) delete gbia;
	if (gde)  delete gde;
	if (gerp)  delete gerp;
	if (gupd)  delete gupd;
	if (gifcb)  delete gifcb;
	if (data)  delete data;
	if (gimu) delete gimu;

	auto tic_end = system_clock::now();
	auto duration = duration_cast<microseconds>(tic_end - tic_start);
	cout << "Spent" << double(duration.count()) * microseconds::period::num / microseconds::period::den << " seconds." << endl;

}