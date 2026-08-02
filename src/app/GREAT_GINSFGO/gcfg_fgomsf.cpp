#include "gcfg_fgomsf.h"

t_gcfg_fgomsf::t_gcfg_fgomsf():
	t_gcfg_ppp(),
	t_gsetins(),	
	t_gsetign(),
	t_gsetfgo()
{
	_IFMT_supported.insert(IMU_INP);

	_OFMT_supported.insert(INS_OUT);
}

t_gcfg_fgomsf::~t_gcfg_fgomsf()
{
}

void t_gcfg_fgomsf::check()
{
	t_gcfg_ppp::check();
    t_gsetins::check();	
	t_gsetign::check();
	t_gsetfgo::check();

}

void t_gcfg_fgomsf::help()
{
	t_gcfg_ppp::help();
	t_gsetins::help();	
	t_gsetign::help();
	t_gsetfgo::help();
}
