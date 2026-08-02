#include "gcfg_ign.h"

t_gcfg_ign::t_gcfg_ign() :
    t_gcfg_ppp(),
    t_gsetins(),
    t_gsetign()
{
    _IFMT_supported.insert(IFMT::IMU_INP);

    _OFMT_supported.insert(INS_OUT);
}

t_gcfg_ign::~t_gcfg_ign()
{
}

void t_gcfg_ign::check()
{
    t_gcfg_ppp::check();
    t_gsetins::check();
    t_gsetign::check();
}

void t_gcfg_ign::help()
{
    t_gcfg_ppp::help();
    t_gsetins::help();
    t_gsetign::help();
}
