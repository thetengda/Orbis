/**
 * @file         gcfg_ign.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        control set from XML for main
 * @version      1.0
 * @date         2025-01-01
 *
 * @copyright Copyright (c) 2025, Wuhan University. All rights reserved.
 *
 */

#ifndef GCFG_IGN_H
#define GCFG_IGN_H


#include "gdata/gimudata.h"
#include "gdata/gifcb.h"
#include "gcoders/ifcb.h"
#include "gcoders/imufile.h"
#include "gset/gsetins.h"
#include "gset/gsetign.h"
#include "gset/gcfg_ppp.h"
#include "gmsf/gintegration.h"


using namespace std;
using namespace gnut;
using namespace great;

class t_gcfg_ign : public virtual t_gcfg_ppp,
                   public virtual t_gsetins,
                   public virtual t_gsetign
{
    public:
        t_gcfg_ign();
        ~t_gcfg_ign();
        
        void check();
        void help();

    protected:
        
};





#endif