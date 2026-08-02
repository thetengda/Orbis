/**
 * @file         gintegration.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        loosely/tightly coupled integration
 * @version      1.0
 * @date         2025-01-01
 *
 * @copyright Copyright (c) 2025, Wuhan University. All rights reserved.
 *
 */

#ifndef GINTEGRATION_H
#define GINTEGRATION_H

#include "gins/gins.h"
#include "gutils/gmutex.h"
#include "gset/gsetbase.h"
#include "gdata/gimudata.h"
#include "gio/grtlog.h"
#include "gio/giof.h"
#include "gset/gsetign.h"
#include "gproc/gpvtflt.h"
#include "gdata/gposdata.h"
#include "gmsf/gpublish.h"
#include "gexport/ExportLibGREAT.h"
#include <set>
#include <string>

using namespace gnut;

namespace great
{

    /**
    * @class t_gintegration
    * @brief t_gintegration Class for loosely-coupled and tightly-coupled integration
    */
    class LibGREAT_LIBRARY_EXPORT t_gintegration : virtual public t_gsinskf,
                                              virtual public t_gpvtflt
    {
    public:
        /**
        * @brief Constructor of integration navigation class
        * set parameter value by rover, base, gset, data
        * @param[in]  site            rover's name.
        * @param[in]  site_base        base's name.
        * @param[in]  gset            pointer of user settings from xml file
        * @param[in]  allproc        pointer of data to process
        */
        explicit t_gintegration(string site, string site_base, t_gsetbase* gset, std::shared_ptr<spdlog::logger> spdlog, t_gallproc* allproc);

        /**
        * @brief process main function
        * support two-way filter
        * @param[in]  beg            begin time in xml "gen"
        * @param[in]  end            end time in xml "gen"
        * @param[in]  beg_end        filter direction matches xml "fltmode"
        * @return
            @retval <0        process fail
            @retval =0        normal end
        */
        virtual int processBatchFB(const t_gtime& beg, const t_gtime& end, bool beg_end);

    protected:

        /**
        * @brief init parameters, values and times in integration navigation filter
        * @return
            @retval -1        init fail
            @retval  1        normal end
        *
        */
        virtual int _init();

        /**
        * @brief init parameters, values in satellite navigation of integration navigation
        * subprocess of _init(), combine ins states with spp states (crd parameters are common)
        * @return
            @retval    false    init fail
            @retval     true    normal end
        *
        */
        virtual bool _gnss_init();

        /**
        * @brief process one epoch in GNSS/INS TCI
        * virual function, inherited from t_gsinskf, it is overridden in derived class t_gintegration
        * @param[in]  runEpoch        reference of GNSS process current time
        * @return
            @retval    MEAS_TYPE    indicates which type obs used in current observation update
        */
        virtual int _processEpoch(const t_gtime& runEpoch) override;

        /**
        * @brief TC GNSS/INS data fusion observation matrix A 
        * @param[in]  A  observation design matrix in GNSS filter
        * @return
            @retval -1        abnormal end
            @retval  1        normal end
        */
        virtual int _merge_pose(Matrix& A);

        /**
        * @brief TC GNSS/INS init data fusion
        * covariance of pos and other GNSS parameters are negated due to pos definition
        * @return
            @retval  1        normal end
        */
        virtual int _merge_init();

        /**
        * @brief GNSS observation update in GNSS/INS LCI, STCI and TCI obsevation update
        * automatic decide which method to take according to _ign_type
        * @return
            @retval -1        abnormal end
            @retval  1        normal end
        *
        */
        virtual int _GNSS_Update();


        /**
        * @brief get GNSS filter results
        * maybe only used in align, lci, stci
        * @return
            @retval    MEAS_TYPE    indicates which type obs used in observation update
                                (related to whether doppler is used or not)
        */
        virtual MEAS_TYPE _getPOS(t_gposdata::data_pos& pos);

        /**
        * @brief acquire all type obs at current epoch and insert into _Meas_Type
        * @return
            @retval        int        size of _Meas_Type, indicates how many types obs used in observation update
        */
        virtual int _getMeas();

        /**
        * @brief judge whether the current time is an effective data fusion time
        * related to xml "delay_t", care about ins type (for example fsas 0.026 and starneto 0.001)
        * @param[in]  gt            GNSS data process time
        * @param[in]  inst            INS data process time
        * @return
            @retval        bool    whether this epoch is valid
        *
        */
        virtual bool _time_valid(t_gtime gt, t_gtime inst);

        /**
        * @brief feedback GNSS parameters in filter
        * @note only used in TCI
        * @param[in]  dx            solution after filter obs update
        * @return
            @retval -1        abnormal end
            @retval  1        normal end
        *
        */
        virtual int _gnss_feedback(const ColumnVector& dx);

        virtual void _feedback(MEAS_TYPE type, bool successed);
        
        /**
        * @brief feedback GNSS parameters in filter
        * @param[in]  info            output info
        */
        virtual void _write(string info = "");

        /**
        * @brief write trajectory to KML file (xml need to set)
        */
        virtual int _prt_ins_kml();

    protected:
        t_gtime _ins_beg, _gnss_beg;
        t_gtime _ins_end, _gnss_end;
        t_gtime _gnss_crt;
        bool _initial_merge;

        t_gpublish _publisher;
        IMUState _imu_state;


    };

}




#endif
