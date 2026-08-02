/**
*
* @verbatim
    History
    2011-03-25  JD: created

  @endverbatim
* Copyright (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
*
* @file        gdata.h
* @brief       about GNSS data
*.
* @author      JD
* @version     1.0.0
* @date        2011-03-25
*
*/

#ifndef GDATA_H
#define GDATA_H

#include <sstream>
#include <iostream>

#include "gutils/gmutex.h"
#include "gio/gnote.h"
#include "gutils/gcommon.h"
#include "gio/grtlog.h"
using namespace std;

namespace gnut
{

    /**
    *@brief       basic class for data storing derive from t_gmoint
    */
    class LibGnut_LIBRARY_EXPORT t_gdata
    {

    public:
        explicit t_gdata();
        /** @brief default constructor. */
        explicit t_gdata(t_spdlog spdlog);

        /** @brief copy constructor. */
        explicit t_gdata(const t_gdata &data);

        /** @brief default destructor. */
        virtual ~t_gdata();

        /** @brief override operator =. */
        t_gdata &operator=(const t_gdata &data);

        /** @brief data type */
        enum ID_TYPE
        {
            NONE,     ///< = 0,  none
            OBJ,     ///< = 1,  object
            TRN,     ///< = 2,  transmitter
            REC,     ///< = 3,  receiver

            OBS,     ///< = 10, obseravation base
            OBSGNSS, ///< = 11, gnss observations
            SATDATA, ///< = 12, gnss observations + satellite data


            EPH,     ///< = 20, navigation base
            EPHGPS,     ///< = 21, navigation
            EPHGLO,     ///< = 22, navigation
            EPHGAL,     ///< = 23, navigation
            EPHQZS,     ///< = 24, navigation
            EPHBDS,     ///< = 24, navigation
            EPHSBS,     ///< = 25, navigation
            EPHIRN,     ///< = 26, navigation
            EPHPREC, ///< = 27, sp3/clocks
           

            
            ALLNAV,     ///< = 28, all navigation all
            ALLPREC, ///< = 29, all sp3 + rinexc
            
            ALLOBS,             ///< = 31, all observations
            ALLOBJ,             ///< = 32, all objects
            ALLPCV,             ///< = 33, all PCV
            ALLOTL,             ///< = 34, all OTL
            
            ALLOPL,             ///< = xx, all ocean pole looad
            
            ALLPROD,         ///< = xx, all PROD
            ALLBIAS,         ///< = xx, all PROD
            
            ALLPOLEUT1,         ///< = xx, poleut1
            

            
            POS,     ///< = 35, XYZT position/time
            
            CLK,     ///< = 37, clocks
            
            IONEX,     ///< = xx, ionospheric delay from tec grid products (GIM)

            PCV,  ///< = 40, PCV model
            OTL,  ///< = 41, ocean loading model
            
            OPL,
            BIAS, ///< = 42, code & phase biases
           

           
            

            LCI_POS,     ///< = xx, lci position
            
            UPD,         ///< add for upd
            
            //UPD_EPOCH,
            IFCB,    ///< = xx, Inter-Frequency Clock Bias
            

            ALLDE,        ///< = xx, all planeteph
            

            
            LEAPSECOND,          ///< = xx, leap Second

            ALLPCVNEQ,          ///< = xx, all PCV NEG
            
            RSSIMAP,          ///< = xx, rssi map
           
            IMUDATA,     ///< = xx, imu data
            LAST              ///< = xx, last
        };

        /** @brief group type */
        enum ID_GROUP
        {
            GRP_NONE,    ///<   none
            GRP_OBSERV,  ///<   observations
            GRP_EPHEM,   ///<   ephemerides
            GRP_PRODUCT, ///<   positions/products
            GRP_MODEL,   ///<   models
            GRP_OBJECT,  ///<   objects
            GRP_LAST     ///<   last
        };

        /** @brief set glog pointer */
        void spdlog(t_spdlog spdlog);

        /** @brief get glog pointer */
        t_spdlog spdlog() const { return _spdlog; }

        /** @brief set gnote pointer */
        void gnote(t_gallnote *n) { _gnote = n; }

        /** @brief get gnote pointer */
        t_gallnote *gnote() const { return _gnote; }
        //void glog(t_glog* l) { _log = l; }

        /** @brief get glog pointer */
        //t_glog* glog() const { return _log; }

        /** @brief get data type */
        const ID_TYPE &id_type() const { return _type; }

        /** @brief get group type */
        const ID_GROUP &id_group() const { return _group; }

        /** @brief get group type in string format */
        string str_group() const;

        /** @brief convert data type to string */
        static string type2str(ID_TYPE type);

        /** @brief lock mutex */
        void lock() const { this->_gmutex.lock(); }

        /** @brief unlock mutex */
        void unlock() const { this->_gmutex.unlock(); }

    protected:
        /**
         * @brief data type
         * @param t 
         * @return int 
         */
        int id_type(const ID_TYPE &t); 

        /**
         * @brief group type
         * @param g 
         * @return int 
         */
        int id_group(const ID_GROUP &g);

        mutable t_gmutex _gmutex; ///< mutex
        t_gallnote *_gnote;       ///< gnote
        ID_TYPE _type;            ///< type_ID
        ID_GROUP _group;          ///< group_ID
        t_spdlog _spdlog;         ///< spdlog pointer
        //t_glog* _log;              ///< glog
    private:
    };
}
#endif
