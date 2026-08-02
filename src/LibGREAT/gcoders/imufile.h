/**
 * @file         imufile.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        read and decode imu file 
 * @version      1.0
 * @date         2025-01-01
 * @verbatim
 *   imufile is the INS observation file,including Gyro and Acce observation.
 *   The example of th imu file( No File Head ):
 *   ===================================================================================
 *       Sow     Gyro_X        Gyro_Y        Gyro_Z        Acce_X        Acce_Y        Acce_Z
 *    11111.01   0.2500        0.2500        0.2500        0.2500        0.2500        10.000
 *    11111.02   0.2500        0.2500        0.2500        0.2500        0.2500        10.000
 *    ...
 *    11200.90   0.2500        0.2500        0.2500        0.2500        0.2500        10.000
 *   ===================================================================================
 * @endverbatim
 * @detials
 *  The order of data record is different in different file, so provide some orders for different file,
 *  including a->g,g->a; rfu,flu,frd and so on.You can add extra order if need.
 *
 * @copyright Copyright (c) 2025, Wuhan University. All rights reserved.
 *
 */
#ifndef IMUFILE_H
#define IMUFILE_H

#include "gins/gutility.h"
#include "gcoders/gcoder.h"
#include "gutils/gstring.h"
#include <Eigen/Eigen>
#include "gexport/ExportLibGREAT.h"
#include "gset/gsetins.h"

using namespace std;
using namespace gnut;

namespace great
{

    /**
    * @class t_imufile
    * @brief t_imufile Class for decoding and encoding the IMU file
    *
    * imufile is used for sins mechanical arrangement.
    * The document contains information about gyro and acce observation.
    * The gcoder t_imufile corresponding to the gdata t_gimudata.
    */
    class LibGREAT_LIBRARY_EXPORT t_imufile : public t_gcoder
    {
    public:

        /**
        * @brief constructor.
        *
        * @param[in]  s        setbase control
        * @param[in]  version  version of the gcoder
        * @param[in]  sz       size of the buffer
        */
        explicit t_imufile(t_gsetbase* s, string version = "", int sz = DEFAULT_BUFFER_SIZE);

        /** 
        * @brief destructor. 
        */
        ~t_imufile() {}


        /**
        * @brief decode the header of the IMU data file.
        * imu file doesn't include head block
        *
        * @param[in]  buff        buffer of the data
        * @param[in]  sz          buffer size of the data
        * @param[in]  errmsg      error message of the data decoding
        * @return
            @retval >=0 consume size of header decoding
            @retval <0  finish reading
        */
        virtual int decode_head(char* buff, int sz, vector<string>& errmsg);

        /**
        * @brief decode the data body of the IMU data file.
        * the main entry of decode IMU file and call the corresponding sub function according to 'AxisOrder' in xml file 
        *
        * decode data body of IMU file, all the data read will store in the t_gimudata
        * imu data final result -> t(s) gyro(rad) acce(m/s) 
        *
        * @param[in]  buff        buffer of the data
        * @param[in]  sz          buffer size of the data
        * @param[in]  cnt         size of lines successfully decoded
        * @param[in]  errmsg      error message of the data decoding
        * @return
            @retval >=0 consume size of data decoding
            @retval <0  finish reading
        */
        virtual int decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg);

        /**
        * @brief decode the data body of the IMU data file with starneto ascii format.
        * order="starneto"
        * decode data body of IMU file, all the data read will store in the t_gimudata
        *
        * @param[in]  line        line in imu file
        * @param[in]  t              time of data epoch
        * @param[out] v           imu data(gyro & acce)

        * @return
            @retval 1   success
            @retval -1  fail
        */
        virtual int decode_starneto(const string& line, double& t, vector<double>& v);
        
        /**
        * @brief decode the data body of the IMU data file with starneto binary format.
        * decode data body of IMU file, all the data read will store in the t_gimudata
        *
        * @param[in]  line        line in imu file
        * @param[in]  t           time of data epoch
        * @param[out] v           imu data(gyro & acce)

        * @return
            @retval 1   success
            @retval -1  fail
        */
        virtual int decode_starneto(const char* block, int sz, vector<pair<double, vector<double>>>& v);

        /**
        * @brief judge whether IMU data is available
        * @param[in]  now         current epoch

        * @return
            @retval true   available
            @retval false  unavailable
        */
        virtual bool available(const t_gtime& now);

        /**
        * @brief encode the header of the IMU data file.
        * imu file doesn't include head block
        *
        * @param[in]  buff        buffer of the data
        * @param[in]  sz          buffer size of the data
        * @param[in]  errmsg      error message of the data decoding
        * @return
            @retval >=0 consume size of header decoding
            @retval <0  finish reading
        */
        virtual  int encode_head(char* buff, int sz, vector<string>& errmsg);

        /**
        * @brief encode the data body of the IMU data file.
        * the main entry of encode IMU file.
        *
        * encode data body of IMU file.
        *
        * @param[in]  buff        buffer of the data
        * @param[in]  sz          buffer size of the data
        * @param[in]  cnt         size of lines successfully encoded
        * @param[in]  errmsg      error message of the data encoding
        * @return
            @retval >=0 consume size of data encoding
            @retval <0  finish reading
        */
        virtual  int encode_data(char* buff, int sz, int& cnt, vector<string>& errmsg);


    private:
        double _tt;
        bool _complete;        ///< whether the reading completed 
        string _order;         ///< the imu record order
        double _ts;            ///< data interval
        double _freq;          ///< data interval
        UNIT_TYPE _GyroUnit;   ///< Gyro data Unit
        UNIT_TYPE _AcceUnit;   ///< Acce data Unit
        UNIT_TYPE _MagUnit;    ///< Mag data Unit
    };
}






#endif