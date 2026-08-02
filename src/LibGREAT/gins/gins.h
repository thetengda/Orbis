/**
 * @file         gins.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        INS core algorithm
 * @version      1.0
 * @date         2025-01-01
 * @details
                 SINS update algorithm,including AlignCoarse, Mechanization.
                 SINS integration class inskf,provide some interface to LC/TC.

 * @note         Acknowledgement to Prof. Gongmin Yan, Northwestern Polytechnical University
 * @cite         PSINS, Precise Strapdown Inertial Navigation System (https://psins.org.cn/)
 *
 * @copyright Copyright (c) 2025, Wuhan University. All rights reserved.
 *
 */

#ifndef GINS_H
#define GINS_H

#include "gimu.h"
#include "gset/gsetign.h"
#include "gdata/gimudata.h"
#include "gutils/gmutex.h"
#include "gall/gallpar.h"
#include "gio/gxml.h"
//#include "gins/gkalman.h"
using namespace gnut;

namespace great
{
    /**
    * @class t_gsins
    * @brief t_gsins Class is ins core class for ins updating.
    */
    class LibGREAT_LIBRARY_EXPORT t_gsins
    {
    public:

        /**
        * @brief Constructor
        * @note set parameter value
        */
        explicit t_gsins(const t_gquat& qnb0 = t_gquat(), const Eigen::Vector3d& vn0 = Eigen::Vector3d::Zero(), const Eigen::Vector3d& pos0 = Eigen::Vector3d::Zero(), double t0 = 0.0);


        /**
        * @brief Constructor
        * @note set parameter value by gset
        */
        explicit t_gsins(gnut::t_gsetbase* set);


        /**
        * @brief Constructor
        * @note set parameter value by gset
        */
        void set_posvel(const Eigen::Vector3d& pos0 = Eigen::Vector3d::Zero(), const Eigen::Vector3d& vn0 = Eigen::Vector3d::Zero());

        /**
        * @brief coarse alignment
        *
        * @attention data must be static
        * @param[in] wmm    angle increment
        * @param[in] vmm    velocity increment
        * @return
            initial attitude
        */
        Eigen::Vector3d align_coarse(const Eigen::Vector3d& wmm, const Eigen::Vector3d& vmm);

        /**
        * @brief sins mechanical arrangement
        *
        * @param[in] wm        angle increment
        * @param[in] vm        velocity increment
        * @param[in] scm    scheme of processing
        *
        */
        void Update(const vector<Eigen::Vector3d>& wm, const vector<Eigen::Vector3d>& vm, const t_scheme& scm);

        /**
        * @brief print ins file header into otringstream
        */
        void prt_header(ostringstream& os, bool imu_scale = false, bool use_odo = false);


        /**
        * @brief print ins parameter into otringstream
        *
        * @verbatim
            position is exported to (X,Y,Z),(B,L,H),....
            velocity is exported to (Vx,Vy,Vz),(Ve,Vn,VU),...
            attitude is exported to (pitch,roll,yaw[-180-180]),...
            eb is exported to (ebx,eby,ebz)[deg/h]
            db is exported to (dbx,dby,dbz)[mg]
          @endverbatim
        *
        * @param[in] os   out os
        *
        */
        void prt_sins(ostringstream& os);
        //void int_sec(bool _beg_end = true);
        void debug_ins_info();

    public:
        double nts, t;                                /// ultimate interval,now
        t_gimu imu;                                    /// imu data
        t_gearth eth;                                /// earth info
        Eigen::Vector3d att, vn, pos, an;            /// attitude,velocity,position,accelarate
        Eigen::Vector3d fb, fn, wib, web, wnb, vb;    /// special force,angular
        Eigen::Vector3d eb, db, Kg, Ka;                /// gyro and acce bias and scale factors
        Eigen::Vector3d _tauG, _tauA, 
            _tauGScale, _tauAScale;                    /// Markov time
        Eigen::Matrix3d Cnb, Cbn;                    /// DCM,Scale
        t_gquat qnb;                                /// quat

        Eigen::Vector3d ve, ae, fe, pos_ecef, xyz_out;
        Eigen::Matrix3d Ceb, Cbe;                    /// DCM
        t_gquat qeb;
        Eigen::Matrix3d Cvb, Cbv;                    ///< IMU installation Rotation 
        t_gquat qvb;



        double _tauOdo;
        double pure_ins_time;
    };



    /**
    * @class t_gsinskf
    * @brief t_gsinskf Class is SINS integration base class for coupled integration
    */
    class LibGREAT_LIBRARY_EXPORT t_gsinskf
    {
    public:

        /**
        * @brief Constructor
        * @note set parameter value by gset
        */


        t_gsinskf(gnut::t_gsetbase* gset, t_spdlog spdlog, string name = "");

        /** @brief default destructor. */
        virtual ~t_gsinskf();

        /**
        * @brief add imudata to integration navigation class
        * @param[in]  imu            pointer if t_gimudata class
        * @return                    void
        */
        virtual void Add_IMU(t_gimudata* imu);

        virtual void init_par(string site);
        virtual int get_last_idx();
        /**
        * @brief Initial Alignment
        */
        virtual bool align_coarse(const vector<Eigen::Vector3d>& wm, const vector<Eigen::Vector3d>& vm, const t_scheme& scm);

        /**
        * @brief Initial Alignment by position vector
        */
        virtual bool align_pva(const Eigen::Vector3d& pos, const t_scheme& scm);

        /**
        * @brief Initial Alignment by velocity vector
        */
        virtual bool align_vva(const Eigen::Vector3d& vel, const t_scheme& scm);

        /**
        * @brief print Pk into otringstream
        *
        * @param[in] t				time
        * @param[in] Pk				vector of towrite
        * @param[out] os			ostringstream
        *
        */
        /**
        * @brief Initial Alignment.
        * @details fast cascaded initial alignment method, including sa, pva, vva.
        */
        virtual bool cascaded_align(const Eigen::Vector3d& pos, const Eigen::Vector3d& vel);

        /**
        * @brief set Ft
        */
        virtual void set_Ft();

        /**
        * @brief set measurement matrix Hk
        */
        virtual void set_Hk();

        /** @brief body motion state. */
        MOTION_TYPE motion_state();

        /** @brief body measurement state. */
        MEAS_TYPE meas_state();

        /**
        * @brief write info into file
        */
        virtual void write();


		virtual int _prt_ins_kml() { return -1; }

        /** @brief setting. */
        void gset(gnut::t_gsetbase* set) { _setkf = set; };

        /** @brief output setting. */
        void set_out();

        double kftk;                                         /// now
        int nq, nr;                                          /// dimension of state and meas
        int measflag;                                        /// flag of measurement index
        Eigen::MatrixXd Ft, Pk, Hk, Rk, Phik;                /// Kalman Matrix
        Eigen::VectorXd Xk, Zk, Qt;                          /// Kalman Vector
        // Read-only access to a diagnostic-only process-noise accumulator.
        // The storage is kept outside this exported class so the diagnostic
        // does not change its ABI or the actual filter covariance Pk.
        Eigen::MatrixXd diagnostic_process_covariance() const;
        double diagnostic_process_dt() const;
        int diagnostic_process_steps() const;
        void reset_diagnostic_process_covariance();
        MEAS_TYPE Flag;                                /// Measurement Flag,see also MEAS_TYPE
        IGN_TYPE _ign_type;

        Eigen::Vector3d lever, odo_lever, uwb_lever;            /// arm lever (IMU->GNSS, IMU->DMI, IMU->UWB)
        Eigen::Vector3d    MeasVel, MeasPos, MeasAtt;  /// meas velosity and position  & attitude addwh
        Eigen::Vector3d    _Cov_MeasVn, _Cov_MeasPos, _Cov_MeasAtt;  /// meas velosity and position convariance
        Eigen::Vector3d _Cov_MeasNHC, _Cov_MeasZUPT;
        double tmeas, MeasYaw, MeasVf, MeasHgt;        /// meas time,meas yaw,meas velocity foward,meas height
        double _Cov_MeasZIHR, _Cov_MeasOdo, _Cov_MeasYaw;
        map<double, double> _map_yaw;
        t_gsins sins;               
        t_gallpar param_of_sins;
    
    protected:

        /**
        * @brief init parameters, values in inertial navigation of integration navigation
        * @note a subprocess of _init()
        * @return                    function running state (true or false)
        */
        virtual bool _ins_init();

        /**
        * @brief add inertial position constraint to GNSS filter
        * @note only used in stci or tci
        * @return                    function running state (true or false)
        */
        virtual bool _valid_ins_constraint();

        /**
        * @brief time update of kalman filter
        * @param[in]  kfts            interval (for example 0.005 in 200Hz imu)
        * @param[in]  fback            whether feedback (0 represent no, 1 represent yes)
        * @param[in]  inflation        generally set to 1
        * @return                    void
        */
        virtual void time_update(double kfts, double inflation = 1.0);

        /**
        * @brief means update of kalman filter
        * @note only used in LCI, STCI
        * @param[in]  fading        generally set to 1
        * @return                    function running state (1 represents normal end, -1 represents abnormal)
        */
        virtual int _meas_update();

        /**
        * @brief Hk, Zk and Rk are set according to the type of observation
        * @note only used in LCI, STCI
        * @return                    true or false
        */
        virtual bool _set_meas();

        /**
        * @brief automatic adjust the dimension of matrix Hk, Zk and Rk according to the type of observation
        * @note only used in LCI, STCI
        * @return                    true or false
        */
        virtual bool _resize();


        /**
        * @brief normalized posterior residual test
        * @note only used in LCI
        * @return                    function running state (0 represents normal, 1 represents residual overrun)
        */
        virtual int outlier_detect();

        /**
        * @brief kalman filter feedback(full closed loop correction)
        * @note definition of IMU bias is different
        * @param[in]  kfts            interval (for example 0.005 in 200Hz imu)
        * @return                    void
        */
        virtual void feedback();

        int _align_count;             // aligncoarse data number
        bool _aligned;                // aligned
        double _yaw0 = 0.0;
        bool _first_align;
        double _resample_intv;
        string _name = "";
        Eigen::Vector3d _first_pos, _pre_pos;
        Eigen::Vector3d wmm, vmm;     // total wm/vm increment (in alignment)
        vector<Eigen::Vector3d> _wm, _vm, _mm;//  wm/vm increment (every epoch)
        Eigen::MatrixXd _global_variance, _gv_sav;      // store the global variance infomation
        Eigen::VectorXd _v_post;                            /// post residual

        t_gimudata* _imudata;
        t_gtime _ins_crt;
        t_scheme _shm;                // processing scheme (important!)
        gnut::t_gsetbase* _setkf;    // set
        gnut::t_giof* _fins;        // ins t_giof
	
        t_spdlog _spdlogkf;
        //gnut::t_glog* _logkf;

        set<MEAS_TYPE> _Meas_Type;
        map<MEAS_TYPE, double> _map_maxnorm;

    };
}


#endif
