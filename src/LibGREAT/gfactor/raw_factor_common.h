#ifndef FGO_RAW_FACTOR_COMMON_H
#define FGO_RAW_FACTOR_COMMON_H

#include "gfactor/raw_equation.h"
#include "gfgo/gprecisebiasFGO.h"

#include <array>
#include <cmath>
#include <mutex>
#include <vector>

namespace gfgo
{
namespace raw_factor_detail
{
    struct RawLinearization
    {
        double crd[3] = {0.0, 0.0, 0.0};
        double clk = 0.0;
        double trp = 0.0;
        double sion = 0.0;
        double isb = 0.0;
        double ifb = 0.0;
        double amb = 0.0;
        double coefficient[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        double isb_coefficient = 0.0;
        double ifb_coefficient = 0.0;
        double amb_coefficient = 0.0;
        double residual = 0.0;
        double sqrt_info = 0.0;
    };

    struct RawEvaluationCache
    {
        std::mutex mutex;
        bool valid = false;
        bool frozen = false;
        std::array<double, 9> state{{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
        RawLinearization linearization;

        void reset()
        {
            std::lock_guard<std::mutex> lock(mutex);
            valid = false;
            frozen = false;
        }

        bool freeze()
        {
            std::lock_guard<std::mutex> lock(mutex);
            frozen = valid;
            return frozen;
        }
    };

    class RawPreparableFactor
    {
    public:
        virtual ~RawPreparableFactor() = default;
        virtual bool prepare(const std::vector<double *> &blocks) = 0;
    };

    inline par_type isbType(GSYS system)
    {
        switch (system)
        {
        case GSYS::GAL:
            return par_type::GAL_ISB;
        case GSYS::BDS:
            return par_type::BDS_ISB;
        case GSYS::GLO:
            return par_type::GLO_ISB;
        case GSYS::QZS:
            return par_type::QZS_ISB;
        default:
            return par_type::NO_DEF;
        }
    }

    inline par_type ifbType(GSYS system, FREQ_SEQ frequency)
    {
        switch (system)
        {
        case GSYS::GPS:
            return frequency == FREQ_3 ? par_type::IFB_GPS : par_type::NO_DEF;
        case GSYS::GAL:
            if (frequency == FREQ_3) return par_type::IFB_GAL;
            if (frequency == FREQ_4) return par_type::IFB_GAL_2;
            if (frequency == FREQ_5) return par_type::IFB_GAL_3;
            return par_type::NO_DEF;
        case GSYS::BDS:
            if (frequency == FREQ_3) return par_type::IFB_BDS;
            if (frequency == FREQ_4) return par_type::IFB_BDS_2;
            if (frequency == FREQ_5) return par_type::IFB_BDS_3;
            return par_type::NO_DEF;
        case GSYS::QZS:
            return frequency == FREQ_3 ? par_type::IFB_QZS : par_type::NO_DEF;
        default:
            return par_type::NO_DEF;
        }
    }

    inline void setOrAdd(t_gallpar &params, const string &site, par_type type,
                         const string &prn, double value)
    {
        int index = params.getParam(site, type, prn);
        if (index < 0)
        {
            t_gpar par(site, type, params.parNumber() + 1, prn);
            par.value(value);
            par.apriori(value);
            params.addParam(par);
            params.reIndex();
        }
        else
        {
            params[index].value(value);
        }
    }

    inline bool linearize(const RAWEquMsg &message,
                          t_gallpar &params,
                          t_gprecisebiasFGO *bias_model,
                          const double crd[3],
                          double clk,
                          double trp,
                          double sion,
                          double isb,
                          double ifb,
                          double amb,
                          bool use_isb,
                          bool use_ifb,
                          bool use_amb,
                          RawLinearization &out)
    {
        if (!bias_model || message.obs == GOBS::X)
            return false;

        out = RawLinearization();
        out.crd[0] = crd[0];
        out.crd[1] = crd[1];
        out.crd[2] = crd[2];
        out.clk = clk;
        out.trp = trp;
        out.sion = sion;
        out.isb = isb;
        out.ifb = ifb;
        out.amb = amb;

        const string site = message.site.empty() ? message.satdata.site() : message.site;
        const string sat = message.sat_id.empty() ? message.satdata.sat() : message.sat_id;

		// t_gprecisebiasFGO runs with estimator=FGO and therefore adds the
		// receiver eccentricity in _apply_rec_RTK().  Convert the graph's ARP
		// coordinate to its marker-coordinate input; the offset is constant, so
		// the coordinate Jacobian remains unchanged.
        setOrAdd(params, site, par_type::CRD_X, "",
				 out.crd[0] - message.receiver_eccentricity[0]);
        setOrAdd(params, site, par_type::CRD_Y, "",
				 out.crd[1] - message.receiver_eccentricity[1]);
        setOrAdd(params, site, par_type::CRD_Z, "",
				 out.crd[2] - message.receiver_eccentricity[2]);
        setOrAdd(params, site, par_type::CLK, "", out.clk);
        setOrAdd(params, site, par_type::TRP, "", out.trp);
        setOrAdd(params, site, par_type::SION, sat, sion);
        if (use_isb)
        {
            const par_type type = isbType(message.satdata.gsys());
            if (type != par_type::NO_DEF)
                setOrAdd(params, site, type, "", isb);
        }
        if (use_ifb)
        {
            const par_type type = ifbType(message.satdata.gsys(), message.freq);
            if (type == par_type::NO_DEF)
                return false;
            setOrAdd(params, site, type, "", ifb);
        }

        t_gsatdata obsdata = message.satdata;
        t_gobs gobs(message.obs);
        t_gbaseEquation equation;
        t_gtime epoch = message.time;
        if (!bias_model->cmb_equ(false, true, epoch, params, obsdata, gobs, equation, false) ||
            equation.B.empty() || equation.P.empty() || equation.l.empty() ||
            equation.P.front() <= 0.0)
        {
            return false;
        }

        for (const auto &item : equation.B.front())
        {
            const int index = item.first - 1;
            if (index < 0 || index >= static_cast<int>(params.parNumber()))
                continue;

            const par_type type = params.getPar(index).parType;
            switch (type)
            {
            case par_type::CRD_X:
                out.coefficient[0] += item.second;
                break;
            case par_type::CRD_Y:
                out.coefficient[1] += item.second;
                break;
            case par_type::CRD_Z:
                out.coefficient[2] += item.second;
                break;
            case par_type::CLK:
                out.coefficient[3] += item.second;
                break;
            case par_type::TRP:
                out.coefficient[4] += item.second;
                break;
            case par_type::SION:
                out.coefficient[5] += item.second;
                break;
            case par_type::GLO_ISB:
            case par_type::GAL_ISB:
            case par_type::BDS_ISB:
            case par_type::QZS_ISB:
                out.isb_coefficient += item.second;
                break;
            case par_type::IFB_GPS:
            case par_type::IFB_GAL:
            case par_type::IFB_GAL_2:
            case par_type::IFB_GAL_3:
            case par_type::IFB_BDS:
            case par_type::IFB_BDS_2:
            case par_type::IFB_BDS_3:
            case par_type::IFB_QZS:
                out.ifb_coefficient += item.second;
                break;
            default:
                break;
            }
        }

        out.residual = equation.l.front() + message.additive_correction;
        if (use_amb)
        {
            out.residual -= amb;
            out.amb_coefficient = 1.0;
        }
        out.sqrt_info = sqrt(equation.P.front());
        return std::isfinite(out.sqrt_info) && out.sqrt_info > 0.0;
    }

    inline bool linearizeCached(const RAWEquMsg &message,
                                t_gallpar &params,
                                t_gprecisebiasFGO *bias_model,
                                const double crd[3],
                                double clk,
                                double trp,
                                double sion,
                                double isb,
                                double ifb,
                                double amb,
                                bool use_isb,
                                bool use_ifb,
                                bool use_amb,
                                RawEvaluationCache &cache,
                                RawLinearization &out)
    {
        const std::array<double, 9> state{{crd[0], crd[1], crd[2], clk, trp,
                                           sion, use_isb ? isb : 0.0,
                                           use_ifb ? ifb : 0.0,
                                           use_amb ? amb : 0.0}};
        // prepare() freezes this cache before Ceres starts worker threads.
        // The solve then evaluates an immutable first-order observation model,
        // so no shared precise-model, navigation or windup state is touched.
        if (cache.frozen)
        {
            out = cache.linearization;
            double delta = 0.0;
            for (int i = 0; i < 3; ++i)
                delta += out.coefficient[i] * (state[i] - cache.state[i]);
            delta += out.coefficient[3] * (state[3] - cache.state[3]);
            delta += out.coefficient[4] * (state[4] - cache.state[4]);
            delta += out.coefficient[5] * (state[5] - cache.state[5]);
            delta += out.isb_coefficient * (state[6] - cache.state[6]);
            delta += out.ifb_coefficient * (state[7] - cache.state[7]);
            delta += out.amb_coefficient * (state[8] - cache.state[8]);
            out.residual -= delta;
            return std::isfinite(out.residual);
        }

        std::lock_guard<std::mutex> lock(cache.mutex);
        if (cache.valid && cache.state == state)
        {
            out = cache.linearization;
            return true;
        }

        RawLinearization linearization;
        if (!linearize(message, params, bias_model, crd, clk, trp, sion,
                       isb, ifb, amb, use_isb, use_ifb, use_amb, linearization))
        {
            cache.valid = false;
            return false;
        }

        cache.state = state;
        cache.linearization = linearization;
        cache.valid = true;
        out = linearization;
        return true;
    }

    template <typename Factor>
    inline bool prepareFactor(Factor &factor, RawEvaluationCache &cache,
                              const std::vector<double *> &blocks)
    {
        if (blocks.size() != factor.parameter_block_sizes().size())
            return false;
        cache.reset();
        std::vector<const double *> parameters;
        parameters.reserve(blocks.size());
        for (double *block : blocks)
        {
            if (!block)
                return false;
            parameters.push_back(block);
        }
        double residual = 0.0;
        if (!factor.Evaluate(parameters.data(), &residual, nullptr) ||
            !std::isfinite(residual))
            return false;
        return cache.freeze();
    }
}
}

#endif
