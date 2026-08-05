#ifndef FGO_RAW_FACTOR_COMMON_H
#define FGO_RAW_FACTOR_COMMON_H

#include "gfactor/raw_equation.h"
#include "gfgo/gprecisebiasFGO.h"

#include <cmath>

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
        double amb = 0.0;
        double coefficient[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        double isb_coefficient = 0.0;
        double amb_coefficient = 0.0;
        double residual = 0.0;
        double sqrt_info = 0.0;
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
                          const t_gallpar &base_params,
                          t_gprecisebiasFGO *bias_model,
                          const double crd[3],
                          double clk,
                          double trp,
                          double sion,
                          double isb,
                          double amb,
                          bool use_isb,
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
        out.amb = amb;

        t_gallpar params = base_params;
        const string site = message.site.empty() ? message.satdata.site() : message.site;
        const string sat = message.sat_id.empty() ? message.satdata.sat() : message.sat_id;

        setOrAdd(params, site, par_type::CRD_X, "", out.crd[0]);
        setOrAdd(params, site, par_type::CRD_Y, "", out.crd[1]);
        setOrAdd(params, site, par_type::CRD_Z, "", out.crd[2]);
        setOrAdd(params, site, par_type::CLK, "", out.clk);
        setOrAdd(params, site, par_type::TRP, "", out.trp);
        setOrAdd(params, site, par_type::SION, sat, sion);
        if (use_isb)
        {
            const par_type type = isbType(message.satdata.gsys());
            if (type != par_type::NO_DEF)
                setOrAdd(params, site, type, "", isb);
        }

        t_gsatdata obsdata = message.satdata;
        t_gobs gobs(message.obs);
        t_gbaseEquation equation;
        t_gtime epoch = message.time;
        if (!bias_model->cmb_equ(false, true, epoch, params, obsdata, gobs, equation) ||
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
}
}

#endif
