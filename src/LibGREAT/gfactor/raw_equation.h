#ifndef FGO_RAW_EQUATION_H
#define FGO_RAW_EQUATION_H

#include "gutils/gtime.h"
#include "gutils/gnss.h"
#include "gdata/gsatdata.h"

namespace gfgo
{
    using namespace gnut;

    /**
     * One uncombined PPP observation equation.
     *
     * The selected observation is stored explicitly.  RAW_ALL must not
     * select a different observable again while the factor is being
     * evaluated, otherwise receiver-code changes and OSB application can
     * silently make the equation non-reproducible.
     */
    struct RAWEquMsg
    {
        t_gtime time;
        t_gsatdata satdata;
        GOBSTYPE obs_type = TYPE;
        GOBS obs = GOBS::X;
        FREQ_SEQ freq = FREQ_X;
        GOBSBAND band = BAND;
        string site;
        string sat_id;
        string amb_id;
        string ion_id;
        int sat_global_id = -1;
        int amb_index = -1;
    };
}

#endif
