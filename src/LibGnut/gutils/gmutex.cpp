
/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 * 
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  This file is part of the G-Nut C++ library.
 
-*/

#include <iostream>
#include <iomanip>

#include "gutils/gmutex.h"

namespace gnut
{

    t_gmutex::t_gmutex()
    {
    }

    t_gmutex::t_gmutex(const t_gmutex &Other)
    {
        // A copied data object owns a new, initially unlocked mutex.
    }

    t_gmutex::~t_gmutex()
    {
    }

    t_gmutex t_gmutex::operator=(const t_gmutex &Other)
    {
        return t_gmutex();
    }

    void t_gmutex::lock()
    {
        _mutex.lock();
    }

    void t_gmutex::unlock()
    {
        _mutex.unlock();
    }

} // namespace
