
/**
*
* @verbatim
    History
    2013-08-14  JD: created

  @endverbatim
* Copyright (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
*
* @file        gmutex.h
* @brief       Purpose: mutual exclusion
* @author      JD
* @version     1.0.0
* @date        2013-08-14
*
*/

#ifndef MUTEX_H
#define MUTEX_H

#include "gexport/ExportLibGnut.h"

#include <thread>
#include <mutex>

namespace gnut
{
    /** @brief class for t_gmutex. */
    class LibGnut_LIBRARY_EXPORT t_gmutex
    {
    public:
        /** @brief default constructor. */
        t_gmutex();

        /** @brief copy constructor. */
        t_gmutex(const t_gmutex &Other);

        /** @brief default destructor. */
        ~t_gmutex();

        /** @brief override operator =. */
        t_gmutex operator=(const t_gmutex &Other);

        /** @brief lock. */
        void lock();

        /** @brief unlock. */
        void unlock();

        // Retained for source compatibility only. The old value was racy and
        // must not be used to decide whether another thread may enter.
        bool isLock = false;

    protected:
        // Several legacy call paths re-enter the same data-object lock. A
        // recursive mutex preserves that behavior while still blocking other
        // threads; the former shared bool bypass was not mutual exclusion.
        std::recursive_mutex _mutex;
    };
} // namespace

#endif // MUTEX_H
