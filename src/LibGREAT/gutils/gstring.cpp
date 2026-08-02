/**
 * @file         gstring.cpp
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        for string
 * @version      1.0
 * @date         2024-08-29
 * 
 * @copyright Copyright (c) 2024, Wuhan University. All rights reserved.
 * 
 */
#include "gutils/gstring.h"

#include <iostream>
#include <sstream>
#include <stdarg.h>
#include <math.h>
#include <algorithm>

namespace great
{
    std::string format(const char *fmt, ...)
    {
        va_list ap;
        va_start(ap, fmt);
        int len = vsnprintf(nullptr, 0, fmt, ap);
        va_end(ap);
        std::string buf(len + 1, '\0');
        va_start(ap, fmt);
        vsnprintf(&buf[0], buf.size(), fmt, ap);
        va_end(ap);
        buf.pop_back();
        return buf;
    }

    void split(const std::string& s, std::string delim, std::vector<std::string>& ret)
    {
        size_t last = 0;
        size_t index = s.find_first_of(delim, last);

        while (index != std::string::npos)
        {
            ret.push_back(s.substr(last, index - last));
            last = index + 1;
            index = s.find_first_of(delim, last);
        }
        if (index - last > 0)
        {
            ret.push_back(s.substr(last, index - last));
        }
    }
}