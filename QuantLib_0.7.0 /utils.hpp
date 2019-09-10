#ifndef UTILS_HPP
#define UTILS_HPP

#include "defines.hpp"
#include "exception.hpp"
#include <string>
#include <algorithm>
#include <ctype.h>
#include <cmath>

inline
string trim(const string& s) {
    auto wsfront = find_if_not(s.cbegin(), s.cend(), ::isspace);
    auto wsback = find_if_not(s.crbegin(), s.crend(), ::isspace).base();
    return (wsfront >= wsback ? string() : string(wsfront, wsback)); 
}
    
inline
double toContCmpd(double rate, size_t annfreq) {
    ASSERT(annfreq >= 1.0, "compounding frequency less than 1 not allowed");
    double tmp = pow(1.0 + rate / annfreq, annfreq);
    return log(tmp);
}

inline
double fromContCmpd(double rate, size_t annfreq) {
    ASSERT(annfreq >= 1.0, "compounding frequency less than 1 not allowed");
    double tmp = exp(rate / annfreq);
    return (tmp - 1.0) * annfreq;
}

#endif // UTILS_HPP
