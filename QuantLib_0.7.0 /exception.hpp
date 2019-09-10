#ifndef EXCEPTION_HPP
#define EXCEPTION_HPP

#include <exception>
#include <sstream>
#include <string>

using namespace std;

class Exception;

/** Streaming operator for class Exception
*   T can be ostream, string
*/


template <typename T>
Exception & operator<<( Exception & ex, T const & msg);


/** Exception class */

class Exception: public exception {
public:
    /** initializing ctor */
    explicit Exception(const string & msg): what_(msg) {}
    /** Returns the error msg */
    virtual const char* what() const noexcept override { return what_.c_str(); }
    /** dtor */
    virtual ~Exception() {}

protected:
    mutable string what_;

    template <typename T>
    friend Exception & operator<<( Exception & ex, T const & msg);
};

template <typename T>
Exception & operator<<( Exception & ex, T const & msg) {
    ostringstream s;
    s << msg;
    ex.what_ += s.str();
    return ex;
}


#define ASSERT(condition, errmsg) \
    if (!(condition)) {           \
    if ((string(errmsg)).empty()) \
        throw Exception("error: " #condition); \
    else \
        throw Exception(errmsg);  \
    }

#endif // EXCEPTION_HPP
