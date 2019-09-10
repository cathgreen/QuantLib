#ifndef SPTR_HPP
#define SPTR_HPP

#include "defines.hpp"
#include <memory>
#include <string>

using namespace std;

using SPtrInt = shared_ptr<int>;
using SPtrLong = shared_ptr<long>;
using SPtrDouble = shared_ptr<double>;
using SPtrString = shared_ptr<string>;

#endif // SPTR_HPP
