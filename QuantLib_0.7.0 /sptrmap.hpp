#ifndef SPTRMAP_HPP
#define SPTRMAP_HPP

#include "sptr.hpp"
#include "utils.hpp"
#include <map>
#include <string>
#include <vector>
#include <cctype>
#include <algorithm>
#include <memory>


/**
   SPtrMap is used for Market class, to store states such as yieldcurves. 
   Usage example:
   market().yieldCurves().get(“USD”)->discount(T);
   string is the name of the yieldCurve, such as "USD"
   returns shared_ptr<T>: the ptr to YieldCurve object
   unsigned long: the current version of the map
*/

template<typename T>
class SPtrMap : public map<string, pair<shared_ptr<T>, unsigned long>> {
public:
    using ptr_type = shared_ptr<T>;
    using pair_type = pair<ptr_type, unsigned long>;
    using map_type = map<string, pair_type>;
    
    /** Returns a list of names of the contained objects */
    vector<string> list() const;

    /** Returns true if the map contains an entry under this name */
    bool contains(const string& name) const;

    /** Retrieves the smart pointer by name, return by value of ptr_type */
    ptr_type get(const string& name) const;

    /** Stores the smart pointer to object using the passed-in name 
        Returns the name and the version number
        pass by value of ptr_type
    */
    pair<string, unsigned long> set(const string& name, ptr_type sp);

    /** Stores the raw pointer to object using the passed-in name 
        Returns the name and the version number
    */
    pair<string, unsigned long> set(const string& name, T* p);

    /** Returns the version of the pointed object */
    unsigned long version(const string& name) const;

    /** Clears the map and resets the current version to 0 */
    void clear();
    
private:
    
    // Removes leading and trailing blanks and upper cases the passed in string.
    // Throws an exception if the string has internal blanks.
    string processName(const string& name) const;

    // Helper function to retrieve the pointer and the version number
    // Call processName() before calling this.
    pair_type get_pair(const string& name) const;

    // Helper function to remove the pointer and the version number;
    // returns true if the pair was found and removed.
    // Call processName() before calling this.
    bool remove_pair(const string& name);
    
    // state
    unsigned long ver_;        // the current version of the pointed object
};


///////////////////////////////////////////////////////////////////////////////
// Inline definitions

template<typename T>
inline string
SPtrMap<T>::processName(const string& name) const {
    string ret = trim(name);
    ASSERT(!ret.empty(), "empty object names not allowed");
    transform(ret.begin(), ret.end(), ret.begin(), ::toupper);
    auto b = find_if(ret.cbegin(), ret.cend(), ::isspace);
    ASSERT(b == ret.cend(), "blanks not allowed in object names");
    return ret;
}

template<typename T>
inline typename SPtrMap<T>::pair_type
SPtrMap<T>::get_pair(const string& name) const {
    auto it = map_type::find(name.c_str());
    return (it == map_type::end()? pair_type() : it-> second);
}

template<typename T>
inline bool
SPtrMap<T>::remove_pair(const string& name) {
     auto it = map_type::find(name.c_str());
     if (it == map_type::end())
         return false;
     else {
         map_type::erase(it);
     }
     return true;
}

template<typename T>
inline vector<string>
SPtrMap<T>::list() const {
    vector<string> keys;
    for (auto it = map_type::cbegin(); it != map_type::cend(); ++it)
        keys.push_back(it->first);
    return keys;
}

template<typename T>
inline bool
SPtrMap<T>::contains(const string& name) const {
    string nm = processName(name);
    auto it = map_type::find(nm.c_str());
    return it != map_type::end();
}

template<typename T>
inline typename SPtrMap<T>::ptr_type
SPtrMap<T>::get(const string& name) const {
    string nm = processName(name);
    pair_type pr = get_pair(nm);
    return pr.first;
}

template<typename T>
inline pair<string, unsigned long>
SPtrMap<T>::set(const string& name, ptr_type sp) {
     string nm = processName(name);
     bool b = remove_pair(nm);
     ver_++;
     map_type::insert(make_pair(nm, make_pair(sp, ver_)));
     return make_pair(nm, ver_);
}

template<typename T>
inline pair<string, unsigned long>
SPtrMap<T>::set(const string& name, T* p) {
    return set(name, ptr_type(p));   // p will be deleted after this function call
}

template<typename T>
inline unsigned long
SPtrMap<T>::version(const string& name) const {
    string nm = processName(name);
    auto pr = get_pair(nm);
    return pr.second;
}

template<typename T>
inline void SPtrMap<T>::clear() {
    map_type::clear();
    ver_ = 0;
    return;
}

#endif // SPTRMAP_HPP

