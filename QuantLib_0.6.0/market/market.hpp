#ifndef MARKET_HPP
#define MARKET_HPP

#include "yieldcurve.hpp"
#include "volatilitytermstructure.hpp"
#include "../sptrmap.hpp"

class Market {
public:

    /** Returns the unique instance */
    static Market& instance();

    /** Clears the market of all objects */
    void clear();

    /** Returns the yield curves map, for read-write */
    SPtrMap<YieldCurve>& yieldCurves() {
        return ycmap_;
    }

    /** Returns the volatility term structures map, for read-write */
    SPtrMap<VolatilityTermStructure>& volatilities() {
        return volmap_;
    }

private:

    /** allow private default ctor */
    Market() {}

    /** forbid copy ctor and copy-assignment operator */
    Market(const Market & rhs) = delete;
    Market& operator=(const Market & rhs) = delete;
    
    //state
    SPtrMap<YieldCurve> ycmap_;
    SPtrMap<VolatilityTermStructure> volmap_;
};

/** Free function returning the market singleton */
Market& market();

#endif // MARKET_HPP
