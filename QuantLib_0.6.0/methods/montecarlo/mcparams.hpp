#ifndef MCPARAMS_HPP
#define MCPARAMS_HPP

#include "../../defines.hpp"
#include "../../exception.hpp"

struct McParams {
    /** The known URNG types */
    enum class UrngType
        {
             MINSTDRAND,
             MT19937,
             RANLUX3,
             RANLUX4
        };
    
    /** The known path generator types */
    enum class PathGenType
        {
            EULER
        };

    /** Default ctor */
    McParams(UrngType u = UrngType::MT19937, PathGenType p = PathGenType::EULER)
        : urngType(u), pathGenType(p) {}

    // state
    UrngType urngType;
    PathGenType pathGenType;
};


#endif // MCPARAMS_HPP
