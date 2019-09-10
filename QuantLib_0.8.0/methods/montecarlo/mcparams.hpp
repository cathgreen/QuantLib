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
      RANLUX4,
      SOBOL
    };
    
  /** The known path generator types */
  enum class PathGenType
    {
      EULER
    };

  /** Control variate types */
  enum class ControlVarType
    {
      NONE,
      ANTITHETIC
    };
  

  /** Default ctor */
  McParams(UrngType u = UrngType::MT19937,
           PathGenType p = PathGenType::EULER,
           ControlVarType c = ControlVarType::NONE)
    : urngType(u), pathGenType(p), controlVarType(c) {}

  // state
  UrngType urngType;
  PathGenType pathGenType;
  ControlVarType controlVarType;
};


#endif // MCPARAMS_HPP
