#ifndef RNG_HPP
#define RNG_HPP

#include "normalrng.hpp"
#include "sobolurng.hpp"

/** Linear congruential */
using NormalRngMinStdRand = NormalRng<std::minstd_rand>;

/** Mersenne Twister */
using NormalRngMt19937 = NormalRng<std::mt19937>;

/** RanLux level 3 */
using NormalRngRanLux3 = NormalRng<std::ranlux24>;

/** RanLux level 4 */
using NormalRngRanLux4 = NormalRng<std::ranlux48>;

/** Sobol */
using NormalRngSobol = NormalRng<SobolURng>;


#endif // RNG_HPP
