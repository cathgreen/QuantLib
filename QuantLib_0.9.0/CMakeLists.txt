cmake_minimum_required (VERSION 3.14)

set (CMAKE_CXX_STANDARD 11)
project (QuantLib)


find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

set(CMAKE_BUILD_TYPE Release)

set(MySourcesList 
math/stats/errorfunction.cpp 
math/interpol/piecewisepolynomial.cpp
pricer/simplepricers.cpp
market/yieldcurve.cpp
market/volatilitytermstructure.cpp
market/market.cpp
pricer/bsmcpricer.cpp
pricer/multiassetbsmcpricer.cpp
math/random/sobolurng.cpp 
methods/montecarlo/pathgenerator.cpp
math/linalg/choldcmp.cpp
math/linalg/eigensym.cpp
math/linalg/spectrunc.cpp
)

set (VER_TAR 090)
string(APPEND SUFFIX_TAR ${VER_TAR} ".dylib")

#Generate the shared library from the sources
add_library(QuantLib SHARED ${MySourcesList})

target_link_libraries(QuantLib ${ARMADILLO_LIBRARIES})

if(MSVC)
  target_compile_options(QuantLib PRIVATE /W4 /WX)
else()
  target_compile_options(QuantLib PRIVATE -Wall -Wextra -pedantic)
endif()

set_target_properties(QuantLib
        PROPERTIES SUFFIX ${SUFFIX_TAR}
        )

