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

#Generate the shared library from the sources
add_library(QuantLib080 SHARED ${MySourcesList})

#add_executable(test main.cpp ${MySourcesList})

target_link_libraries(QuantLib080 ${ARMADILLO_LIBRARIES})

if(MSVC)
  target_compile_options(QuantLib080 PRIVATE /W4 /WX)
else()
  # target_compile_options(QuantLib080 PRIVATE -Wall -Wextra -pedantic -Werror)
   target_compile_options(QuantLib080 PRIVATE -Wall -Wextra -pedantic)
endif()

# add_custom_command(
#    OUTPUT ${CMAKE_SOURCE_DIR}/dependencies/library/lib.o
#    WORKING_DIR ${CMAKE_SOURCE_DIR}/dependencies/library
#    COMMAND make
#)
