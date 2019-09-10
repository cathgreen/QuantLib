cmake_minimum_required (VERSION 3.14)

set (CMAKE_CXX_STANDARD 11)
project (QuantLib)

# include_directories(
# math/stats
# math/interpol
# pricer)

find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})


add_executable(test 
main.cpp
math/stats/errorfunction.cpp 
math/interpol/piecewisepolynomial.cpp
pricer/simplepricers.cpp
market/yieldcurve.cpp
market/volatilitytermstructure.cpp
market/market.cpp
pricer/bsmcpricer.cpp
math/random/sobolurng.cpp 
)

target_link_libraries(test ${ARMADILLO_LIBRARIES})


if(MSVC)
  target_compile_options(test PRIVATE /W4 /WX)
else()
  # target_compile_options(test PRIVATE -Wall -Wextra -pedantic -Werror)
   target_compile_options(test PRIVATE -Wall -Wextra -pedantic)
endif()

add_custom_command(
    OUTPUT ${CMAKE_SOURCE_DIR}/dependencies/library/lib.o
    WORKING_DIR ${CMAKE_SOURCE_DIR}/dependencies/library
    COMMAND make
)
