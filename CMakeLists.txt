cmake_minimum_required (VERSION 3.14)

set (CMAKE_CXX_STANDARD 11)
project (TestQuantLib)
 
#For the shared library:
set (VER 9)
string(APPEND SharedLib "libQuantLib0" ${VER} "0.dylib")
string(APPEND LibDir "QuantLib_0." ${VER} ".0")
string(APPEND VER_TAR "0" ${VER} "0")


option (DIRECT_BUILD ON)

if (DIRECT_BUILD)
  if (DIRECT_BUILD STREQUAL "OFF")
	set ( PROJECT_LINK_LIBS ${SharedLib} )
	link_directories( ${CMAKE_CURRENT_LIST_DIR}/${LibDir}/build )	
  else()
	set ( PROJECT_LINK_LIBS "")
  endif()		
endif (DIRECT_BUILD)


include_directories(${CMAKE_CURRENT_LIST_DIR}/${LibDir})

find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})


if (DIRECT_BUILD)
  if (DIRECT_BUILD STREQUAL "OFF")
	set(MySourcesList "")
  else() 
	set(MySourcesList ${CMAKE_CURRENT_LIST_DIR}/${LibDir}/math/stats/errorfunction.cpp 
			${CMAKE_CURRENT_LIST_DIR}/${LibDir}/math/interpol/piecewisepolynomial.cpp
			${CMAKE_CURRENT_LIST_DIR}/${LibDir}/pricer/simplepricers.cpp
			${CMAKE_CURRENT_LIST_DIR}/${LibDir}/market/yieldcurve.cpp
			${CMAKE_CURRENT_LIST_DIR}/${LibDir}/market/volatilitytermstructure.cpp
			${CMAKE_CURRENT_LIST_DIR}/${LibDir}/market/market.cpp
			${CMAKE_CURRENT_LIST_DIR}/${LibDir}/pricer/bsmcpricer.cpp
			${CMAKE_CURRENT_LIST_DIR}/${LibDir}/pricer/multiassetbsmcpricer.cpp
			${CMAKE_CURRENT_LIST_DIR}/${LibDir}/math/random/sobolurng.cpp 
			${CMAKE_CURRENT_LIST_DIR}/${LibDir}/methods/montecarlo/pathgenerator.cpp
			${CMAKE_CURRENT_LIST_DIR}/${LibDir}/math/linalg/choldcmp.cpp
			${CMAKE_CURRENT_LIST_DIR}/${LibDir}/math/linalg/eigensym.cpp
			${CMAKE_CURRENT_LIST_DIR}/${LibDir}/math/linalg/spectrunc.cpp
			${CMAKE_CURRENT_LIST_DIR}/${LibDir}/methods/pde/pdebase.cpp
			${CMAKE_CURRENT_LIST_DIR}/${LibDir}/methods/pde/pde1dsolver.cpp
			)
  endif()
endif (DIRECT_BUILD)

add_executable(test testQL.cpp ${MySourcesList})
target_link_libraries(test ${PROJECT_LINK_LIBS} ${ARMADILLO_LIBRARIES})


if(MSVC)
  target_compile_options(test PRIVATE /W4 /WX)
else()
   target_compile_options(test PRIVATE -Wall -Wextra -pedantic)
endif()

set_target_properties(test
        PROPERTIES SUFFIX ${VER_TAR}
        )
