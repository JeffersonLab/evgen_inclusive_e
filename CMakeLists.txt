#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(InclusiveElectron CXX Fortran)
set(EXE "evgen_inclusive_e") 

#----------------------------------------------------------------------------
# Setup GNUInstallDirs for installation
include(GNUInstallDirs)

#----------------------------------------------------------------------------
# The following lines are required so that the executable can dynamically locate 
# the library dependencies independent of use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_FULL_LIBDIR}" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}")
ENDIF("${isSystemDir}" STREQUAL "-1")

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

#----------------------------------------------------------------------------
# Force use of c++11 standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#----------------------------------------------------------------------------
# ROOT path
# Need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
# - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
# - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED)

# Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

#----------------------------------------------------------------------
# get LHAPDF6 LIB and DIR
set(LHAPDF "$ENV{LHAPDF}")
link_directories(${LHAPDF}/lib)

#----------------------------------------------------------------------------
# add boost libraries
find_package(Boost REQUIRED) 

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS}) 
endif()

#----------------------------------------------------------------------------
# debug flags
# Fortran 
# set(CMAKE_Fortran_FLAGS "-g -O2 -Wall -ffixed-line-length-none -x f77-cpp-input ")
# CXX 
# set(CMAKE_CXX_FLAGS "-g -Wall")
set(CMAKE_BUILD_TYPE=RelWithDebInfo) 

#----------------------------------------------------------------------------
# Generate dictionaries as well (so we can write objects to ROOT files) 
include_directories(
${ROOT_INCLUDE_DIR} 
${LHAPDF}/include 
${PROJECT_SOURCE_DIR}/include 
${PROJECT_SOURCE_DIR}/src/
${CMAKE_CURRENT_BINARY_DIR}/include 
)

#----------------------------------------------------------------------------
# copy input directory to the install location 
set(INPUT_FILES
input/parameters_d2n_He3.txt         
input/parameters_solid_PVDIS_LD2.txt 
input/parameters_solid_SIDIS_He3.txt
input/parameters_solid_JPsi_LH2.txt  
input/parameters_solid_PVDIS_LH2.txt 
input/parameters_solid_SIDIS_NH3.txt
)

foreach(_infile ${INPUT_FILES})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_infile}
    ${PROJECT_BINARY_DIR}/${_infile}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: HEADERS are included so they will show up in IDEs
file(GLOB SOURCES ${PROJECT_SOURCE_DIR}/src/*.C ${PROJECT_SOURCE_DIR}/src/*.f)
file(GLOB HEADERS ${PROJECT_SOURCE_DIR}/include/*.h)

#----------------------------------------------------------------------------
# Add the executable
add_executable(${EXE} evgen_inclusive_e.C ${SOURCES} ${HEADERS})
target_link_libraries(${EXE} ${ROOT_LIBRARIES} LHAPDF ${Boost_LIBRARIES})

#----------------------------------------------------------------------------
# Install items 
install(TARGETS ${EXE} DESTINATION ${INSTALL_DIR})
install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(DIRECTORY ${PROJECT_SOURCE_DIR}/input DESTINATION ${CMAKE_INSTALL_PREFIX})

