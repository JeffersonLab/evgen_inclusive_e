#----------------------------------------------------------------------------

# Setup the project
#
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(InclusiveElectron)

#------------------------------------------------------------------------
#ROOT path
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")
find_package(ROOT REQUIRED)

# this is to ensure we find the header we need
include_directories(${CMAKE_CURRENT_BINARY_DIR}/include)

#----------------------------------------------------------------------
# get LHAPDF6 LIB and DIR

####################### LHAPDF ##############################
#LHAPDF=/home/yxzhao/workarea/solid_g1_g5/LHAPDF6/installation
#LHAPDF_LIB=$(LHAPDF)/lib
#PDFLIBS=-L$(LHAPDF_LIB) -lLHAPDF

set(LHAPDF "$ENV{LHAPDF}")
#set(LHAPDF "/work/halla/solid/apps/LHAPDF/LHAPDF-6.1.6/") 

#debugging symbols, complains
set(CMAKE_CXX_FLAGS "-g -Wall")

#----------------------------------------------------------------------------
# Generate dictionaries as well (so we can write objects to ROOT files

include_directories(${LHAPDF}/include ${PROJECT_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/include ${PROJECT_SOURCE_DIR}/src/)

link_directories(${LHAPDF}/lib)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.C)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.h)

#----------------------------------------------------------------------------
# Add the executable
#
add_executable(solid_inclusive_e solid_inclusive_e.C ${sources} ${headers})
target_link_libraries(solid_inclusive_e ${ROOT_LIBRARIES} LHAPDF)





