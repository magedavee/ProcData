#----------------------------------------------------------------------------
# Compile command, from inside build directory:
# cmake ../PROSPECT-G4/
# make -j4

#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(PROSPECT-G4)

add_compile_options(-std=c++11)

#----------------------------------------------------------------------------
#For Mac OS X Use @rpath on OS X
if(APPLE) 
  if(POLICY CMP0042)
      CMAKE_POLICY(SET CMP0042 NEW)
   endif()
endif()

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc ${PROJECT_SOURCE_DIR}/src/*/*.cc )
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh ${PROJECT_SOURCE_DIR}/include/*/*.hh)
include_directories(${PROJECT_SOURCE_DIR}/include
                    ${PROJECT_SOURCE_DIR}/include/Physics
                    ${PROJECT_SOURCE_DIR}/include/Geometry
                    ${PROJECT_SOURCE_DIR}/include/EventGen
                    ${PROJECT_SOURCE_DIR}/include/Output
                    ${PROJECT_SOURCE_DIR}/include/Tracking
                    ${PROJECT_SOURCE_DIR}/include/aCORN
                    ${PROJECT_SOURCE_DIR}/include/FaNS
                    $ENV{CRYHOME}/src)
link_directories($ENV{CRYHOME}/lib)

#-------------------------------
# ROOT dependencies
#-------------------------------

######## locate ROOT install, headers, includes
include(FindROOT.cmake)
include_directories(${ROOT_INCLUDE_DIR})
link_directories(${ROOT_LIBRARY_DIR})

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()
include(${Geant4_USE_FILE})
Message("Geant4_USE_FILE = ${Geant4_USE_FILE}")

##############
# Output paths
##############

set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")


###############
# TObject subclasses dictionary & ROOT objects shared library
###############

set(SB_G4_LINKDEF EventLinkDef.h)
set(SB_G4_DICTIONARY ${CMAKE_CURRENT_BINARY_DIR}/EventLib_Dict.cxx)
set(DICTFILES "Event.hh")
ROOT_GENERATE_DICTIONARY("${SB_G4_DICTIONARY}" "${DICTFILES}" LINKDEF "${SB_G4_LINKDEF}" OPTIONS "")
LIST(APPEND sources ${SB_G4_DICTIONARY})

add_library(EventLib SHARED ${sources} ${headers})
target_link_libraries(EventLib CRY ${ROOT_LIBRARIES} ${Geant4_LIBRARIES})
SET(SB_G4_LIBRARY_PROPERTIES
    ${SB_G4_LIBRARY_PROPERTIES}
    VERSION "0"
    SOVERSION "0"
    SUFFIX ".so"
)
set_target_properties(EventLib PROPERTIES ${SB_G4_LIBRARY_PROPERTIES})
install(TARGETS EventLib DESTINATION ${LIBRARY_OUTPUT_PATH})
link_directories(${LIBRARY_OUTPUT_PATH})

#################
# Executables
#################

add_executable(PROSPECT-G4 PROSPECT-G4.cc)
target_link_libraries(PROSPECT-G4 EventLib ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})
if(WITH_HDF5)
    message("HDF5 compile dependencies enabled.")
    FIND_PACKAGE (HDF5 NAMES hdf5 hdf5_hl REQUIRED)
    if(HDF5_FOUND)
       message("Including HDF5 directories '${HDF5_INCLUDE_DIRS}' and libraries '${HDF5_LIBRARIES}'")
       include_directories(${HDF5_INCLUDE_DIRS})
       target_link_libraries (PROSPECT-G4 ${HDF5_LIBRARIES})
       target_link_libraries (EventLib ${HDF5_LIBRARIES})
    endif (HDF5_FOUND)
    add_definitions(-DWITH_HDF5)
endif()

install(TARGETS PROSPECT-G4 DESTINATION ${EXECUTABLE_OUTPUT_PATH})

#################
# Documentation with Doxygen
#################
find_package(Doxygen)
if(DOXYGEN_FOUND)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    add_custom_target(doc
        ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen" VERBATIM
    )
endif(DOXYGEN_FOUND)

#################
# Build settings
#################

#set(CMAKE_BUILD_TYPE release)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
set(CMAKE_BUILD_TYPE debug)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O3 -Wextra -Wpedantic")
