message("${BoldGreen}----------------------------")
message("Preparing MonteCarlo libraries...")
message("----------------------------${ColourReset}")

include_directories(${CMAKE_SOURCE_DIR})

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#

option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif(WITH_GEANT4_UIVIS)

message("Geant4 use file variable is pointing at: ${Geant4_USE_FILE}")
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project


# Find ROOT
find_package (ROOT REQUIRED)
include_directories (${ROOT_INCLUDE_DIR})

#Xerces-C support
find_package (Xerces REQUIRED)
include_directories (${XERCESC_INCLUDE})

#CRMC support
if (USE_CRMC)
  find_package (CRMC)
  if (CRMC_FOUND)
      add_compile_definitions( CRMC='true' )
      include_directories (${CRMC_INCLUDE_DIRS})
      link_directories(${CRMC_INCLUDE_DIRS})
  else (NOT CRMC_FOUND)
    message("${BoldRed} CRMC library was not found! Please set CRMC_SRC and CRMC_INSTALL environment variables. ${ColourReset}")
  endif (CRMC_FOUND)
endif(USE_CRMC)


# Add ROOTs header paths

#message(${CMAKE_MODULE_PATH} )

option(ANALYSIS_HH_ "Build with ROOT" 1)
add_definitions(-DG4UI_USE -DG4VIS_USE -DANALYSIS_HH_)


if(useROOT)
	EXECUTE_PROCESS(COMMAND root-config --cflags OUTPUT_VARIABLE ROOT_CXX_FLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}")
endif(useROOT)


include_directories(${PROJECT_SOURCE_DIR}/MonteCarlo/include ${ROOT_INCLUDE_DIR})


#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/MonteCarlo/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/MonteCarlo/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(zdc zdc.cc ${sources} ${headers})

if(useROOT)
	EXECUTE_PROCESS(COMMAND root-config --libs OUTPUT_VARIABLE ROOT_LD_FLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
endif(useROOT)

target_link_libraries(zdc ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} ${XERCESC_LIBRARY} )

if (USE_CRMC AND CRMC_FOUND)
    target_link_libraries(zdc  ${CRMC_LIBRARIES} )
endif ()

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build . This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(EXAMPLE_SCRIPTS
  #zdc.in
  #zdc.out
  init_vis.mac
  run1.mac
  vis.mac
  geometry.mac
  )

foreach(_script ${EXAMPLE_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/MonteCarlo/${_script}
    ${PROJECT_BINARY_DIR}/MonteCarlo/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# For internal Geant4 use - but has no effect if you build this
# example standalone
#
add_custom_target( DEPENDS zdc)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#

if (USE_CRMC AND CRMC_FOUND)
    install(FILES $ENV{CRMC_INSTALL}/crmc.param DESTINATION bin)
endif ()
install(TARGETS zdc DESTINATION bin)
