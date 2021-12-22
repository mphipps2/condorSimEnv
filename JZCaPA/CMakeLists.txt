#Basic Cmake structure created by R.Longo on 12/13/2018
project(JZCaPA)
cmake_minimum_required(VERSION 2.8)
#aux_source_directory(. SRC_LIST)

#Colors definition for messages
if(NOT WIN32)
  string(ASCII 27 Esc)
  set(ColourReset "${Esc}[m")
  set(ColourBold  "${Esc}[1m")
  set(Red         "${Esc}[31m")
  set(Green       "${Esc}[32m")
  set(Yellow      "${Esc}[33m")
  set(Blue        "${Esc}[34m")
  set(White       "${Esc}[37m")
  set(BoldRed     "${Esc}[1;31m")
  set(BoldGreen   "${Esc}[1;32m")
  set(BoldYellow  "${Esc}[1;33m")
  set(BoldBlue    "${Esc}[1;34m")
  set(BoldWhite   "${Esc}[1;37m")
endif()

message("${BoldBlue}
       # #######  #####         ######     #
       #      #  #     #   ##   #     #   # #
       #     #   #        #  #  #     #  #   #
       #    #    #       #    # ######  #     #
 #     #   #     #       ###### #       #######
 #     #  #      #     # #    # #       #     #
  #####  #######  #####  #    # #       #     #
${ColourReset}")

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR})

#To be commented out on old linux...
add_compile_options(-std=c++11)
#add_compile_options(-arch x86_64)

if (${JZCaPA_ON_LXPLUS} MATCHES "YES")
message("${BoldWhite} LXPLUS modern compiler support enabled! ${ColourReset}")
add_compile_options(-std=c++1y)
endif (${JZCaPA_ON_LXPLUS} MATCHES "YES")

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmakeModules/)

#Root support
find_package (ROOT REQUIRED)
include_directories (${ROOT_INCLUDE_DIR})

#Xerces-C support
find_package (Xerces REQUIRED)
include_directories (${XERCESC_INCLUDE})



#Test executable - does nothing
add_executable(test_JZCaPA  main.cpp)
TARGET_LINK_LIBRARIES(test_JZCaPA ${ROOT_LIBRARIES})
TARGET_LINK_LIBRARIES(test_JZCaPA ${XERCESC_LIBRARY})



#Subdirectory structure
#Analysis libraries and user functions
add_subdirectory(Analysis)

#MonteCarlo folder. The installation will be conditional and disabled by default (not everyone needs to work with it and requires dependencies)
#The starting point will be the work done by Mike Phipps

if (${JZCaPA_ENABLE_MC} MATCHES "YES")
    message("${BoldGreen} Monte Carlo extension enabled! ${ColourReset}")
    add_subdirectory(MonteCarlo)
endif (${JZCaPA_ENABLE_MC} MATCHES "YES")

#Install test
install(TARGETS test_JZCaPA RUNTIME DESTINATION bin)
#Install doxygen files
install(FILES doxygen/JZCaPA_doxy.cnf DESTINATION doxygen)
install(FILES Utils/ConfigFile2018.xml DESTINATION Utils)
install(FILES Utils/Alignment_2018.xml DESTINATION Utils)
install(FILES Utils/Survey_2018.xml DESTINATION Utils)
install(FILES Utils/ConfigFileMC.xml DESTINATION Utils)
install(FILES Utils/Alignment_MC.xml DESTINATION Utils)
install(FILES Utils/PurgMag3D.TABLE DESTINATION Utils)
install(FILES Utils/eposPTneut.root DESTINATION Utils)
install(FILES Utils/ZDC_AvsC.root DESTINATION Utils)
install(FILES Analysis/ZDCBeamTestRun99.root DESTINATION bin)
install(FILES MonteCarlo/geometry.mac DESTINATION bin)
install(FILES MonteCarlo/run1.mac DESTINATION bin)
install(FILES MonteCarlo/init_vis.mac DESTINATION bin)
install(FILES MonteCarlo/vis.mac DESTINATION bin)

install(FILES Utils/jcapa_env_lxplus.sh DESTINATION Utils)
install(FILES Utils/jzcapa_env_uiuc.sh DESTINATION Utils)
install(DIRECTORY DESTINATION results)
install(DIRECTORY Utils/Timing_data DESTINATION Utils)
