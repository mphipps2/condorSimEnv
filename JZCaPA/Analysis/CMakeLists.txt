message("${BoldYellow}----------------------------")
message("Preparing Analysis libraries...")
message("----------------------------${ColourReset}")

include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/Analysis/include)

#In the following the list of source files to be linked
set (ANALYSIS_SRC_DIR "${CMAKE_SOURCE_DIR}/Analysis/src")

#### Get all header files
FILE(GLOB_RECURSE PA_INC
    "${CMAKE_SOURCE_DIR}/Analysis/include/*.hh"
)
FILE(GLOB_RECURSE DEPENDENCIES
    "${CMAKE_SOURCE_DIR}/Analysis/include/*.hh"
    "${CMAKE_SOURCE_DIR}/Analysis/src/*.cpp"
)

add_executable(AnalysisExample ${DEPENDENCIES} userFunctions/AnalysisExample.cpp)
TARGET_LINK_LIBRARIES(AnalysisExample ${ROOT_LIBRARIES})
TARGET_LINK_LIBRARIES(AnalysisExample ${XERCESC_LIBRARY})

install(TARGETS AnalysisExample RUNTIME DESTINATION bin)

add_executable(rpdMLTrainingAnalysis ${DEPENDENCIES} userFunctions/rpdMLTrainingAnalysis.cpp)
TARGET_LINK_LIBRARIES(rpdMLTrainingAnalysis ${ROOT_LIBRARIES})
TARGET_LINK_LIBRARIES(rpdMLTrainingAnalysis ${XERCESC_LIBRARY})

install(TARGETS rpdMLTrainingAnalysis RUNTIME DESTINATION bin)

add_executable(myAnalysis ${DEPENDENCIES} userFunctions/myAnalysis.cpp)
TARGET_LINK_LIBRARIES(myAnalysis ${ROOT_LIBRARIES})
TARGET_LINK_LIBRARIES(myAnalysis ${XERCESC_LIBRARY})

install(TARGETS myAnalysis RUNTIME DESTINATION bin)


add_executable(postAnalysis ${DEPENDENCIES} userFunctions/postAnalysis.cpp)
TARGET_LINK_LIBRARIES(postAnalysis ${ROOT_LIBRARIES})
TARGET_LINK_LIBRARIES(postAnalysis ${XERCESC_LIBRARY})

install(TARGETS postAnalysis RUNTIME DESTINATION bin)

#add_executable(postAnalysisZDC ${DEPENDENCIES} userFunctions/postAnalysisZDC.cpp)
#TARGET_LINK_LIBRARIES(postAnalysisZDC ${ROOT_LIBRARIES})
#TARGET_LINK_LIBRARIES(postAnalysisZDC ${XERCESC_LIBRARY})

#install(TARGETS postAnalysisZDC RUNTIME DESTINATION bin)

#add_executable(postAnalysisRPD ${DEPENDENCIES} userFunctions/postAnalysisRPD.cpp)
#TARGET_LINK_LIBRARIES(postAnalysisRPD ${ROOT_LIBRARIES})
#TARGET_LINK_LIBRARIES(postAnalysisRPD ${XERCESC_LIBRARY})

#install(TARGETS postAnalysisRPD RUNTIME DESTINATION bin)

add_executable(RPD_Matrix_Analysis ${DEPENDENCIES} userFunctions/RPD_Matrix_Analysis.cpp)
TARGET_LINK_LIBRARIES(RPD_Matrix_Analysis ${ROOT_LIBRARIES})
TARGET_LINK_LIBRARIES(RPD_Matrix_Analysis ${XERCESC_LIBRARY})

install(TARGETS RPD_Matrix_Analysis RUNTIME DESTINATION bin)
