#[=======================================================================[.rst:
FindCRMC
-------

Finds the CRMC library.

Imported Targets
^^^^^^^^^^^^^^^^

This module provides the following imported targets, if found:

``CRMC::CRMC``
The CRMC library

Result Variables
^^^^^^^^^^^^^^^^

This will define the following variables:

``CRMC_FOUND``
True if the system has the CRMC library.
``CRMC_VERSION``
The version of the CRMC library which was found.
``CRMC_INCLUDE_DIR``
Include directories needed to use CRMC.
``CRMC_LIBRARIES``
Libraries needed to link to CRMC.


#]=======================================================================]

MESSAGE(STATUS "Looking for CRMC...")

SET(CRMC_INCLUDE_DIRS)

FIND_PATH(CRMC_BUILD_INCLUDE_DIR NAMES CRMCconfig.h
 PATHS $ENV{CRMC_INSTALL}
 PATH_SUFFIXES src
)

list( APPEND CRMC_INCLUDE_DIRS ${CRMC_BUILD_INCLUDE_DIR} )


FIND_PATH(CRMC_SRC_INCLUDE_DIR NAMES CRMC.h
 PATHS $ENV{CRMC_SRC}
 PATH_SUFFIXES src
)

list( APPEND CRMC_INCLUDE_DIRS ${CRMC_SRC_INCLUDE_DIR} )


FIND_PATH(CRMC_LIBRARY_DIR NAMES libCrmcBasic.so
  PATHS $ENV{CRMC_INSTALL}
  PATH_SUFFIXES lib
)

IF(CRMC_LIBRARY_DIR)
  SET(CRMC_LIBRARIES "-L${CRMC_LIBRARY_DIR} -lCrmcBasic -lHepEvtDummy${DEFINED_MODELS}")
ELSE()
  SET(CRMC_LIBRARIES "")
ENDIF(CRMC_LIBRARY_DIR)

FILE(READ ${CRMC_BUILD_INCLUDE_DIR}/CRMCconfig.h CRMCCFGH)

STRING(REGEX MATCHALL "\n *#define CRMC_VERSION_MAJOR +[0-9]+" CRMCVERMAJ
  ${CRMCCFGH})
STRING(REGEX MATCH "\n *#define CRMC_VERSION_MINOR +[0-9]+" CRMCVERMIN
  ${CRMCCFGH})

STRING(REGEX REPLACE "\n *#define CRMC_VERSION_MAJOR +" ""
  CRMCVERMAJ ${CRMCVERMAJ})
STRING(REGEX REPLACE "\n *#define CRMC_VERSION_MINOR +" ""
  CRMCVERMIN ${CRMCVERMIN})

SET(CRMC_VERSION ${CRMCVERMAJ}.${CRMCVERMIN})

MESSAGE(STATUS "Looking for CRMC... - found ${CRMC_INCLUDE_DIRS}")
MESSAGE(STATUS "Looking for CRMC... - found ${CRMC_LIBRARIES}")
MESSAGE(STATUS "Looking for CRMC... - version ${CRMC_VERSION} ")

IF (CRMC_INCLUDE_DIRS AND CRMC_LIBRARIES)
SET(CRMC_FOUND TRUE)
ELSE (CRMC_INCLUDE_DIRS AND CRMC_LIBRARIES)
SET(CRMC_FOUND FALSE)
ENDIF (CRMC_INCLUDE_DIRS AND CRMC_LIBRARIES)


include( FindPackageHandleStandardArgs )
find_package_handle_standard_args( CRMC DEFAULT_MSG CRMC_INCLUDE_DIRS
   CRMC_LIBRARIES )
