CMAKE_MINIMUM_REQUIRED (VERSION 2.6)
PROJECT (libdai)

# The version number.
SET (libdai_VERSION "0.3.0")

STRING (REGEX MATCHALL "[0-9]" libdai_VERSION_PARTS "${libdai_VERSION}")

LIST (GET libdai_VERSION_PARTS 0 libdai_VERSION_MAJOR)
LIST (GET libdai_VERSION_PARTS 1 libdai_VERSION_MINOR)
LIST (GET libdai_VERSION_PARTS 2 libdai_VERSION_PATCH)

SET (libdai_SOVERSION "${libdai_VERSION_MAJOR}.${libdai_VERSION_MINOR}")

SET(BUILD_SHARED_LIBS ON)

IF(CMAKE_COMPILER_IS_GNUCXX)
  SET (EXTRA_C_FLAGS "-Wno-deprecated -Wall -W -Wextra -fPIC")
  SET (EXTRA_C_FLAGS_RELEASE "-O3")
  SET (EXTRA_C_FLAGS_DEBUG "-O3 -g -DDAI_DEBUG")

  IF(APPLE)
      ADD_DEFINITIONS(-DMACOSX)
  ENDIF(APPLE)
  IF(CYGWIN)
    ADD_DEFINITIONS(-DCYGWIN)
  ENDIF(CYGWIN)
ENDIF()
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${EXTRA_C_FLAGS}")
SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${EXTRA_C_FLAGS_RELEASE}")
SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${EXTRA_C_FLAGS_DEBUG}")

# Searching for additional libraries
SET (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMakeModules/;${CMAKE_MODULE_PATH}")

FIND_PACKAGE(Boost COMPONENTS REQUIRED program_options unit_test_framework)
FIND_PACKAGE(GMP REQUIRED)

# Add include directories
INCLUDE_DIRECTORIES ("${CMAKE_SOURCE_DIR}/include")
INCLUDE_DIRECTORIES ("${Boost_INCLUDE_DIRS}")
INCLUDE_DIRECTORIES ("${GMP_INCLUDE_DIRS}")

# Add libraries
SET (EXTRA_LIBS "${GMP_LIBRARIES}")

SET(WITH_BP ON CACHE BOOL "Build with Belief Propagation (BP)")
SET(WITH_FBP ON CACHE BOOL "Build with Fractional Belief Propagation (FBP)")
SET(WITH_TRWBP ON CACHE BOOL "Build with Tree-Reweighted Belief Propagation (TRWBP)")
SET(WITH_MF ON CACHE BOOL "Build with Mean Field (MF)")
SET(WITH_HAK ON CACHE BOOL "Build with algorithms by Heskes, Albers and Kappen (HAK)")
SET(WITH_LC ON CACHE BOOL "Build with loop corrections for approximate inference (LC)")
SET(WITH_TREEEP ON CACHE BOOL "Build with Tree Expectation Propagation (TREEEP)")
SET(WITH_JTREE ON CACHE BOOL "Build with Junction Tree algorithm (JT)")
SET(WITH_MR ON CACHE BOOL "Build with loop corrections as proposed by Montanari and Rizzo (MR)")
SET(WITH_GIBBS ON CACHE BOOL "Build with Gibbs sampling (GIBBS)")
SET(WITH_CBP ON CACHE BOOL "Build with Conditioned Belief Propagation (CBT)")
SET(WITH_DECMAP ON CACHE BOOL "Build with MAP state by decimation (DECMAP)")

SET(WITH_MATLAB OFF CACHE BOOL "Build MatLab interface? (MatLab needs to be installed)")
SET(WITH_CIMG OFF CACHE BOOL "Build image segmentation example? (CImg needs to be installed)")
IF(WITH_CIMG)
    SET(NEW_CIMG ON CACHE BOOL "CImg version 1.3.0 or newer?")
ENDIF(WITH_CIMG)

SET(WITH_TESTS ON CACHE BOOL "Build and run tests")

IF(WITH_BP)
    ADD_DEFINITIONS ("-DDAI_WITH_BP")
ENDIF(WITH_BP)
IF(WITH_FBP)
    ADD_DEFINITIONS ("-DDAI_WITH_FBP")
ENDIF(WITH_FBP)
IF(WITH_TRWBP)
    ADD_DEFINITIONS ("-DDAI_WITH_TRWBP")
ENDIF(WITH_TRWBP)
IF(WITH_MF)
    ADD_DEFINITIONS ("-DDAI_WITH_MF")
ENDIF(WITH_MF)
IF(WITH_HAK)
    ADD_DEFINITIONS ("-DDAI_WITH_HAK")
ENDIF(WITH_HAK)
IF(WITH_LC)
    ADD_DEFINITIONS ("-DDAI_WITH_LC")
ENDIF(WITH_LC)
IF(WITH_TREEEP)
    ADD_DEFINITIONS ("-DDAI_WITH_TREEEP")
ENDIF(WITH_TREEEP)
IF(WITH_JTREE)
    ADD_DEFINITIONS ("-DDAI_WITH_JTREE")
ENDIF(WITH_JTREE)
IF(WITH_MR)
    ADD_DEFINITIONS ("-DDAI_WITH_MR")
ENDIF(WITH_MR)
IF(WITH_GIBBS)
    ADD_DEFINITIONS ("-DDAI_WITH_GIBBS")
ENDIF(WITH_GIBBS)
IF(WITH_CBP)
    ADD_DEFINITIONS ("-DDAI_WITH_CBP")
ENDIF(WITH_CBP)
IF(WITH_DECMAP)
    ADD_DEFINITIONS ("-DDAI_WITH_DECMAP")
ENDIF(WITH_DECMAP)

IF (WITH_TESTS)
    ENABLE_TESTING()
ENDIF()

# Build target lib
ADD_SUBDIRECTORY (src)
ADD_SUBDIRECTORY (include/dai)

SET (EXTRA_LIBS dai)
# SET (EXTRA_LIBS ${EXTRA_LIBS} dai)

# Build target tests
ADD_SUBDIRECTORY (tests)

# Build target utils
ADD_SUBDIRECTORY (utils)

# Build target examples
ADD_SUBDIRECTORY (examples)

# TODOs:
# * note that matlab is not working out of the box with CMake

# add swig
# add python

FIND_PACKAGE(Doxygen)
IF(DOXYGEN_FOUND)
    SET(WITH_DOC OFF CACHE BOOL "Generate HTML docs using Doxygen")
ELSE(DOXYGEN_FOUND)
    MESSAGE(STATUS "Doxygen not found. Cannot generate HTML docs using Doxygen. Please get a copy http://www.doxygen.org")
ENDIF(DOXYGEN_FOUND)

IF(WITH_DOC AND DOXYGEN_FOUND)
    # prepare doxygen configuration file
    CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

    # add doxygen as target
    ADD_CUSTOM_TARGET(doc ALL ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen" VERBATIM
    )

    # cleanup $build/api-doc on "make clean"
    SET_PROPERTY(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES api-doc)

    # install HTML API documentation and manual pages
    SET(DOC_PATH "share/doc/libdai-${libdai_VERSION}")

    INSTALL(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/api-doc/html
        DESTINATION ${DOC_PATH}
    )
ENDIF(WITH_DOC AND DOXYGEN_FOUND)

#README : doc scripts/makeREADME Makefile
#	DAI_VERSION=$(DAI_VERSION) DAI_DATE=$(DAI_DATE) scripts/makeREADME
#
#TAGS :
#	etags src/*.cpp include/dai/*.h tests/*.cpp utils/*.cpp
#	ctags src/*.cpp include/dai/*.h tests/*.cpp utils/*.cpp
