cmake_minimum_required (VERSION 2.6)

#
# проверка Unicode (UTF-8):
# 偉大和強大的俄語。
# اللغة الروسية العظمى والأقوياء.
# Η μεγάλη και πανίσχυρη ρωσική γλώσσα.
# השפה הרוסיה הגדולה והעצומה.
# 大いなるロシア語。
# 위대하고 강력한 러시아어 언어입니다.
# Wielki i potężny język rosyjski.
#

set( PROJECT_NAME "CTChestBoneSegmenter" )
project ( ${PROJECT_NAME} )

###############################################
#set(CMAKE_BUILD_TYPE "Debug" )
set(CMAKE_BUILD_TYPE "Release" )
if(NOT CMAKE_BUILD_TYPE)
#    set(CMAKE_BUILD_TYPE "Release")
    set(CMAKE_BUILD_TYPE "Debug")
endif()

set(PLATFORM "x32")
if("${CMAKE_SIZEOF_VOID_P}" STREQUAL "8")
    set(PLATFORM "x64")
endif()

set(BUILD_PREFIX "${CMAKE_BUILD_TYPE}_${PLATFORM}")
message("build prefix=[${BUILD_PREFIX}]")

###############################################
FIND_PACKAGE( OpenCV REQUIRED )

INCLUDE_DIRECTORIES(
	${PROJECT_SOURCE_DIR}
	${PROJECT_SOURCE_DIR}/src
	${PROJECT_SOURCE_DIR}/src/analyzecv
)

###############################################
set( TARGET_1		"CTChestSegmenterAnalz" )
set( TARGET_2		"CTChestSegmenterNifti" )

###############################################
set( ALL_SRC
        ${PROJECT_SOURCE_DIR}/src/semicommon.cpp
        ${PROJECT_SOURCE_DIR}/src/chsegm.cpp
        ${PROJECT_SOURCE_DIR}/src/fft.cpp
        ${PROJECT_SOURCE_DIR}/src/analyzecv/analyzecv.cpp
        ${PROJECT_SOURCE_DIR}/src/common.cpp
        ${PROJECT_SOURCE_DIR}/src/srfio.cpp
)

set( ALL_HDR
        ${PROJECT_SOURCE_DIR}/src/fft.h
        ${PROJECT_SOURCE_DIR}/src/srfio.h
        ${PROJECT_SOURCE_DIR}/src/semicommon.h
        ${PROJECT_SOURCE_DIR}/src/analyzecv/analyzecv.h
        ${PROJECT_SOURCE_DIR}/src/stdafx.h
        ${PROJECT_SOURCE_DIR}/src/common.h
        ${PROJECT_SOURCE_DIR}/src/chsegm.h
)


set( TARGET_1_SRC
        ${PROJECT_SOURCE_DIR}/CTChestSegmenter.cpp
)

set( TARGET_2_SRC
        ${PROJECT_SOURCE_DIR}/CTChestSegmenterNifti.cpp
)

###### Find Nifti
if (UNIX)
    include_directories(/usr/include/nifti)
endif()

set(LIB_NIFTI
        niftiio
        znz
)

###############################################
ADD_EXECUTABLE("${TARGET_1}_${BUILD_PREFIX}" ${ALL_SRC} ${ALL_HDR} ${TARGET_1_SRC})
TARGET_LINK_LIBRARIES("${TARGET_1}_${BUILD_PREFIX}" ${OpenCV_LIBS})

ADD_EXECUTABLE("${TARGET_2}_${BUILD_PREFIX}" ${ALL_SRC} ${ALL_HDR} ${TARGET_2_SRC})
TARGET_LINK_LIBRARIES("${TARGET_2}_${BUILD_PREFIX}" ${LIB_NIFTI} ${OpenCV_LIBS})

##get_cmake_property(_variableNames VARIABLES)
##foreach (_variableName ${_variableNames})
##    message(STATUS "${_variableName}=${${_variableName}}")
##endforeach()
