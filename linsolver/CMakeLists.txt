PROJECT (linsolver)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

FILE(GLOB LIB_SRC "src/*.cpp")
if(UNIX)
  ADD_DEFINITIONS(--std=c++0x)
endif (UNIX)
if (${PARDISO_AVAILABLE} MATCHES "true")
  add_definitions(-DPARDISO_AVAILABLE)
  message("Has pardiso")
endif()

INCLUDE_DIRECTORIES ("${PROJECT_SOURCE_DIR}/include")
INCLUDE_DIRECTORIES ("${PROJECT_SOURCE_DIR}/../femlib/include")
INCLUDE_DIRECTORIES ("${PROJECT_SOURCE_DIR}/../include")
INCLUDE_DIRECTORIES ("${PROJECT_SOURCE_DIR}/include/pcgsolver")
INCLUDE_DIRECTORIES ("${PROJECT_SOURCE_DIR}/../include/coin")

ADD_LIBRARY (linsolver  ${LIB_SRC})

SET(LIB_TYPE STATIC)	#set the lib type
