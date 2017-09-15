PROJECT (RunLinux)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

#find_package(Eigen3)
#include_directories (${EIGEN3_INCLUDE_DIR})

if(APPLE)
  include_directories(/opt/local/include/)
  include_directories(/opt/local/include/eigen3)
  include_directories("${PROJECT_SOURCE_DIR}/../include/osx")
  ADD_DEFINITIONS(--std=c++0x)
elseif(UNIX)
  include_directories("${PROJECT_SOURCE_DIR}/../include/unix")
  include_directories("${PROJECT_SOURCE_DIR}/../include/")
  ADD_DEFINITIONS(--std=c++0x)
elseif(WIN32)
	include_directories("${PROJECT_SOURCE_DIR}/../include/")
	LINK_DIRECTORIES("$(ProjectDir)/../../libs")
endif()

set(PARDISO_AVAILABLE ON)

if (PARDISO_AVAILABLE)
  add_definitions(-DPARDISO_AVAILABLE)
  message("Has pardiso")
  if(APPLE)
    include_directories(/opt/intel/mkl/include/)
    LINK_DIRECTORIES(/opt/intel/mkl/lib/)
    LINK_DIRECTORIES(/opt/intel/lib/)
  elseif(UNIX)
    LINK_DIRECTORIES("${PROJECT_SOURCE_DIR}/../libs")
  elseif(WIN32)
    link_directories($ENV{PARDISO_LIB_PATH})
	message($ENV{PARDISO_LIB_PATH})
   endif()
endif()

set(PROFILING "false" CACHE BOOL "Chosen By User")
if(${PROFILING} MATCHES "true")
  add_definitions("-pg")
endif()

add_subdirectory(../femlib       "${CMAKE_CURRENT_BINARY_DIR}/femlib")
add_subdirectory(../linsolver   "${CMAKE_CURRENT_BINARY_DIR}/linsolver")
add_subdirectory(../StepperCPU   "${CMAKE_CURRENT_BINARY_DIR}/StepperCPU")
add_subdirectory(../mccd       "${CMAKE_CURRENT_BINARY_DIR}/mccd")

FILE(GLOB EXE_SRC "*.cpp")

INCLUDE_DIRECTORIES ("${PROJECT_SOURCE_DIR}/../femlib/include")
INCLUDE_DIRECTORIES ("${PROJECT_SOURCE_DIR}/../StepperCPU/include")
INCLUDE_DIRECTORIES ("${PROJECT_SOURCE_DIR}/../linsolver/include/pcgsolver")
INCLUDE_DIRECTORIES ("${PROJECT_SOURCE_DIR}/../linsolver/include")
INCLUDE_DIRECTORIES ("${PROJECT_SOURCE_DIR}/../MCCD/include")

INCLUDE_DIRECTORIES(include)

if(APPLE)
#  link_directories("${PROJECT_SOURCE_DIR}/../lib/osx")
  link_directories("/opt/local/lib")
endif()

ADD_EXECUTABLE(Run ${EXE_SRC})

TARGET_LINK_LIBRARIES(Run StepperCPU)
TARGET_LINK_LIBRARIES(Run FEMLib)
TARGET_LINK_LIBRARIES(Run linsolver)
TARGET_LINK_LIBRARIES(Run mccd)

if (PARDISO_AVAILABLE)
  if(APPLE)
    TARGET_LINK_LIBRARIES(Run mkl_intel mkl_intel_thread mkl_core)
  elseif(UNIX)
    TARGET_LINK_LIBRARIES(Run pardiso openblas gfortran)
  elseif(WIN32)
	target_link_libraries(Run libpardiso500-WIN-X86-64)
  endif()
endif()

if(APPLE)
  target_link_libraries(Run ${COCOA_LIBRARY} ${OpenGL_LIBRARY} ${IOKIT_LIBRARY} ${COREVIDEO_LIBRARY})
elseif(UNIX)
  TARGET_LINK_LIBRARIES(Run gomp)
  TARGET_LINK_LIBRARIES(Run pthread)
elseif(WIN32)
  add_custom_command(TARGET Run POST_BUILD
					   COMMAND set "var=${PROJECT_SOURCE_DIR}"
					   COMMAND set "var=%var:/=\\%"
					   COMMAND set "var1=${PROJECT_BINARY_DIR}"
					   COMMAND set "var1=%var1:/=\\%"
					   COMMAND copy "%var%\\..\\dll\\*.dll" "%var1%\\$<CONFIGURATION>" 
					   COMMENT "copying dlls ..."
				   )
endif()

if(${PROFILING} MATCHES "true")
  SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -pg")
endif()
