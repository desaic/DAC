PROJECT (RunOpt)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

#find_package(Eigen3)
#include_directories (${EIGEN3_INCLUDE_DIR})

if(APPLE)
  include_directories(/opt/local/include/)
  include_directories(/opt/local/include/eigen3)
  include_directories("${PROJECT_SOURCE_DIR}/../include/osx")
elseif(UNIX)
  include_directories("${PROJECT_SOURCE_DIR}/../include/unix")
  include_directories("${PROJECT_SOURCE_DIR}/../include/")
  include_directories("${PROJECT_SOURCE_DIR}/../include/lbfgs")
  include_directories("${PROJECT_SOURCE_DIR}/../include/coin")
elseif(WIN32)
    add_definitions("-DHAVE_CONFIG_H") #for ipopt
    include_directories("${PROJECT_SOURCE_DIR}/../include/win/coin")
	include_directories("${PROJECT_SOURCE_DIR}/../include/lbfgs")
	include_directories("${PROJECT_SOURCE_DIR}/../include/")
	LINK_DIRECTORIES("$(ProjectDir)/../../libs/ipopt")
endif()

set(PARDISO_AVAILABLE ON)

if (PARDISO_AVAILABLE)
  add_definitions(-DPARDISO_AVAILABLE)
  message("Has pardiso")
  if(APPLE)
    include_directories(/opt/intel/mkl/include/)
    LINK_DIRECTORIES(/opt/intel/mkl/lib/)
  elseif(UNIX)
    LINK_DIRECTORIES("${PROJECT_SOURCE_DIR}/../libs")
  elseif(WIN32)
    link_directories($ENV{PARDISO_LIB_PATH})
	message($ENV{PARDISO_LIB_PATH})
  endif()
endif()


add_subdirectory(../femlib       "${CMAKE_CURRENT_BINARY_DIR}/femlib")
add_subdirectory(../linsolver   "${CMAKE_CURRENT_BINARY_DIR}/linsolver")
add_subdirectory(../StepperCPU   "${CMAKE_CURRENT_BINARY_DIR}/StepperCPU")

FILE(GLOB EXE_SRC "*.cpp" "src/*.cpp")
ADD_DEFINITIONS(--std=c++0x)

INCLUDE_DIRECTORIES ("${PROJECT_SOURCE_DIR}/../femlib/include")
INCLUDE_DIRECTORIES ("${PROJECT_SOURCE_DIR}/../StepperCPU/include")
INCLUDE_DIRECTORIES ("${PROJECT_SOURCE_DIR}/../linsolver/include/pcgsolver")
INCLUDE_DIRECTORIES ("${PROJECT_SOURCE_DIR}/../linsolver/include")

INCLUDE_DIRECTORIES(include)

ADD_EXECUTABLE(Opt ${EXE_SRC})

TARGET_LINK_LIBRARIES(Opt StepperCPU)
TARGET_LINK_LIBRARIES(Opt FEMLib)
TARGET_LINK_LIBRARIES(Opt linsolver)

if (PARDISO_AVAILABLE)
  if(APPLE)
    TARGET_LINK_LIBRARIES(Opt mkl_intel mkl_intel_thread mkl_core)
  elseif(UNIX)
    TARGET_LINK_LIBRARIES(Opt pardiso openblas gfortran)
  elseif(WIN32)
	target_link_libraries(Opt libpardiso500-WIN-X86-64)
  endif()
endif()

if (WIN32)
  add_custom_command(TARGET Opt POST_BUILD
					   COMMAND set "var=${PROJECT_SOURCE_DIR}"
					   COMMAND set "var=%var:/=\\%"
					   COMMAND set "var1=${PROJECT_BINARY_DIR}"
					   COMMAND set "var1=%var1:/=\\%"
					   COMMAND copy "%var%\\..\\dll\\*.dll" "%var1%\\$<CONFIGURATION>" 
					   COMMENT "copying dlls ..."
				    )
endif()
if (WIN32)
TARGET_LINK_LIBRARIES(Opt IpOpt-vc10)
endif ()

