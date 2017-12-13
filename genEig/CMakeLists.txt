PROJECT (eig)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

set(PARDISO_AVAILABLE "true" CACHE BOOL "Chosen By User")

if (${PARDISO_AVAILABLE} MATCHES "true")
  add_definitions(-DPARDISO_AVAILABLE)
  message("Has pardiso")
  if(APPLE)
    include_directories(/opt/intel/mkl/include/)
    LINK_DIRECTORIES(/opt/intel/mkl/lib/)
  elseif(UNIX)
    LINK_DIRECTORIES("${PROJECT_SOURCE_DIR}/../libs")
  endif()
endif()
LINK_DIRECTORIES("${PROJECT_SOURCE_DIR}/../libs")
FILE(GLOB EXE_SRC "src/*.cpp")
ADD_DEFINITIONS(--std=c++0x)

INCLUDE_DIRECTORIES ("${PROJECT_SOURCE_DIR}/include")

INCLUDE_DIRECTORIES(include)

if(APPLE)
  link_directories("${PROJECT_SOURCE_DIR}/../lib/osx")
  link_directories("/opt/local/lib")
endif()

ADD_EXECUTABLE(eig ${EXE_SRC})

if (${PARDISO_AVAILABLE} MATCHES true)
  if(APPLE)
    TARGET_LINK_LIBRARIES(eig mkl_intel mkl_intel_thread mkl_core)
  elseif(UNIX)
    TARGET_LINK_LIBRARIES(eig arpack pardiso openblas gfortran gomp)
  else()
    TARGET_LINK_LIBRARIES(eig arpack pardiso openblas gfortran gomp)
  endif()
endif()

#TARGET_LINK_LIBRARIES(eig pthread )
