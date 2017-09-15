PROJECT (RunViewer)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

#find_package(Eigen3)
#include_directories (${EIGEN3_INCLUDE_DIR})
set(FREEIMAGE_AVAILABLE ON)
if(FREEIMAGE_AVAILABLE)
  add_definitions(-DFREEIMAGE_AVAILABLE)
endif()
if(APPLE)
  include_directories(/opt/local/include/)
  include_directories(/opt/local/include/eigen3)
  include_directories("${PROJECT_SOURCE_DIR}/../include/osx")
  LINK_DIRECTORIES("${PROJECT_SOURCE_DIR}/../libs")
elseif(UNIX)
  include_directories("${PROJECT_SOURCE_DIR}/../include/unix")
  include_directories("${PROJECT_SOURCE_DIR}/../include/")
elseif(WIN32)
	include_directories("${PROJECT_SOURCE_DIR}/../include/")
	include_directories("${PROJECT_SOURCE_DIR}/../include/glfw")
	include_directories("${PROJECT_SOURCE_DIR}/../include/FreeImage")
	LINK_DIRECTORIES("$(ProjectDir)/../../libs")
endif()

add_subdirectory(../femlib       "${CMAKE_CURRENT_BINARY_DIR}/femlib")
add_subdirectory(../viewer       "${CMAKE_CURRENT_BINARY_DIR}/viewer")

FILE(GLOB EXE_SRC "*.cpp")
ADD_DEFINITIONS(--std=c++0x)

INCLUDE_DIRECTORIES ("${PROJECT_SOURCE_DIR}/../femlib/include")
INCLUDE_DIRECTORIES ("${PROJECT_SOURCE_DIR}/../viewer/include")
INCLUDE_DIRECTORIES ("${PROJECT_SOURCE_DIR}/../MCCD/include")

INCLUDE_DIRECTORIES(include)


ADD_EXECUTABLE(Run ${EXE_SRC})

TARGET_LINK_LIBRARIES(Run FEMLib)
TARGET_LINK_LIBRARIES(Run FEMViewer)

if(APPLE)
  FIND_LIBRARY(COCOA_LIBRARY Cocoa)
  FIND_LIBRARY(OpenGL_LIBRARY OpenGL)
  FIND_LIBRARY(IOKIT_LIBRARY iokit)
  FIND_LIBRARY(COREVIDEO_LIBRARY corevideo)
  target_link_libraries(Run glfw3 ${COCOA_LIBRARY} ${OpenGL_LIBRARY} ${IOKIT_LIBRARY} ${COREVIDEO_LIBRARY})
elseif(UNIX)
  TARGET_LINK_LIBRARIES(Run X11  Xxf86vm Xrandr Xi Xinerama Xcursor gomp)
  TARGET_LINK_LIBRARIES(Run png glfw3 GL GLU pthread)
elseif(WIN32)
  TARGET_LINK_LIBRARIES(Run FreeImage glfw3 opengl32 glu32)
  add_custom_command(TARGET Run POST_BUILD
					   COMMAND set "var=${PROJECT_SOURCE_DIR}"
					   COMMAND set "var=%var:/=\\%"
					   COMMAND set "var1=${PROJECT_BINARY_DIR}"
					   COMMAND set "var1=%var1:/=\\%"
					   COMMAND copy "%var%\\..\\dll\\*.dll" "%var1%\\$<CONFIGURATION>" 
					   COMMENT "copying dlls ..."
				   )
endif()
