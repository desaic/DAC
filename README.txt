Design Dynamics

The code is a library mainly for FEM simulation with statics, dynamics and frictional contact.
To demonstrate the usage of the library, an example executable project in runLinux is provided.
A simple viewer that works together with the example executable is also provided in runViewer.
The code is mostly used in Windows.

Compile on OS x (Not well tested)

Dependencies:
  Eigen3
    Assumed to be in /opt/local/include/eigen3
    Since Eigen is headers only, it's easy to change the include directory to other places.
  mkl
    Assumed to be in default install location
    may not work with older version because function names changed from upper to lower cases and
    number of arguments reduced by 1.
  glfw3 for rendering
    Expects a libglfw3.a under "libs" in project root directory.
    
Compile the simulation
>cd runLinux
>mkdir debug
>cd debug

make debug or release build
>cmake ../ -DCMAKE_BUILD_TYPE=Debug

Make a directory for the results by default it's
>mkdir ../../models/simpleTest/0
Run the sim
>./Run ../../models/simpleTest/conf3d.txt
The result directory is specified by "resultdir" in conf3d.txt and confreplay3d.txt. 

Compile the viewer
>cd ../../runViewer
>mkdir release
>cd release
>cmake ../ -DCMAKE_BUILD_TYPE=Release
>./Run ../../models/simpleTest/confreplay3d.txt

Compiling in Windows with Visual Studio 2015
The binary dependencies can be downloaded on Github in Releases.
Download libs.zip, extract the content and place under root project directory.
run CMake and open the solution files.
Compile the main example program by pointing Cmake to the CMakeLists.txt file in runLinux.
Compile the viewer in runViewer.

Press ']' to start simulation
'p' to pause
'[' to single step
'o' to reset
'c' to toggle between contacts when paused.
'w' 's' 'a' 'd' 'r' 'f' to move forward backward left right up down
Click mouse to capture mouse. Move mouse to rotate camera after capture.
Click again to release mouse.

