# DesignDynamics
Design Dynamics

Compile on OS x

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

Press ']' to start simulation
'p' to pause
'[' to single step
'o' to reset
'c' to toggle between contacts when paused.
'w' 's' 'a' 'd' 'r' 'f' to move forward backward left right up down
Click mouse to capture mouse. Move mouse to rotate camera after capture.
Click again to release mouse.

