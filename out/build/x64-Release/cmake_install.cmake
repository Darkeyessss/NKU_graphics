# Install script for directory: C:/Users/18176/Desktop/fluid-sim/code

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "C:/Users/18176/Desktop/fluid-sim/code/out/install/x64-Release")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "RelWithDebInfo")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("C:/Users/18176/Desktop/fluid-sim/code/out/build/x64-Release/third_party/imgui/cmake_install.cmake")
  include("C:/Users/18176/Desktop/fluid-sim/code/out/build/x64-Release/third_party/glad/cmake_install.cmake")
  include("C:/Users/18176/Desktop/fluid-sim/code/out/build/x64-Release/third_party/glm/cmake_install.cmake")
  include("C:/Users/18176/Desktop/fluid-sim/code/out/build/x64-Release/common/cmake_install.cmake")
  include("C:/Users/18176/Desktop/fluid-sim/code/out/build/x64-Release/fluid2d/Lagrangian/cmake_install.cmake")
  include("C:/Users/18176/Desktop/fluid-sim/code/out/build/x64-Release/fluid2d/Eulerian/cmake_install.cmake")
  include("C:/Users/18176/Desktop/fluid-sim/code/out/build/x64-Release/fluid3d/Lagrangian/cmake_install.cmake")
  include("C:/Users/18176/Desktop/fluid-sim/code/out/build/x64-Release/fluid3d/Eulerian/cmake_install.cmake")
  include("C:/Users/18176/Desktop/fluid-sim/code/out/build/x64-Release/ui/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "C:/Users/18176/Desktop/fluid-sim/code/out/build/x64-Release/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
