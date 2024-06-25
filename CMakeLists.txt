cmake_minimum_required(VERSION 3.15)
project(sub_pixel_edge_detect LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
#if (WIN32)
#    MESSAGE(STATUS "Now is windows")
#    set(OpenCV_DIR "../libs/opencv/build/x64/vc16/lib")
#    set(spdlog_DIR "../libs/spdlog/build")
#endif ()
if(MSVC)
    add_compile_options(/FS)
    #    set(CMAKE_C_COMPILER "C:/Program Files/Microsoft Visual Studio/2022/Community/VC/Tools/MSVC/14.39.33519/bin/Hostx86/x86/cl.exe")
    message(STATUS "MSVC")
    set(OpenCV_DIR "C:/tools/opencv/build/x64/vc16/lib")
endif()

#
#find_package(Eigen3 REQUIRED)
#find_package(xlnt REQUIRED)
#find_package(fmt REQUIRED)
#find_package(spdlog REQUIRED)
find_package(OpenCV REQUIRED)

include_directories( ${OpenCV_INCLUDE_DIRS}  )
MESSAGE(STATUS "OpenCV_INCLUDE_DIRS: ${OpenCV_INCLUDE_DIRS}")

add_executable(sub_pixel_edge_detect main.cpp WriteFile.h edgeTest.cpp edgeTest.h)
target_link_libraries(sub_pixel_edge_detect ${OpenCV_LIBS})
#if(NOT TARGET spdlog)
#    # Stand-alone build
#
#endif()

#add_subdirectory(edge-profile)
#add_subdirectory(test)