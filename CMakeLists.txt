cmake_minimum_required(VERSION 3.8)

if (UNIX)
    set(CMAKE_CXX_COMPILER g++)
    set(CMAKE_CXX_FLAGS "-pthread")
endif ()

project(Twolayer)


set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "-O3 -Wall -pedantic")
set(CMAKE_CXX_STANDARD 14)
include_directories(src bin)
if(UNIX)
    set(CMAKE_CXX_FLAGS "-pthread")
endif()


include_directories(processData)
include_directories(/home/wwc129/c_tools/my_include)
link_directories(/home/wwc129/c_tools/my_lib)

set(source main.cpp src/curvature.cpp src/property.cpp src/IC.cpp src/sync.h
        src/sync.cpp model.cpp processData/init_cur.cpp 
        processData/write.h processData/init_IC.cpp processData/init_IC.h processData/init_total.cpp processData/init_total.h
        processData/init_sync.cpp  src/wwcFit.cpp
        )
    

add_executable(main ${source})


target_link_libraries(main wwcReadFile)
target_link_libraries(main wwcString)
target_link_libraries(main Print)
#target_link_libraries(main wwcFit)
target_link_libraries(main ${Boost_LIBRARIES})
