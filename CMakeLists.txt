SET (VTK_DIR "/data_space/Program_files/VTK-build")

cmake_minimum_required(VERSION 2.8)

PROJECT(testVTK CXX)







SET (METIS_LIBRARIES 
    "/data_space/france/metis-5.1.0/build/Linux-x86_64/libmetis/libmetis.a" 
    CACHE PATH "mydoc")
SET (METIS_INCLUDEDIRS 
    "/data_space/france/metis-5.1.0/include/" 
    CACHE PATH "mydoc")

SET (VTK_DIR "/data_space/Program_files/VTK-build")
SET (CMAKE_BUILD_TYPE "Debug")

include_directories(${METIS_INCLUDEDIRS} ${CMAKE_CURRENT_SOURCE_DIR})


find_package(VTK 7.1.0 REQUIRED NO_MODULE)
include(${VTK_USE_FILE})
add_executable(metisVTK metisVTK.cxx)


target_link_libraries(metisVTK ${VTK_LIBRARIES} ${METIS_LIBRARIES})
