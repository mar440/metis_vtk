
cmake_minimum_required(VERSION 2.8)

PROJECT(testVTK CXX)







SET (METIS_LIBRARIES 
    "/data_space/france/metis-5.1.0/build/Linux-x86_64/libmetis/libmetis.a" 
    CACHE PATH "mydoc")
SET (METIS_INCLUDEDIRS 
    "/data_space/france/metis-5.1.0/include/" 
    CACHE PATH "mydoc")


include_directories(${METIS_INCLUDEDIRS} ${CMAKE_CURRENT_SOURCE_DIR})


find_package(VTK 7.1.0 REQUIRED NO_MODULE)
include(${VTK_USE_FILE})
add_executable(metisVTK metisVTK.cxx)


target_link_libraries(metisVTK ${VTK_LIBRARIES} ${METIS_LIBRARIES})
