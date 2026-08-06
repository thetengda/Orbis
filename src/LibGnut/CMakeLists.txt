#Minimum requirement of CMake version 
cmake_minimum_required(VERSION 3.15)

#Project name and version number
project(${LibGnut} VERSION 0.1.0)

# add define
# ==========================================================================================================================================
add_definitions(-D GREAT_LibGnut_LIBRARY)
# ==========================================================================================================================================

# link third party (vcpkg-managed, see src/vcpkg.json)
find_package(Eigen3 REQUIRED NO_MODULE)
find_package(spdlog REQUIRED)
find_package(pugixml REQUIRED)

# set include path
include_directories(${LibGnutSrc})

message(STATUS " include path for ${LibGnut} is : ${LibGnutSrc}")
#message(STATUS "                                                                     ")

# for all the header files and source_files
file(GLOB_RECURSE header_files ${PROJECT_SOURCE_DIR} *.h *.hpp)
file(GLOB_RECURSE source_files ${PROJECT_SOURCE_DIR} *.cpp)

# grouping
if(CMAKE_SYSTEM_NAME MATCHES "Windows")
    # add files to project
    source_group("CMake Files"  FILES CMakeLists.txt)
    source_group("Header Files" FILES header_files)
    source_group("Soruce Files" FILES source_files)

    # find files and put them into different groups
    foreach(fold ${file_fold})  
    file(GLOB_RECURSE tmp_header_files 
        ${PROJECT_SOURCE_DIR}/${fold}/*.h 
        ${PROJECT_SOURCE_DIR}/${fold}/*.hpp)
    file(GLOB_RECURSE tmp_source_files 
        ${PROJECT_SOURCE_DIR}/${fold}/*.cpp)
    source_group("Header Files\\${fold}"        FILES ${tmp_header_files})
    source_group("Source Files\\${fold}"        FILES ${tmp_source_files})
endforeach(fold)      
else()
    #message(STATUS "For Linux / unix, there is no need to group so detailed")
endif()

# grouping in VS
# ==========================================================================================================================================
if(${PrintFile})
    message(STATUS " ====================================================================================================================== ")
    message(STATUS " header_files_for_${PROJECT_NAME}" )
    foreach(list ${header_files})   
        message(STATUS " ==> ${list}")
    endforeach(list)   

    message(STATUS " source_files_for_${PROJECT_NAME}" )
    foreach(list ${source_files})   
        message(STATUS " ==> ${list}")
    endforeach(list)   
    message(STATUS " ====================================================================================================================== ")
endif()

# link project
# ==========================================================================================================================================

add_library(${PROJECT_NAME} SHARED ${header_files} ${source_files})
# ==========================================================================================================================================

# link third
# ==========================================================================================================================================
if(CMAKE_SYSTEM_NAME MATCHES "Windows")
    link_directories(${BUILD_DIR}/Lib/Debug)
    link_directories(${BUILD_DIR}/Lib/Release)
    link_directories(${BUILD_DIR}/Lib/RelWithDebInfo)
    link_directories(${BUILD_DIR}/Lib/MinSizeRel)
	link_directories(${BUILD_DIR}/Lib)
else()
    link_directories(${BUILD_DIR}/Lib)
endif()
# ==========================================================================================================================================

# link vcpkg-provided third-party targets (PUBLIC so consumers inherit the include paths)
target_link_libraries(${PROJECT_NAME} PUBLIC Eigen3::Eigen spdlog::spdlog pugixml::pugixml)






    

