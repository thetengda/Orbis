# CONFIGURE_DEPENDS source discovery requires CMake 3.12 or newer.
cmake_minimum_required(VERSION 3.12)

#Project name and version number
project(${LibGREAT} VERSION 0.1.0)

# add define
# ==========================================================================================================================================
add_definitions(-D GREAT_LibGREAT_LIBRARY)
add_definitions(-D GLOG_NO_ABBREVIATED_SEVERITIES)

# ==========================================================================================================================================

# link third party (vcpkg-managed, see src/vcpkg.json)
find_package(Eigen3 REQUIRED NO_MODULE)
find_package(glfw3 REQUIRED)
find_package(Ceres REQUIRED)
find_package(pugixml REQUIRED)
find_package(spdlog REQUIRED)

# set include path
include_directories(${LibGnutSrc})
include_directories(${LibGREATSrc})

message(STATUS " include path for ${LibGREAT} is : ${LibGREATSrc}")

# for all the header files and source files.  CONFIGURE_DEPENDS makes a
# subsequent build re-run CMake when a new RAW factor/manager source is
# added to the tree.
file(GLOB_RECURSE header_files CONFIGURE_DEPENDS
    "${PROJECT_SOURCE_DIR}/*.h" "${PROJECT_SOURCE_DIR}/*.hpp")
file(GLOB_RECURSE source_files CONFIGURE_DEPENDS
    "${PROJECT_SOURCE_DIR}/*.cpp")

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
# vcpkg-provided imported targets (ceres carries eigen3 + glog; glfw links opengl32)
target_link_libraries(${PROJECT_NAME}    PUBLIC Eigen3::Eigen ceres glfw pugixml::pugixml
                                        PRIVATE spdlog::spdlog)
if(WIN32)
    # gviewer uses legacy OpenGL (glBegin/glClear/...) resolved from opengl32.lib
    target_link_libraries(${PROJECT_NAME} PUBLIC opengl32)
endif()
target_link_libraries(${PROJECT_NAME}   PUBLIC ${LibGnut})
add_dependencies(${PROJECT_NAME}       ${LibGnut})

set(pvt     GREAT_PVT)
add_subdirectory(${ROOT}/app/${pvt}    ${BUILD_DIR}/${pvt})
SET_PROPERTY(TARGET ${pvt}     PROPERTY FOLDER "app")

set(msf       GREAT_MSF)
add_subdirectory(${ROOT}/app/${msf}   ${BUILD_DIR}/${msf})
SET_PROPERTY(TARGET ${msf}    PROPERTY FOLDER "app")

set(pvtfgo       GREAT_PVTFGO)
add_subdirectory(${ROOT}/app/${pvtfgo}   ${BUILD_DIR}/${pvtfgo})
SET_PROPERTY(TARGET ${pvtfgo}    PROPERTY FOLDER "app")


set(gins       GREAT_GINSFGO)
add_subdirectory(${ROOT}/app/${gins}   ${BUILD_DIR}/${gins})
SET_PROPERTY(TARGET ${gins}    PROPERTY FOLDER "app")
