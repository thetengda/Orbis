#Minimum requirement of CMake version 
cmake_minimum_required(VERSION 3.15)

#Project name and version number
project(${msf})

file(GLOB header_files     *.h)
file(GLOB source_files    *.cpp)

source_group("CMake Files" FILES CMakeLists.txt)
source_group("Header Files" FILES header_files)
source_group("Soruce Files" FILES source_files)

set(include_path
    ${LibGnutSrc}
    ${LibGREATSrc})
    
include_directories(${include_path})

add_definitions(-D GLOG_NO_ABBREVIATED_SEVERITIES)

add_executable(${PROJECT_NAME} ${header_files} ${source_files})

if(CMAKE_SYSTEM_NAME MATCHES "Windows")
    set(link_path
        ${BUILD_DIR}/Lib/Debug
        ${BUILD_DIR}/Lib/Release
        ${BUILD_DIR}/Lib/RelWithDebInfo
        ${BUILD_DIR}/Lib/MinSizeRel
    )
    link_directories(${link_path})
else()
    set(link_path
        ${BUILD_DIR}/Lib)
    link_directories(${link_path})
endif()

set(lib_list 
    ${LibGnut}
    ${LibGREAT})
    
target_link_libraries(${PROJECT_NAME} ${lib_list})

# Stage vcpkg runtime DLLs (ceres/glog/fmt/glfw/pugixml/spdlog) next to this EXE.
add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_RUNTIME_DLLS:${PROJECT_NAME}>
            $<TARGET_FILE_DIR:${PROJECT_NAME}>
    COMMAND_EXPAND_LISTS)

add_dependencies(${PROJECT_NAME} ${lib_list})




