cmake_minimum_required(VERSION 3.0)
project(mazebuilder)

#===========================================================================================
# ASSIMP CONFIGURATION

set(ENABLE_ASSIMP   ON CACHE BOOL "Add Open Asset Import Library (assimp) to the project" FORCE)
set(ASSIMP_BUILD_ASSIMP_TOOLS   OFF CACHE BOOL "Build ASSIMP tools" FORCE)
set(ASSIMP_BUILD_TESTS          OFF CACHE BOOL "Build tests" FORCE)
set(ASSIMP_INSTALL              OFF CACHE BOOL "Configure an install for ASSIMP" FORCE)

#===========================================================================================
# GLAD CONFIGURATION

set(GLAD_OUT_DIR        "${CMAKE_CURRENT_BINARY_DIR}" CACHE STRING "Output directory" FORCE)
set(GLAD_PROFILE        "core" CACHE STRING "OpenGL profile" FORCE)
set(GLAD_API            "gl=3.3" CACHE STRING "API type/version pairs, like \"gl=3.2,gles=\", no version means latest" FORCE)
set(GLAD_GENERATOR      "c" CACHE STRING "Language to generate the binding for" FORCE)
set(GLAD_EXTENSIONS     "" CACHE STRING "Path to extensions file or comma separated list of extensions, if missing all extensions are included" FORCE)
set(GLAD_SPEC           "gl" CACHE STRING "Name of the spec" FORCE)
set(GLAD_ALL_EXTENSIONS OFF CACHE BOOL "Include all extensions instead of those specified by GLAD_EXTENSIONS" FORCE)
set(GLAD_NO_LOADER      OFF CACHE BOOL "No loader" FORCE)
set(GLAD_REPRODUCIBLE   OFF CACHE BOOL "Reproducible build" FORCE)

#===========================================================================================

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Set directory paths
set(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(LIBS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/libs)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

set(sources ${SOURCE_DIR}/main.cpp)

add_executable(mazebuilder ${sources})

function(update_submodule module)

    find_package(Git REQUIRED)

    if(NOT EXISTS ${LIBS_DIR}/${module})
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive -- ${LIBS_DIR}/${module}
                WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                COMMAND_ERROR_IS_FATAL ANY)
    endif()

    if(EXISTS ${LIBS_DIR}/${module}/CMakeLists.txt)
        add_subdirectory(${LIBS_DIR}/${module})
    endif()

endfunction(update_submodule)

update_submodule(assimp)
update_submodule(glad)
update_submodule(glfw)
update_submodule(glm)
update_submodule(stb)

target_include_directories(mazebuilder PUBLIC ${LIBS_DIR}/assimp/include)
target_link_libraries(mazebuilder PUBLIC assimp)

target_include_directories(mazebuilder PUBLIC ${LIBS_DIR}/glm/include)
target_link_libraries(mazebuilder PUBLIC glm)

target_include_directories(mazebuilder PUBLIC ${LIBS_DIR}/glfw/include)
target_link_libraries(mazebuilder PUBLIC glfw)

target_include_directories(mazebuilder PUBLIC ${LIBS_DIR}/glad/include)
target_link_libraries(mazebuilder PUBLIC glad)

target_include_directories(mazebuilder PUBLIC ${LIBS_DIR}/assimp/include)
target_link_libraries(mazebuilder PUBLIC assimp)

target_include_directories(mazebuilder PUBLIC ${LIBS_DIR}/stb/include)

target_include_directories(mazebuilder PUBLIC ${SOURCE_DIR}/headers)

configure_file(${SOURCE_DIR}/headers/conf/root_directory.h.in ${SOURCE_DIR}/headers/root_directory.h)

find_package(OpenGL REQUIRED)
if (OpenGL_FOUND)
    target_include_directories(mazebuilder PRIVATE ${OPENGL_INCLUDE_DIRS})
    target_link_libraries(mazebuilder PRIVATE ${OPENGL_LIBRARIES})
else()
    message(FATAL_ERROR "OpenGL was not found on the system")
endif()

set_target_properties(mazebuilder PROPERTIES
        CXX_STANDARD 20
        CXX_STANDARD_REQUIRED YES
        CXX_EXTENSIONS NO)

# Set project folders
set_target_properties(mazebuilder PROPERTIES FOLDER ${PROJECT_NAME})
