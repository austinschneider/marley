cmake_minimum_required(VERSION 3.10)
project(marley)

# Include GNUInstallDirs to get standard install directories like CMAKE_INSTALL_BINDIR
include(GNUInstallDirs)

# Set the default C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Read the MARLEY version from the .VERSION file if it exists
if(EXISTS "${PROJECT_SOURCE_DIR}/.VERSION")
    file(READ "${PROJECT_SOURCE_DIR}/.VERSION" MARLEY_VERSION)
    string(STRIP "${MARLEY_VERSION}" MARLEY_VERSION)
else()
    # Attempt to get the GIT revision
    find_package(Git QUIET)
    if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
        execute_process(
            COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
            OUTPUT_VARIABLE GIT_COMMIT_HASH
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
        )
        if(GIT_COMMIT_HASH)
            # Check if the working tree is dirty
            execute_process(
                COMMAND ${GIT_EXECUTABLE} diff --quiet
                WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                RESULT_VARIABLE GIT_DIFF_EXIT_CODE
            )
            if(NOT GIT_DIFF_EXIT_CODE EQUAL 0)
                set(GIT_COMMIT_HASH "${GIT_COMMIT_HASH}-dirty")
            endif()
        else()
            set(GIT_COMMIT_HASH "unknown version")
        endif()
    else()
        set(GIT_COMMIT_HASH "unknown version")
    endif()
    set(MARLEY_VERSION "${GIT_COMMIT_HASH}")
endif()

# Define the MARLEY_VERSION macro
add_compile_definitions(MARLEY_VERSION="${MARLEY_VERSION}")

# Option to ignore ROOT
option(IGNORE_ROOT "Ignore any ROOT installations that may be present." OFF)

# Find GSL using gsl-config
find_program(GSL_CONFIG_EXECUTABLE NAMES gsl-config PATHS ENV PATH)
if(NOT GSL_CONFIG_EXECUTABLE)
    message(FATAL_ERROR "Could not find 'gsl-config'. Please ensure that GSL is installed and 'gsl-config' is in your PATH.")
endif()

execute_process(COMMAND ${GSL_CONFIG_EXECUTABLE} --version
    OUTPUT_VARIABLE GSL_VERSION_STRING
    OUTPUT_STRIP_TRAILING_WHITESPACE)
message(STATUS "Found GSL version ${GSL_VERSION_STRING}")

execute_process(COMMAND ${GSL_CONFIG_EXECUTABLE} --cflags
    OUTPUT_VARIABLE GSL_CFLAGS
    OUTPUT_STRIP_TRAILING_WHITESPACE)
message(STATUS "GSL CFLAGS: ${GSL_CFLAGS}")

execute_process(COMMAND ${GSL_CONFIG_EXECUTABLE} --libs
    OUTPUT_VARIABLE GSL_LIBS
    OUTPUT_STRIP_TRAILING_WHITESPACE)
message(STATUS "GSL LIBS: ${GSL_LIBS}")

# Separate the flags into lists
separate_arguments(GSL_CFLAGS)
separate_arguments(GSL_LIBS)

# Include GSL headers
include_directories(${GSL_CFLAGS})

# Optionally find ROOT
if(NOT IGNORE_ROOT)
    find_package(ROOT COMPONENTS RIO Hist Tree REQUIRED)
    if(ROOT_FOUND)
        set(USE_ROOT ON)
        include(${ROOT_USE_FILE})
        message(STATUS "Found ROOT ${ROOT_VERSION}: ${ROOT_EXECUTABLE}")

        # Construct the path to root-config
        if(DEFINED ROOT_BINDIR)
            set(ROOT_CONFIG "${ROOT_BINDIR}/root-config")
        else()
            message(FATAL_ERROR "ROOT_BINDIR is not defined. Cannot locate root-config.")
        endif()

        # Verify that root-config exists
        if(NOT EXISTS "${ROOT_CONFIG}")
            message(FATAL_ERROR "Cannot find root-config at ${ROOT_CONFIG}")
        endif()

        # Get ROOT's cflags using root-config
        execute_process(
            COMMAND ${ROOT_CONFIG} --cflags
            OUTPUT_VARIABLE ROOT_CXXFLAGS
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )

        # Extract the C++ standard used by ROOT
        string(REGEX MATCH "-std=c\\+\\+[0-9a-z]+" ROOT_STD_FLAG "${ROOT_CXXFLAGS}")
        if(ROOT_STD_FLAG)
            # Remove '-std=' to get the standard version
            string(REPLACE "-std=" "" ROOT_CXX_STANDARD_STRING "${ROOT_STD_FLAG}")
            # Extract the numeric part (e.g., c++17 -> 17)
            string(REGEX REPLACE "c\\+\\+([0-9]+)[a-z]*" "\\1" ROOT_CXX_STANDARD_VERSION "${ROOT_CXX_STANDARD_STRING}")
            message(STATUS "ROOT was built with C++ standard: c++${ROOT_CXX_STANDARD_VERSION}")

            # Set the ROOT_CXX_STANDARD variable
            set(ROOT_CXX_STANDARD ${ROOT_CXX_STANDARD_VERSION})
        else()
            message(WARNING "Could not detect the C++ standard used by ROOT. Defaulting to C++14.")
            set(ROOT_CXX_STANDARD 14)
        endif()

        # Remove any '-std=' flags from ROOT_CXXFLAGS to avoid conflicts
        string(REGEX REPLACE "-std=c\\+\\+[0-9a-z]+" "" ROOT_CXXFLAGS "${ROOT_CXXFLAGS}")
        # Also remove any duplicate spaces
        string(REGEX REPLACE "[ ]+" " " ROOT_CXXFLAGS "${ROOT_CXXFLAGS}")

        # No need to use ROOT_CXXFLAGS directly; ROOT's imported targets handle this
        add_compile_definitions(USE_ROOT)
    else()
        set(USE_ROOT OFF)
        message(WARNING "ROOT not found. Building without ROOT support.")
    endif()
else()
    set(USE_ROOT OFF)
    message(STATUS "Ignoring ROOT as per user request.")
endif()

# Set compiler warnings and flags
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall -Wextra -Wpedantic -Wno-error=unused-parameter -Wcast-align)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5.0)
        add_compile_options(-Wno-shadow)
    endif()
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        add_compile_options(-Wno-keyword-macro -Wno-missing-braces)
    endif()
endif()

# Source files
file(GLOB MARLEY_SOURCES CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/src/*.cc")

list(REMOVE_ITEM MARLEY_SOURCES
    "${PROJECT_SOURCE_DIR}/src/marley.cc"
    "${PROJECT_SOURCE_DIR}/src/marley_root.cc"
    "${PROJECT_SOURCE_DIR}/src/marsum.cc"
    "${PROJECT_SOURCE_DIR}/src/RootJSONConfig.cc"
    "${PROJECT_SOURCE_DIR}/src/RootOutputFile.cc"
    "${PROJECT_SOURCE_DIR}/src/RootEventFileReader.cc"
    "${PROJECT_SOURCE_DIR}/src/MacroEventFileReader.cc"
)

# Print the list of source files
message(STATUS "MARLEY_SOURCES contains the following files:")
foreach(SOURCE_FILE ${MARLEY_SOURCES})
    message(STATUS "${SOURCE_FILE}")
endforeach()

# Build the MARLEY shared library
message(STATUS "Adding MARLEY library target.")
add_library(MARLEY SHARED ${MARLEY_SOURCES})
set_target_properties(MARLEY PROPERTIES OUTPUT_NAME "MARLEY")
target_compile_options(MARLEY PRIVATE ${GSL_CFLAGS})
target_link_libraries(MARLEY PRIVATE ${GSL_LIBS})
set_target_properties(MARLEY PROPERTIES VERSION ${MARLEY_VERSION})
target_include_directories(MARLEY PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/marley>
)

# Print properties of the MARLEY target
get_target_property(MARLEY_INCLUDE_DIRS MARLEY INCLUDE_DIRECTORIES)
get_target_property(MARLEY_COMPILE_OPTIONS MARLEY COMPILE_OPTIONS)
get_target_property(MARLEY_LINK_LIBRARIES MARLEY LINK_LIBRARIES)
message(STATUS "MARLEY target include directories: ${MARLEY_INCLUDE_DIRS}")
message(STATUS "MARLEY target compile options: ${MARLEY_COMPILE_OPTIONS}")
message(STATUS "MARLEY target link libraries: ${MARLEY_LINK_LIBRARIES}")

# Build MARLEY_ROOT shared library if ROOT is found
if(USE_ROOT)
    set(ROOT_SOURCES
        "${PROJECT_SOURCE_DIR}/src/marley_root.cc"
        "${PROJECT_SOURCE_DIR}/src/RootJSONConfig.cc"
        "${PROJECT_SOURCE_DIR}/src/OutputFile.cc"
        "${PROJECT_SOURCE_DIR}/src/RootOutputFile.cc"
        "${PROJECT_SOURCE_DIR}/src/RootEventFileReader.cc"
        "${PROJECT_SOURCE_DIR}/src/MacroEventFileReader.cc"
    )

    # Generate ROOT dictionary
    set(DICTIONARY_HEADERS
        "${PROJECT_SOURCE_DIR}/include/marley/Particle.hh"
        "${PROJECT_SOURCE_DIR}/include/marley/Event.hh"
        "${PROJECT_SOURCE_DIR}/include/marley/MacroEventFileReader.hh"
        "${PROJECT_SOURCE_DIR}/include/marley/Parity.hh"
    )
    ROOT_GENERATE_DICTIONARY(G__MARLEY_ROOT ${DICTIONARY_HEADERS} LINKDEF "${PROJECT_SOURCE_DIR}/include/marley/marley_linkdef.hh")
    list(APPEND ROOT_SOURCES "${CMAKE_CURRENT_BINARY_DIR}/G__MARLEY_ROOT.cxx")

    add_library(MARLEY_ROOT SHARED ${ROOT_SOURCES})
    target_include_directories(MARLEY_ROOT PRIVATE ${GSL_CFLAGS})
    target_compile_options(MARLEY_ROOT PRIVATE ${GSL_CFLAGS})
    target_link_libraries(MARLEY_ROOT PUBLIC MARLEY ${ROOT_LIBRARIES} ${GSL_LIBS})
    set_target_properties(MARLEY_ROOT PROPERTIES VERSION ${MARLEY_VERSION})

    # Set the C++ standard for MARLEY_ROOT
    set_target_properties(MARLEY_ROOT PROPERTIES CXX_STANDARD ${ROOT_CXX_STANDARD})
    set_target_properties(MARLEY_ROOT PROPERTIES CXX_STANDARD_REQUIRED YES)
endif()

# Build the 'marley' executable
add_executable(marley_exec "${PROJECT_SOURCE_DIR}/src/marley.cc")
set_target_properties(marley_exec PROPERTIES OUTPUT_NAME "marley")
target_include_directories(marley_exec PRIVATE ${GSL_CFLAGS})
target_compile_options(marley_exec PRIVATE ${GSL_CFLAGS})
target_link_libraries(marley_exec PRIVATE MARLEY ${GSL_LIBS})
target_include_directories(marley_exec PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/marley>
)

if(USE_ROOT)
    target_link_libraries(marley_exec PRIVATE MARLEY_ROOT ${ROOT_LIBRARIES})
    # Set the C++ standard for 'marley'
    set_target_properties(marley_exec PROPERTIES CXX_STANDARD ${ROOT_CXX_STANDARD})
    set_target_properties(marley PROPERTIES CXX_STANDARD_REQUIRED YES)
endif()

# Build 'marsum' executable if ROOT is found
if(USE_ROOT)
    add_executable(marsum "${PROJECT_SOURCE_DIR}/src/marsum.cc")
    target_include_directories(marsum PRIVATE ${GSL_CFLAGS})
    target_compile_options(marsum PRIVATE ${GSL_CFLAGS})
    target_link_libraries(marsum PRIVATE MARLEY MARLEY_ROOT ${ROOT_LIBRARIES} ${GSL_LIBS})
    target_include_directories(marsum PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/marley>
    )
    # Set the C++ standard for 'marsum'
    set_target_properties(marsum PROPERTIES CXX_STANDARD ${ROOT_CXX_STANDARD})
    set_target_properties(marsum PROPERTIES CXX_STANDARD_REQUIRED YES)
endif()

# Build test executable
add_executable(martest)
file(GLOB TEST_SOURCES CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/src/tests/*.cc")
target_sources(martest PRIVATE ${TEST_SOURCES})
target_include_directories(martest PRIVATE ${GSL_CFLAGS})
target_compile_options(martest PRIVATE ${GSL_CFLAGS})
target_link_libraries(martest PRIVATE MARLEY ${GSL_LIBS})
target_include_directories(martest PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/marley>
)

if(USE_ROOT)
    target_link_libraries(martest PRIVATE MARLEY_ROOT ${ROOT_LIBRARIES})
    # Set the C++ standard for 'martest'
    set_target_properties(martest PROPERTIES CXX_STANDARD ${ROOT_CXX_STANDARD})
    set_target_properties(martest PROPERTIES CXX_STANDARD_REQUIRED YES)
endif()

# Install rules
install(TARGETS marley_exec
    EXPORT marleyConfig
    DESTINATION ${CMAKE_INSTALL_BINDIR}
)
install(TARGETS MARLEY
    EXPORT marleyConfig
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/marley
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

if(USE_ROOT)
    install(TARGETS marsum DESTINATION ${CMAKE_INSTALL_BINDIR})
    install(TARGETS MARLEY_ROOT LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

# Install data files
install(DIRECTORY ${PROJECT_SOURCE_DIR}/data/ DESTINATION ${CMAKE_INSTALL_DATADIR}/marley)

# Install examples
install(DIRECTORY ${PROJECT_SOURCE_DIR}/examples DESTINATION ${CMAKE_INSTALL_DATADIR}/marley)

# Install mroot script if ROOT is found
if(USE_ROOT)
    install(PROGRAMS ${PROJECT_SOURCE_DIR}/src/scripts/mroot DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()

# Configure and install marley-config script
configure_file(${PROJECT_SOURCE_DIR}/src/scripts/marley-config.in ${CMAKE_CURRENT_BINARY_DIR}/marley-config @ONLY)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/marley-config DESTINATION ${CMAKE_INSTALL_BINDIR})

# Set RPATH for installed binaries
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_LIBDIR}")

# Handle build types (Debug/Release)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# CMake configuration for MARLEY package
include(CMakePackageConfigHelpers)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/marleyConfigVersion.cmake"
    VERSION ${MARLEY_VERSION}
    COMPATIBILITY AnyNewerVersion
)
export(EXPORT marleyConfig
    FILE "${CMAKE_CURRENT_BINARY_DIR}/marleyConfig.cmake"
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/marleyConfigVersion.cmake"
    DESTINATION share/marley/cmake
)

install(EXPORT marleyConfig
    DESTINATION share/marley/cmake
)
