# Minimum required CMake version
cmake_minimum_required(VERSION 3.10)

# Project name and version
project(MARLEY)

# Include GNUInstallDirs to get standard install directories
include(GNUInstallDirs)

# Set default C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Option to ignore ROOT
option(IGNORE_ROOT "Ignore any ROOT installations that may be present." OFF)
option(IGNORE_HEPMC3 "Ignore any HepMC3 installations that may be present." OFF)
option(BUILD_TESTS "Build MARLEY tests" OFF)
option(BUILD_DEBUG "Build with debugging symbols" OFF)

# Set build type
if(BUILD_DEBUG)
    set(CMAKE_BUILD_TYPE Debug)
    add_compile_options(-O0 -g)
else()
    if(NOT CMAKE_BUILD_TYPE)
        set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    endif()
    add_compile_options(-O3)
endif()

# Read MARLEY version from .VERSION file if it exists
if(EXISTS "${CMAKE_SOURCE_DIR}/.VERSION")
    file(READ "${CMAKE_SOURCE_DIR}/.VERSION" MARLEY_VERSION)
    string(STRIP "${MARLEY_VERSION}" MARLEY_VERSION)
else()
    # Attempt to get the Git revision
    find_package(Git QUIET)
    if(GIT_FOUND AND EXISTS "${CMAKE_SOURCE_DIR}/.git")
        execute_process(
            COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
            OUTPUT_VARIABLE GIT_COMMIT_HASH
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
        )
        if(GIT_COMMIT_HASH)
            # Check if the working tree is dirty
            execute_process(
                COMMAND ${GIT_EXECUTABLE} diff --quiet
                WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
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

# Define preprocessor macros
add_compile_definitions(MARLEY_VERSION="${MARLEY_VERSION}")
add_compile_definitions(MARLEY_GIT_REVISION="${GIT_COMMIT_HASH}")

# Set include directories
include_directories(${CMAKE_SOURCE_DIR}/include)

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

execute_process(COMMAND ${GSL_CONFIG_EXECUTABLE} --libs
    OUTPUT_VARIABLE GSL_LIBS
    OUTPUT_STRIP_TRAILING_WHITESPACE)

# Separate the flags into lists
separate_arguments(GSL_CFLAGS)
separate_arguments(GSL_LIBS)

# Include GSL headers
include_directories(${GSL_CFLAGS})

# Find or use built-in HepMC3
if(NOT IGNORE_HEPMC3)
    find_program(HEPMC3_CONFIG_EXECUTABLE NAMES HepMC3-config PATHS ENV PATH)
    if(HEPMC3_CONFIG_EXECUTABLE)
        execute_process(COMMAND ${HEPMC3_CONFIG_EXECUTABLE} --includedir
            OUTPUT_VARIABLE HEPMC3_INCLUDE_DIR
            OUTPUT_STRIP_TRAILING_WHITESPACE)
        execute_process(COMMAND ${HEPMC3_CONFIG_EXECUTABLE} --cflags
            OUTPUT_VARIABLE HEPMC3_CXXFLAGS
            OUTPUT_STRIP_TRAILING_WHITESPACE)
        execute_process(COMMAND ${HEPMC3_CONFIG_EXECUTABLE} --libs
            OUTPUT_VARIABLE HEPMC3_LIBS
            OUTPUT_STRIP_TRAILING_WHITESPACE)
        separate_arguments(HEPMC3_CXXFLAGS)
        separate_arguments(HEPMC3_LIBS)
        message(STATUS "Found HepMC3: ${HEPMC3_INCLUDE_DIR}")
    else()
        message(STATUS "HepMC3 not found, using built-in HepMC3.")
        set(USE_BUILTIN_HEPMC3 ON)
        set(HEPMC3_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/hepmc3/include")
        include_directories(${HEPMC3_INCLUDE_DIR})
        file(GLOB HEPMC3_SOURCES CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/hepmc3/*.cc")
        add_library(HepMC3 OBJECT ${HEPMC3_SOURCES})
        set_target_properties(HepMC3 PROPERTIES POSITION_INDEPENDENT_CODE ON)
        target_include_directories(HepMC3 PRIVATE ${HEPMC3_INCLUDE_DIR})
    endif()
else()
    message(STATUS "Ignoring HepMC3 as per user request.")
endif()

# Optionally find ROOT
if(NOT IGNORE_ROOT)
    find_package(ROOT COMPONENTS Core RIO Hist Tree REQUIRED)
    if(ROOT_FOUND)
        set(USE_ROOT ON)
        include(${ROOT_USE_FILE})
        message(STATUS "Found ROOT ${ROOT_VERSION}: ${ROOT_EXECUTABLE}")

        # Construct the path to root-config
        if(NOT DEFINED ROOT_CONFIG_EXECUTABLE)
            if(DEFINED ROOT_BINDIR)
                set(ROOT_CONFIG_EXECUTABLE "${ROOT_BINDIR}/root-config")
            else()
                message(FATAL_ERROR "Cannot locate root-config.")
            endif()
        endif()

        # Verify that root-config exists
        if(NOT EXISTS "${ROOT_CONFIG_EXECUTABLE}")
            message(FATAL_ERROR "Cannot find root-config at ${ROOT_CONFIG_EXECUTABLE}")
        endif()

        # Get ROOT's cflags using root-config
        execute_process(
            COMMAND ${ROOT_CONFIG_EXECUTABLE} --cflags
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
            set(ROOT_CXX_STANDARD ${ROOT_CXX_STANDARD_VERSION})

            ## Set the C++ standard to match ROOT
            #set(CMAKE_CXX_STANDARD ${ROOT_CXX_STANDARD_VERSION})
        else()
            message(WARNING "Could not detect the C++ standard used by ROOT. Defaulting to C++17.")
            set(ROOT_CXX_STANDARD 17)
        endif()

        # Remove any '-std=' flags from ROOT_CXXFLAGS to avoid conflicts
        string(REGEX REPLACE "-std=c\\+\\+[0-9a-z]+" "" ROOT_CXXFLAGS "${ROOT_CXXFLAGS}")
        # Also remove any duplicate spaces
        string(REGEX REPLACE "[ ]+" " " ROOT_CXXFLAGS "${ROOT_CXXFLAGS}")

        # Include ROOT headers
        include_directories(${ROOT_INCLUDE_DIRS})

        # Define USE_ROOT macro
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
file(GLOB MARLEY_SOURCES CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/src/*.cc")

# Exclude specific source files
list(REMOVE_ITEM MARLEY_SOURCES
    "${CMAKE_SOURCE_DIR}/src/marley.cc"
    "${CMAKE_SOURCE_DIR}/src/marley_root.cc"
    "${CMAKE_SOURCE_DIR}/src/marsum.cc"
    "${CMAKE_SOURCE_DIR}/src/OutputFileRoot.cc"
    "${CMAKE_SOURCE_DIR}/src/OutputFilePlainRoot.cc"
    "${CMAKE_SOURCE_DIR}/src/MacroEventFileReader.cc"
    "${CMAKE_SOURCE_DIR}/src/mroot.cc"
    "${CMAKE_SOURCE_DIR}/src/JSONConfig.cc"
    "${CMAKE_SOURCE_DIR}/src/OutputFile.cc"
    "${CMAKE_SOURCE_DIR}/src/OutputFileAscii.cc"
    "${CMAKE_SOURCE_DIR}/src/EventFileReader.cc"
)

# Build the MARLEY shared library
add_library(MARLEY SHARED ${MARLEY_SOURCES})
target_include_directories(MARLEY PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_compile_options(MARLEY PRIVATE ${GSL_CFLAGS})
target_link_libraries(MARLEY PRIVATE ${GSL_LIBS})
set_target_properties(MARLEY PROPERTIES VERSION ${MARLEY_VERSION})

if(HEPMC3_CONFIG_EXECUTABLE)
    target_compile_options(MARLEY PUBLIC ${HEPMC3_CXXFLAGS})
    target_include_directories(MARLEY PUBLIC ${HEPMC3_INCLUDE_DIR})
    target_link_libraries(MARLEY PUBLIC ${HEPMC3_LIBS})
else()
    target_link_libraries(MARLEY PUBLIC HepMC3)
endif()

# Build MARLEY_ROOT shared library if ROOT is found
if(USE_ROOT)
    set(ROOT_SOURCES
        "${CMAKE_SOURCE_DIR}/src/marley_root.cc"
        "${CMAKE_SOURCE_DIR}/src/JSONConfig.cc"
        "${CMAKE_SOURCE_DIR}/src/OutputFile.cc"
        "${CMAKE_SOURCE_DIR}/src/OutputFileAscii.cc"
        "${CMAKE_SOURCE_DIR}/src/OutputFileRoot.cc"
        "${CMAKE_SOURCE_DIR}/src/OutputFilePlainRoot.cc"
        "${CMAKE_SOURCE_DIR}/src/EventFileReader.cc"
        #"${CMAKE_SOURCE_DIR}/src/MacroEventFileReader.cc"
    )

    # Generate ROOT dictionary
    set(DICTIONARY_HEADERS
        #"${CMAKE_SOURCE_DIR}/include/marley/Event.hh"
        #"${CMAKE_SOURCE_DIR}/include/marley/Particle.hh"
        "${CMAKE_SOURCE_DIR}/include/marley/Parity.hh"
        #"${CMAKE_SOURCE_DIR}/include/marley/MacroEventFileReader.hh"
        # Add other headers as needed
    )
    set(ROOT_LINKDEF "${CMAKE_SOURCE_DIR}/include/marley/marley_linkdef.hh")
    ROOT_GENERATE_DICTIONARY(G__MARLEY_ROOT ${DICTIONARY_HEADERS} LINKDEF ${ROOT_LINKDEF})
    list(APPEND ROOT_SOURCES "${CMAKE_CURRENT_BINARY_DIR}/G__MARLEY_ROOT.cxx")

    add_library(MARLEY_ROOT SHARED ${ROOT_SOURCES})
    target_include_directories(MARLEY_ROOT PUBLIC ${CMAKE_SOURCE_DIR}/include)
    target_include_directories(MARLEY_ROOT PRIVATE ${ROOT_INCLUDE_DIRS})
    target_compile_options(MARLEY_ROOT PRIVATE ${GSL_CFLAGS})
    target_link_libraries(MARLEY_ROOT PUBLIC MARLEY ${ROOT_LIBRARIES} ${GSL_LIBS})
    set_target_properties(MARLEY_ROOT PROPERTIES VERSION ${MARLEY_VERSION})

    # Set the C++ standard for MARLEY_ROOT
    set_target_properties(MARLEY_ROOT PROPERTIES CXX_STANDARD ${ROOT_CXX_STANDARD})
    set_target_properties(MARLEY_ROOT PROPERTIES CXX_STANDARD_REQUIRED YES)
endif()

# Build the 'marley' executable
add_executable(marley "${CMAKE_SOURCE_DIR}/src/marley.cc")
target_compile_options(marley PRIVATE ${GSL_CFLAGS})
target_include_directories(marley PRIVATE ${GSL_CFLAGS})
target_link_libraries(marley PRIVATE MARLEY ${GSL_LIBS})
if(USE_ROOT)
    target_link_libraries(marley PRIVATE MARLEY_ROOT ${ROOT_LIBRARIES})
    set_target_properties(marley PROPERTIES CXX_STANDARD ${ROOT_CXX_STANDARD})
    set_target_properties(marley PROPERTIES CXX_STANDARD_REQUIRED YES)
endif()

# Build 'marsum' executable if ROOT is found
if(USE_ROOT)
    add_executable(marsum "${CMAKE_SOURCE_DIR}/src/marsum.cc")
    target_compile_options(marsum PRIVATE ${GSL_CFLAGS})
    target_include_directories(marsum PRIVATE ${GSL_CFLAGS} ${ROOT_INCLUDE_DIRS})
    target_link_libraries(marsum PRIVATE MARLEY MARLEY_ROOT ${ROOT_LIBRARIES} ${GSL_LIBS})
    set_target_properties(marsum PROPERTIES CXX_STANDARD ${ROOT_CXX_STANDARD})
    set_target_properties(marsum PROPERTIES CXX_STANDARD_REQUIRED YES)
endif()

# Build 'mroot' script if ROOT is found
if(USE_ROOT)
    install(PROGRAMS ${CMAKE_SOURCE_DIR}/src/scripts/mroot DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()

# Build test executable if BUILD_TESTS is enabled
if(BUILD_TESTS)
    add_executable(martest)
    file(GLOB TEST_SOURCES CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/src/tests/*.cc")
    target_sources(martest PRIVATE ${TEST_SOURCES})
    target_compile_options(martest PRIVATE ${GSL_CFLAGS})
    target_link_libraries(martest PRIVATE MARLEY ${GSL_LIBS})
    if(USE_ROOT)
        target_link_libraries(martest PRIVATE MARLEY_ROOT ${ROOT_LIBRARIES})
        set_target_properties(martest PROPERTIES CXX_STANDARD ${ROOT_CXX_STANDARD})
        set_target_properties(martest PROPERTIES CXX_STANDARD_REQUIRED YES)
    endif()
endif()

# Install rules
install(TARGETS marley DESTINATION ${CMAKE_INSTALL_BINDIR})
install(TARGETS MARLEY LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/marley DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

if(USE_ROOT)
    install(TARGETS marsum DESTINATION ${CMAKE_INSTALL_BINDIR})
    install(TARGETS MARLEY_ROOT LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

# Install data files
install(DIRECTORY ${CMAKE_SOURCE_DIR}/data/ DESTINATION ${CMAKE_INSTALL_DATADIR}/marley)

# Configure and install marley-config script
configure_file(${CMAKE_SOURCE_DIR}/src/scripts/marley-config.in ${CMAKE_CURRENT_BINARY_DIR}/marley-config @ONLY)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/marley-config DESTINATION ${CMAKE_INSTALL_BINDIR})

# Set RPATH for installed binaries
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_LIBDIR}")

# Post-installation commands (if needed)
if(NOT WIN32)
    install(CODE "execute_process(COMMAND ldconfig)")
endif()

