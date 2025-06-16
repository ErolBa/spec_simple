cmake_minimum_required(VERSION 3.20.0)

project(spec LANGUAGES C Fortran)

set(CMAKE_VERBOSE_MAKEFILE ON CACHE BOOL "ON")
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Find OpenMP
find_package(OpenMP)

# Find Threads
set (THREADS_PREFER_PTHREAD_FLAG ON) # Prefer pthreads
find_package(Threads)

# Find libm
set(POW_LIBS "")
include(CheckLibraryExists)
check_library_exists(m pow "" LIBM)
if(LIBM)
    list(APPEND POW_LIBS "m")
endif()


# Find & Configure MPI
# If using the compiler wrapper, there is no need to find the MPI libraries.
get_filename_component (compiler ${CMAKE_Fortran_COMPILER} NAME)
message(STATUS "compiler is ${compiler}")
if (${compiler} STREQUAL mpiifort OR ${compiler} STREQUAL mpifort)
else ()
    find_package (MPI REQUIRED)
    if (MPI_Fortran_FOUND)
        set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MPI_Fortran_COMPILE_FLAGS}")
        set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_Fortran_LINK_FLAGS}")
        include_directories (${MPI_Fortran_INCLUDE_PATH})

        if (NOT DEFINED ${MPI_OVERSUBSCRIBE_FLAG}) 
            execute_process(COMMAND ${MPIEXEC} --version OUTPUT_VARIABLE MPI_VERSION)
            if (${MPI_VERSION} MATCHES "Open MPI" OR
                ${MPI_VERSION} MATCHES "OpenRTE"  OR
                ${MPI_VERSION} MATCHES "slurm"    OR
                ${MPI_VERSION} MATCHES "aprun")
                set (MPI_OVERSUBSCRIBE_FLAG "--oversubscribe")
            endif ()
        endif ()
    endif ()
endif ()


if(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DIFORT")
endif()

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER "10.0")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch")
    endif()
endif()

find_program(AWK awk mawk gawk)

find_package (LAPACK REQUIRED)

configure_file(downloadFindFFTW.cmake.in findFFTW-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/findFFTW-download )
if(result)
    message(FATAL_ERROR "CMake step for findFFTW failed: ${result}")
else()
    message("CMake step for findFFTW completed (${result}).")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/findFFTW-download )
if(result)
    message(FATAL_ERROR "Build step for findFFTW failed: ${result}")
endif()

set(findFFTW_DIR ${CMAKE_CURRENT_BINARY_DIR}/findFFTW-src)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${findFFTW_DIR}")

find_package(FFTW REQUIRED)

find_package(HDF5 COMPONENTS C Fortran HL REQUIRED)
message(STATUS "HDF5 libraries found: ${HDF5_FOUND}")
message(STATUS "HDF5 C include directors : ${HDF5_C_INCLUDE_DIRS}")
message(STATUS "HDF5 Fortran include directors : ${HDF5_Fortran_INCLUDE_DIRS}")

set (EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/build/bin)
set (LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/build/lib)

set (CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/build/modules)
add_subdirectory(src)
if(SKBUILD)
    add_subdirectory("Utilities/python_wrapper")
endif()

add_compile_definitions(OPENMP)
