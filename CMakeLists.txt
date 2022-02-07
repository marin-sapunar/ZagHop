cmake_minimum_required (VERSION 2.8)
project (NAD_program)
enable_language (Fortran)

###################################################################################################
# Options.                                                                                        #
  if(NOT DEFINED CMAKE_BUILD_TYPE OR "${CMAKE_BUILD_TYPE}" STREQUAL "")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
        "Choose the type of build, options are: Release, Debug and DebugAll ..."
        FORCE)
  endif()
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Release Debug DebugAll)

  set(LINALG MKL CACHE STRING "Linear algebra library.")
  set_property(CACHE LINALG PROPERTY STRINGS MKL OpenBLAS)

  option(OMP "Use OMP." ON)
  option(STATIC "Create a static executable." OFF)

  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
  set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)

  include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})
  link_directories(${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

###################################################################################################
# Compiler specific options:                                                                      #
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    set(dialect "-fpp -assume realloc_lhs")
    set(CMAKE_Fortran_FLAGS_RELEASE  "-O3")
    set(CMAKE_Fortran_FLAGS_DEBUG    "-O2 -traceback -g")
    set(CMAKE_Fortran_FLAGS_DEBUGALL "-O0 -traceback -g -check all -debug all")
  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(dialect "-cpp -ffree-line-length-0")
    set(CMAKE_Fortran_FLAGS_RELEASE  "-O3")
    set(CMAKE_Fortran_FLAGS_DEBUG "-O2 -fbacktrace -g")
    set(CMAKE_Fortran_FLAGS_DEBUGALL "-O0 -fbacktrace -Wall -pedantic -fcheck=all -fbounds-check -g")
  endif()
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")

###################################################################################################
# Libraries.                                                                                      #
  add_definitions(-DLINALG=1)
  # CMake FindBLAS/FindLAPACK modules don't work for ilp64 version of MKL or for MKL F95 interface.
  # Set everything manually.
  if(LINALG STREQUAL "MKL")
    set(MKLROOT "" CACHE PATH "MKL root directory.")
    set(F95ROOT "" CACHE PATH "MKL F95 interface root directory.")
    if(MKLROOT STREQUAL "")
      set(MKLROOT $ENV{MKLROOT})
      if(NOT MKLROOT)
        message(FATAL_ERROR
                "When using MKL, set MKLROOT environment variable "
                "or call cmake with -DMKLROOT=/path/to/mklroot_dir")
      endif()
    endif()
    if(F95ROOT STREQUAL "")
      set(F95ROOT $ENV{F95ROOT})
      if(NOT F95ROOT)
        message(FATAL_ERROR
                "When using MKL, set F95ROOT environment variable "
                "or call cmake with -DF95ROOT=/path/to/f95root_dir")
      endif()
    endif()
    message(STATUS "MKLROOT = ${MKLROOT}")
    message(STATUS "F95ROOT = ${F95ROOT}")
    # F95 Interface:
    include_directories(${F95ROOT}/include/intel64/lp64)
    set(F95_LIBS ${F95ROOT}/lib/intel64/libmkl_blas95_lp64.a)
    set(F95_LIBS ${F95_LIBS} ${F95ROOT}/lib/intel64/libmkl_lapack95_lp64.a)
    set(LINKLIBS ${LINKLIBS} ${F95_LIBS})
    # MKL libraries.
    set(MKL_LIBS "-Wl,--no-as-needed -Wl,--start-group")
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
        set(MKL_LIBS ${MKL_LIBS} "-lmkl_intel_lp64 -lmkl_intel_thread")
    elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
        set(MKL_LIBS ${MKL_LIBS} "-lmkl_gf_lp64 -lmkl_gnu_thread")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -m64")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DMKL_LP64 -m64")
    elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
        set(MKL_LIBS ${MKL_LIBS} "-lmkl_intel_lp64 -lmkl_pgi_thread")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -pgf90libs")
    endif()
    set(MKL_LIBS ${MKL_LIBS} "-lmkl_core -Wl,--end-group -lpthread -lm -ldl")
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
        set(MKL_LIBS ${MKL_LIBS} "-liomp5")
    elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
        set(MKL_LIBS ${MKL_LIBS} "-lgomp")
    elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
        set(MKL_LIBS ${MKL_LIBS} "-pgf90libs")
    endif()
    set(LINKLIBS ${LINKLIBS} ${MKL_LIBS})
    add_definitions(-DMKL=1)
    add_definitions(-DLINALG_F95=1)
  elseif(LINALG STREQUAL "OpenBLAS")
    find_package(OpenBLAS)
    include_directories(${OpenBLAS_INCLUDE_DIRS})
    message(STATUS BLAS found: ${OpenBLAS_LIBRARIES})
    set(LINKLIBS ${LINKLIBS} ${OpenBLAS_LIBRARIES})
    add_definitions(-DOPENBLAS=1)
  endif()

  if(OMP)
    find_package(OpenMP REQUIRED)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
    message(STATUS "OpenMP libraries:  ${OpenMP_Fortran_LIB_NAMES}")
  endif()

  if(STATIC)
    set(LINKLIBS ${LINKLIBS} "-static")
  endif()

################################################################################
# Build                                                                        #
 #add_subdirectory(src/chemutils)
 #include_directories(${CMAKE_BINARY_DIR}/mod_chemutil)
 #file(GLOB SOURCES "src/core/*f90" "src/model/Ferretti.f90")
 #file(GLOB SOURCES "src/core/*f90" "src/model/Phenol.f90" "src/model/util.f90")
 #file(GLOB SOURCES "src/core/*f90" "src/model/Pyrazine.f90")
  file(GLOB SOURCES "src/core/*f90" "src/model/Ferretti.f90" "src/util/*f90")


  add_executable(nad.exe ${SOURCES})
  target_link_libraries(nad.exe ${LINKLIBS})
 #target_link_libraries(nad.exe chemutils ${LINKLIBS})
