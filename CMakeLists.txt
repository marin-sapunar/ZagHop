cmake_minimum_required (VERSION 2.8.12)
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
  option(QUANTICS "Compile with Quantics interface." OFF)

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
    message(STATUS "LINALG = ${LINALG}")
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
    include_directories(${MKLROOT}/include/)
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
  elseif(LINALG STREQUAL "Quantics")
    message(STATUS Quantics BLAS: $ENV{QUANTICS_DIR})
    add_definitions(-DQUANTICS=1)

  endif()

  if(QUANTICS)
    add_definitions(-DQUANTICS=1)
    set(QUANTICS_DIR $ENV{QUANTICS_DIR})
    set(QUANTICS_PLATFORM $ENV{QUANTICS_PLATFORM})
    set(QUANTICS_COMPILER $ENV{QUANTICS_COMPILER})
    if(NOT (QUANTICS_DIR AND QUANTICS_PLATFORM AND QUANTICS_COMPILER))
      message(FATAL_ERROR
              "To compile with Quantics interface, set QUANTICS_DIR, QUANTICS_PLATFORM "
              "and QUANTICS_COMPILER environment variables.")
    endif()
    if(QUANTICS_COMPILER STREQUAL "default")
      set(QUANTICS_COMPILER "gfortran")
    endif()
    set(QUANTICS_OBJ ${QUANTICS_DIR}/object/${QUANTICS_PLATFORM}/${QUANTICS_COMPILER})
    set(QUANTICS_LIB ${QUANTICS_DIR}/bin/dyn_libs)
    include_directories(${QUANTICS_OBJ}/include)
    set(QUANTICS_INTER "interface/quantics_inter.f90")
    set(LINKLIBS ${LINKLIBS} ${QUANTICS_OBJ}/versions.o)
    set(LINKLIBS ${LINKLIBS} ${QUANTICS_OBJ}/mctdh.a)
    set(LINKLIBS ${LINKLIBS} ${QUANTICS_OBJ}/propwf.a)
    set(LINKLIBS ${LINKLIBS} ${QUANTICS_OBJ}/geninwf.a)
    set(LINKLIBS ${LINKLIBS} ${QUANTICS_OBJ}/genoper.a)
    set(LINKLIBS ${LINKLIBS} ${QUANTICS_OBJ}/gendvr.a)
    set(LINKLIBS ${LINKLIBS} ${QUANTICS_OBJ}/quanticslib.a)
    set(LINKLIBS ${LINKLIBS} ${QUANTICS_OBJ}/opfuncs.a)
    set(LINKLIBS ${LINKLIBS} ${QUANTICS_OBJ}/libode.a)
    set(LINKLIBS ${LINKLIBS} ${QUANTICS_OBJ}/quanticsmod.a)
    set(LINKLIBS ${LINKLIBS} ${QUANTICS_OBJ}/includes.a)
    set(LINKLIBS ${LINKLIBS} ${QUANTICS_OBJ}/globinc.a)
    set(LINKLIBS ${LINKLIBS} ${QUANTICS_OBJ}/libnum.a)
    set(LINKLIBS ${LINKLIBS} ${QUANTICS_OBJ}/libsys.a)
    set(LINKLIBS ${LINKLIBS} ${QUANTICS_OBJ}/libarpack.a)
    set(LINKLIBS ${LINKLIBS} ${QUANTICS_OBJ}/liblapack.a)
    set(LINKLIBS ${LINKLIBS} ${QUANTICS_OBJ}/libblas.a)
    set(LINKLIBS ${LINKLIBS} ${QUANTICS_OBJ}/libomp.a)
    set(LINKLIBS ${LINKLIBS} "-J ${QUANTICS_OBJ}/include -lsqlite3 -L ${QUANTICS_DIR}/bin/dyn_libs -lsrf -lusrf -lsqlite3")
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
  file(GLOB SOURCES "src/core/*f90" "src/util/*f90")
  if(QUANTICS)
      set(SOURCES ${SOURCES} ${QUANTICS_INTER})
  endif()


  add_executable(nad.exe ${SOURCES})
  target_link_libraries(nad.exe ${LINKLIBS})

  add_custom_command(
        TARGET nad.exe POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_SOURCE_DIR}/interface/*py
                ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/)
