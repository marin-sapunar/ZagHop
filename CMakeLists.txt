cmake_minimum_required (VERSION 3.5)
project(ZagHop LANGUAGES Fortran)

################################################################################
# Options.                                                                     #

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


################################################################################
# Targets:                                                                     #
  file(GLOB SOURCES "src/core/*f90" "src/util/*f90")
  if(QUANTICS)
      set(SOURCES ${SOURCES} "interface/quantics_inter.f90")
  endif()

  add_executable(zaghop ${SOURCES})

################################################################################
# Compiler specific options:                                                   #
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    set(dialect "-fpp -assume realloc_lhs")
  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(dialect "-cpp -ffree-line-length-0")
  endif()
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")

################################################################################
# Libraries.                                                                   #

  add_definitions(-DLINALG=1)
  if(LINALG STREQUAL "MKL")
    set(MKL_INTERFACE "lp64")
    set(ENABLE_BLAS95 ON CACHE BOOL "Enables BLAS Fortran95 API")
    set(ENABLE_LAPACK95 ON CACHE BOOL "Enables LAPACK Fortran95 API")
    find_package(MKL CONFIG REQUIRED)
    target_compile_options(zaghop PRIVATE
        $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
    target_include_directories(zaghop PRIVATE
        $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
    target_link_libraries(zaghop PRIVATE $<LINK_ONLY:MKL::MKL>)
    add_definitions(-DMKL=1)
    add_definitions(-DLINALG_F95=1)
    # MKLConfig.cmake can only set up BLAS95 and LAPACK95 with Intel compilers
    # so they are specified manually here for gfortran.
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      set(F95ROOT "" CACHE PATH "MKL F95 interface root directory.")
      if(F95ROOT STREQUAL "")
        set(F95ROOT $ENV{F95ROOT})
        if(NOT F95ROOT)
          message(FATAL_ERROR
                  "When using MKL, set F95ROOT environment variable "
                  "or call cmake with -DF95ROOT=/path/to/f95root_dir")
        endif()
      endif()
      message(STATUS "F95ROOT = ${F95ROOT}")
      target_include_directories(zaghop PRIVATE
          ${F95ROOT}/include/intel64/lp64)
      target_link_libraries(zaghop PRIVATE
          ${F95ROOT}/lib/intel64/libmkl_blas95_lp64.a)
      target_link_libraries(zaghop PRIVATE
          ${F95ROOT}/lib/intel64/libmkl_lapack95_lp64.a)
    endif()
  elseif(LINALG STREQUAL "Quantics")
    if(NOT QUANTICS)
        message(FATAL_ERROR
                "QUANTICS=ON needs to be set in order to use LINALG=Quantics.")
    endif()
    message(STATUS Quantics BLAS: $ENV{QUANTICS_DIR})
    add_definitions(-DLINALG_QUANTICS=1)
  endif()

  if(QUANTICS)
    add_definitions(-DQUANTICS=1)
    set(QUANTICS_DIR $ENV{QUANTICS_DIR})
    set(QUANTICS_PLATFORM $ENV{QUANTICS_PLATFORM})
    set(QUANTICS_COMPILER $ENV{QUANTICS_COMPILER})
    if(NOT (QUANTICS_DIR AND QUANTICS_PLATFORM AND QUANTICS_COMPILER))
      message(FATAL_ERROR
              "To compile with Quantics interface, set QUANTICS_DIR, "
              "QUANTICS_PLATFORM and QUANTICS_COMPILER environment variables.")
    endif()
    if(QUANTICS_COMPILER STREQUAL "default")
      set(QUANTICS_COMPILER "gfortran")
    endif()
    set(QUANTICS_OBJ 
       ${QUANTICS_DIR}/object/${QUANTICS_PLATFORM}/${QUANTICS_COMPILER})
    target_include_directories(zaghop PRIVATE ${QUANTICS_OBJ}/include)
    target_link_directories(zaghop PRIVATE ${QUANTICS_DIR}/bin/dyn_libs)
    target_link_directories(zaghop PRIVATE ${QUANTICS_OBJ})
    target_link_libraries(zaghop PRIVATE
    ${QUANTICS_OBJ}/versions.o
    ${QUANTICS_OBJ}/mctdh.a
    ${QUANTICS_OBJ}/propwf.a
    ${QUANTICS_OBJ}/geninwf.a
    ${QUANTICS_OBJ}/genoper.a
    ${QUANTICS_OBJ}/gendvr.a
    ${QUANTICS_OBJ}/quanticslib.a
    ${QUANTICS_OBJ}/opfuncs.a
    ${QUANTICS_OBJ}/libode.a
    ${QUANTICS_OBJ}/quanticsmod.a
    ${QUANTICS_OBJ}/includes.a
    ${QUANTICS_OBJ}/globinc.a
    ${QUANTICS_OBJ}/libnum.a
    ${QUANTICS_OBJ}/libsys.a
    ${QUANTICS_OBJ}/libarpack.a
    ${QUANTICS_OBJ}/liblapack.a
    ${QUANTICS_OBJ}/libblas.a
    ${QUANTICS_OBJ}/libomp.a)
    target_link_libraries(zaghop PRIVATE srf usrf sqlite3)
  endif()

  if(STATIC)
    target_link_libraries(zaghop PRIVATE "-static")
  endif()
