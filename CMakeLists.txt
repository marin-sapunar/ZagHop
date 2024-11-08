cmake_minimum_required (VERSION 3.18)
project(ZagHop LANGUAGES Fortran)

################################################################################
# Options.                                                                     #

  set(LINALG LAPACK CACHE STRING "Linear algebra library.")
  set_property(CACHE LINALG PROPERTY STRINGS LAPACK MKL QUANTICS)

  option(OMP "Use OMP." ON)
  option(STATIC "Create a static executable." OFF)
  option(QUANTICS "Compile with Quantics interface." OFF)
  option(BLAS_64 "Use 64bit integer interface for BLAS/LAPACK" OFF)

  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
  set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)

  include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})
  link_directories(${CMAKE_LIBRARY_OUTPUT_DIRECTORY})


################################################################################
# Targets:                                                                     #
  file(GLOB CORE_SOURCE "src/core/*f90")
  file(GLOB UTIL_SOURCE "src/util/*f90")

  add_library(zagutil STATIC ${UTIL_SOURCE})
  if(QUANTICS)
      set(CORE_SOURCE ${CORE_SOURCE} "interface/quantics_inter.f90")
  endif()

  add_executable(zaghop ${CORE_SOURCE})
  target_link_libraries(zaghop PRIVATE zagutil)

################################################################################
# Compiler specific options:                                                   #
  if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(dialect "-fpp")
  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(dialect "-cpp -ffree-line-length-0")
    set_source_files_properties(random_pcg.f90 PROPERTIES
                                COMPILE_FLAGS -fwrapv)
  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
    set(dialect "-e Z -M 878,7206")
  endif()
  set_source_files_properties(./src/util/random_pcg.f90 PROPERTIES
                              COMPILE_FLAGS -O0)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")

################################################################################
# Libraries.                                                                   #

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
    ${QUANTICS_OBJ}/libomp.a)
    target_link_libraries(zaghop PRIVATE srf usrf sqlite3)
  endif()


  add_definitions(-DLINALG=1)
  if(LINALG STREQUAL "LAPACK")
    if (BLAS_64)
      set(BLA_SIZEOF_INTEGER 8)
      add_definitions(-DBLA_INT64=1)
    else()
      set(BLA_SIZEOF_INTEGER 4)
    endif()
    find_package(LAPACK REQUIRED)
    set(LINALG_LIBS LAPACK::LAPACK)
  elseif(LINALG STREQUAL "MKL")
    if (BLAS_64)
      set(MKL_INTERFACE "ilp64")
      add_definitions(-DBLA_INT64=1)
    else()
      set(MKL_INTERFACE "lp64")
    endif()
    find_package(MKL CONFIG REQUIRED)
    set(LINALG_LIBS MKL::MKL)
    add_definitions(-DMKL=1)
  elseif(LINALG STREQUAL "Quantics")
    if(NOT QUANTICS)
        message(FATAL_ERROR
                "QUANTICS=ON needs to be set in order to use LINALG=Quantics.")
    endif()
    if(BLAS_64)
        message(FATAL_ERROR
                "BLAS_64=ON is not compatible with LINALG=QUANTICS.")
    endif()
    message(STATUS Quantics BLAS: $ENV{QUANTICS_DIR})
    set(LINALG_LIBS
          ${QUANTICS_OBJ}/libarpack.a
          ${QUANTICS_OBJ}/liblapack.a
          ${QUANTICS_OBJ}/libblas.a
          )
  endif()
  target_link_libraries(zaghop PRIVATE ${LINALG_LIBS})
  target_link_libraries(zagutil ${LINALG_LIBS})

  if(STATIC)
    target_link_libraries(zaghop PRIVATE "-static")
  endif()

################################################################################
# Tests.                                                                   #
  enable_testing()

  add_executable(util_test tests/unit/util_test.f90)
  target_link_libraries(util_test zagutil ${LINALG_LIBS})
  add_test(NAME util_unit_test COMMAND util_test)

  add_executable(linalg_test tests/unit/linalg_test.f90)
  target_link_libraries(linalg_test zagutil ${LINALG_LIBS})
  add_test(NAME linalg_unit_test COMMAND linalg_test)

  find_package(Python)
  if(Python_FOUND)
    message(STATUS "Python executable: ${Python_EXECUTABLE}")
    add_custom_command(TARGET zaghop POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_directory
              ${CMAKE_SOURCE_DIR}/tests/integration ${CMAKE_BINARY_DIR}/tests)

    add_test(NAME python-zaghop-tests
      COMMAND ${Python_EXECUTABLE} -m unittest -v test_zaghop
      WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/tests")

    set_tests_properties(python-zaghop-tests PROPERTIES ENVIRONMENT_MODIFICATION
      PATH=path_list_append:${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
  else()
     message(STATUS "Python not found.")
  endif()
