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
  option(COVERAGE "Enable code coverage." OFF)

  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
  set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)

  include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})
  link_directories(${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

################################################################################
# Dependencies.                                                                #
  set(JSONFORTRAN_SOURCE_DIR
      ${CMAKE_SOURCE_DIR}/external/json-fortran)
  set(JSONFORTRAN_ENABLE_TESTS OFF CACHE BOOL "" FORCE)
  if(NOT EXISTS ${JSONFORTRAN_SOURCE_DIR}/CMakeLists.txt)
    message(FATAL_ERROR
            "json-fortran submodule not found at ${JSONFORTRAN_SOURCE_DIR}. "
            "Run 'git submodule update --init --recursive'.")
  endif()
  add_subdirectory(${JSONFORTRAN_SOURCE_DIR} ${CMAKE_BINARY_DIR}/jsonfortran)

################################################################################
# Targets:                                                                     #
  file(GLOB CORE_SOURCE "src/core/*f90")
  list(FILTER CORE_SOURCE EXCLUDE REGEX "main\\.f90$")
  file(GLOB UTIL_SOURCE "src/util/*f90")
  file(GLOB POT_SOURCE "src/potential/*f90")
  list(FILTER POT_SOURCE EXCLUDE REGEX "eval_point\\.f90$")
  list(FILTER POT_SOURCE EXCLUDE REGEX "quantics\\.f90$")
  if(QUANTICS)
      set(CORE_SOURCE ${CORE_SOURCE} "src/potential/quantics.f90")
  endif()

  add_library(zagutil STATIC ${UTIL_SOURCE})
  target_link_libraries(zagutil PUBLIC jsonfortran)
  target_include_directories(zagutil PUBLIC
                             $<TARGET_PROPERTY:jsonfortran,BINARY_DIR>)
  add_library(zagpot STATIC ${POT_SOURCE})
  target_link_libraries(zagpot PRIVATE zagutil jsonfortran)
  target_include_directories(zagpot PRIVATE
                             $<TARGET_PROPERTY:jsonfortran,BINARY_DIR>)

  add_library(zagcore STATIC ${CORE_SOURCE})
  target_link_libraries(zagcore PUBLIC zagpot zagutil jsonfortran)
  target_include_directories(zagcore PRIVATE
                             $<TARGET_PROPERTY:jsonfortran,BINARY_DIR>)

  add_executable(zaghop "src/core/main.f90")
  target_link_libraries(zaghop PRIVATE zagcore)

  add_executable(eval_point "src/potential/eval_point.f90")
  target_link_libraries(eval_point PRIVATE zagpot)
  target_link_libraries(eval_point PRIVATE zagutil)
  target_link_libraries(eval_point PRIVATE ${LINALG_LIBS})

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
    if(QUANTICS_COMPILER STREQUAL "default")
      set(QUANTICS_COMPILER "gfortran")
    endif()
    if(NOT (QUANTICS_DIR AND QUANTICS_PLATFORM AND QUANTICS_COMPILER))
      message(FATAL_ERROR
              "To compile with Quantics interface, set QUANTICS_DIR, "
              "QUANTICS_PLATFORM and QUANTICS_COMPILER environment variables.")
    endif()
    set(QUANTICS_OBJ
       ${QUANTICS_DIR}/object/${QUANTICS_PLATFORM}/${QUANTICS_COMPILER})
    target_include_directories(zagcore PRIVATE ${QUANTICS_OBJ}/include)
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
  target_link_libraries(zagcore PUBLIC ${LINALG_LIBS})
  target_link_libraries(zagpot PRIVATE ${LINALG_LIBS})
  target_link_libraries(zagutil PUBLIC ${LINALG_LIBS})

  if(STATIC)
    target_link_libraries(zaghop PRIVATE "-static")
  endif()

  set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
################################################################################
# Tests.                                                                   #
  enable_testing()
  include(CTest)

  if(COVERAGE)
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      set(COVERAGE_FLAGS "--coverage -fprofile-arcs -ftest-coverage")
      target_compile_options(zaghop PRIVATE --coverage -fprofile-arcs -ftest-coverage)
      target_compile_options(zagcore PRIVATE --coverage -fprofile-arcs -ftest-coverage)
      target_link_options(zaghop PRIVATE --coverage)
      target_link_options(zagcore INTERFACE --coverage)
    endif()
  endif()

  add_executable(util_test tests/unit/util_test.f90)
  target_link_libraries(util_test zagutil ${LINALG_LIBS})
  add_test(NAME util_unit_test COMMAND util_test)

  add_executable(linalg_test tests/unit/linalg_test.f90)
  target_link_libraries(linalg_test zagutil ${LINALG_LIBS})
  add_test(NAME linalg_unit_test COMMAND linalg_test)

  # Register a unit test executable as a single ctest test (tests/unit/run_cases.cmake). The
  # executable runs once without arguments and passes with exit status 0, then once for each
  # CASE <argument> <regex>, which passes if its output matches the regular expression (cases
  # which stop with an error or print a warning). The runs are sequential, in the working
  # directory unit_tests/<name> of the build directory. A regex can't contain ';'.
  function(zaghop_add_unit_test name exe)
    set(_dir ${CMAKE_BINARY_DIR}/unit_tests/${name})
    file(MAKE_DIRECTORY ${_dir})
    set(_defs "")
    set(_n 0)
    set(_args ${ARGN})
    while(_args)
      unset(_regex)
      list(POP_FRONT _args _key _arg _regex)
      if(NOT _key STREQUAL "CASE" OR NOT DEFINED _regex)
        message(FATAL_ERROR "zaghop_add_unit_test(${name}): expected CASE <argument> <regex>.")
      endif()
      math(EXPR _n "${_n} + 1")
      list(APPEND _defs "-DCASE_${_n}_ARG=${_arg}" "-DCASE_${_n}_REGEX=${_regex}")
    endwhile()
    add_test(NAME ${name}
      COMMAND ${CMAKE_COMMAND} -DEXE=$<TARGET_FILE:${exe}> -DN_CASE=${_n} ${_defs}
              -P ${CMAKE_SOURCE_DIR}/tests/unit/run_cases.cmake
      WORKING_DIRECTORY ${_dir})
  endfunction()

  add_executable(pot_test tests/unit/pot_test.f90)
  target_link_libraries(pot_test zagpot zagutil jsonfortran ${LINALG_LIBS})
  target_include_directories(pot_test PRIVATE
                             $<TARGET_PROPERTY:jsonfortran,BINARY_DIR>)
  zaghop_add_unit_test(pot_unit_test pot_test
    # Cases which should stop with an error.
    CASE doublet "only implemented for singlets and triplets"
    CASE effective "need a potential with spin_states = components"
    CASE no_coupling "No singlet-triplet spin-orbit coupling at the reference"
    CASE unsigned_projection "The projection model is always signed"
    CASE model_parameter "Unrecognized parameter of model granucci-st: sigma"
    CASE subset_too_many "3 states of multiplicity 1 requested, 2 available"
    CASE subset_multiplicity "There are no states of multiplicity 2"
    # External program (external_json evaluator).
    CASE external_no_response "qm_out.json not found after running: true"
    CASE external_exit "External program failed: false"
    CASE external_point_id "qm_out.json: response for point [0-9]+, expected point [0-9]+"
    CASE external_missing_gradient "Missing gradient of a requested state"
    CASE external_not_computed "Gradient of state 2 was not computed at this point"
    # Transformations between the bases.
    CASE transform_spaces "No transformation from the mch_sf to the adiabatic basis"
    CASE transform_no_diabatic "The diabatic basis is not available for this potential"
    CASE transform_wrong_point "do not link the requested point"
    CASE transform_no_link "No overlaps with a previous point are stored")

  add_executable(sh_test tests/unit/sh_test.f90)
  target_link_libraries(sh_test zagcore)
  zaghop_add_unit_test(sh_unit_test sh_test
    # Cases which should stop with an error.
    CASE ld_wrong_link "do not link the requested point"
    CASE ld_same_point "was already processed"
    CASE initial_no_weight "has no weight in the states of the adiabatic basis"
    # A warning, not an error.
    CASE subset_crossing "Smallest singular value of the overlaps")

  # Register each test_* method of a Python unittest class as a separate ctest test.
  function(zaghop_add_python_tests module class prefix)
    cmake_parse_arguments(ARG "" "" "ENVIRONMENT" ${ARGN})
    set(_file "${CMAKE_SOURCE_DIR}/tests/integration/${module}.py")
    # Reconfigure when the file changes, so new test methods get registered.
    set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${_file}")
    set(_env PATH=path_list_append:${CMAKE_RUNTIME_OUTPUT_DIRECTORY} ${ARG_ENVIRONMENT})
    file(STRINGS "${_file}" _lines REGEX "^    def test_[A-Za-z0-9_]+")
    foreach(_line ${_lines})
      string(REGEX MATCH "test_[A-Za-z0-9_]+" _test "${_line}")
      add_test(NAME "ZagHop:${prefix}::${_test}"
        COMMAND ${Python_EXECUTABLE} -m unittest -v "${module}.${class}.${_test}"
        WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/tests")
      set_tests_properties("ZagHop:${prefix}::${_test}" PROPERTIES
        ENVIRONMENT_MODIFICATION "${_env}")
    endforeach()
  endfunction()

  find_package(Python)
  if(Python_FOUND)
    message(STATUS "Python executable: ${Python_EXECUTABLE}")
    add_custom_target(zagpy ALL
      COMMAND ${Python_EXECUTABLE} -m pip install -e ${CMAKE_SOURCE_DIR}
              --quiet --no-build-isolation
      COMMENT "Installing zagpy Python package (editable)"
    )
    # A target without outputs is always out of date, so the tests are copied on every build,
    # also when only the test files changed.
    add_custom_target(integration_tests ALL
      COMMAND ${CMAKE_COMMAND} -E copy_directory
              ${CMAKE_SOURCE_DIR}/tests/integration ${CMAKE_BINARY_DIR}/tests
      COMMENT "Copying integration tests")

    zaghop_add_python_tests(test_tully       TullyTest      Tully)
    zaghop_add_python_tests(test_granucci    GranucciTest   Granucci)
    zaghop_add_python_tests(test_vc_pyrazine VcPyrazineTest VC:Pyrazine)
    zaghop_add_python_tests(test_vc_so2      VcSo2Test      VC:SO2)
    zaghop_add_python_tests(test_zagpy       ZagpyTest      zagpy)
    if(QUANTICS)
      # The tests run Quantics to generate the restart, dvr and oper files.
      set(_quantics_env
          QUANTICS_DIR=set:${QUANTICS_DIR}
          QUANTICS_PLATFORM=set:${QUANTICS_PLATFORM}
          LD_LIBRARY_PATH=path_list_append:${QUANTICS_DIR}/bin/dyn_libs)
      zaghop_add_python_tests(test_quantics_pyrazine QuanticsPyrazineTest
                              Quantics:Pyrazine ENVIRONMENT ${_quantics_env})
      zaghop_add_python_tests(test_quantics_so2 QuanticsSo2Test
                              Quantics:SO2 ENVIRONMENT ${_quantics_env})
      zaghop_add_python_tests(test_quantics_butatriene_db QuanticsButatrieneDbTest
                              Quantics:ButatrieneDB ENVIRONMENT ${_quantics_env})
    endif()
  else()
     message(STATUS "Python not found.")
  endif()

