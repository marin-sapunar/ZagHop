# Run a unit test executable and its cases as a single ctest test (zaghop_add_unit_test in
# CMakeLists.txt):
#   cmake -DEXE=<path> -DN_CASE=<n> -DCASE_1_ARG=<arg> -DCASE_1_REGEX=<regex> ... -P run_cases.cmake
#
# EXE runs once without arguments and passes with exit status 0. Then it runs once with each
# argument CASE_<i>_ARG (i = 1..N_CASE) and passes if its output (standard output and error)
# matches the regular expression CASE_<i>_REGEX, whatever its exit status (as ctest's
# PASS_REGULAR_EXPRESSION). Every run is reported; the output of a failed run is printed, and
# the script fails if any run failed.
cmake_minimum_required(VERSION 3.18)

get_filename_component(_name "${EXE}" NAME)
set(_failed "")

execute_process(COMMAND "${EXE}" RESULT_VARIABLE _res OUTPUT_VARIABLE _out ERROR_VARIABLE _out)
if(_res STREQUAL "0")
  message("PASS  ${_name}")
else()
  message("FAIL  ${_name}: exit status ${_res}\n${_out}")
  list(APPEND _failed "${_name}")
endif()

if(N_CASE GREATER 0)
  foreach(_i RANGE 1 ${N_CASE})
    set(_arg "${CASE_${_i}_ARG}")
    set(_regex "${CASE_${_i}_REGEX}")
    execute_process(COMMAND "${EXE}" "${_arg}"
                    RESULT_VARIABLE _res OUTPUT_VARIABLE _out ERROR_VARIABLE _out)
    if(_out MATCHES "${_regex}")
      message("PASS  ${_name} ${_arg}")
    else()
      message("FAIL  ${_name} ${_arg}: output doesn't match \"${_regex}\"\n${_out}")
      list(APPEND _failed "${_name} ${_arg}")
    endif()
  endforeach()
endif()

# Not if(_failed): a name such as "false" or "off" is a false constant.
list(LENGTH _failed _n_failed)
if(_n_failed GREATER 0)
  list(JOIN _failed ", " _failed)
  message(FATAL_ERROR "Failed: ${_failed}")
endif()
