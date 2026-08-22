cmake_minimum_required(VERSION 3.16)
project(sorting LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib")
set(CMAKE_BUILD_WITH_INSTALL_RPATH OFF)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_LIST_DIR}/.." CACHE PATH "Install prefix" FORCE)
endif()

find_package(Boost REQUIRED COMPONENTS program_options)

find_program(ROOT_CONFIG_EXECUTABLE root-config REQUIRED)
execute_process(
  COMMAND ${ROOT_CONFIG_EXECUTABLE} --cflags
  OUTPUT_VARIABLE ROOT_CFLAGS
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
execute_process(
  COMMAND ${ROOT_CONFIG_EXECUTABLE} --libs
  OUTPUT_VARIABLE ROOT_LIBS
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
separate_arguments(ROOT_CFLAGS_LIST UNIX_COMMAND "${ROOT_CFLAGS}")
separate_arguments(ROOT_LIBS_LIST UNIX_COMMAND "${ROOT_LIBS}")

add_library(timing_estimator SHARED
  timing/src/TimingEstimator.cxx
)
target_include_directories(timing_estimator PUBLIC
  ${CMAKE_CURRENT_LIST_DIR}/timing/include
)

add_executable(decoder decoder.cc)
add_executable(dcalib dcalib.cc)
add_executable(calibrator calibrator.cc)
add_executable(coordinator coordinator.cc)
add_executable(cleaner cleaner.cc)
add_executable(checker checker.cc)
add_executable(sorter sorter.cc)
add_executable(after-pulse-suppressor after-pulse-suppressor.cc)
add_executable(merger merger.cc)
add_executable(trigger trigger.cc)
add_executable(timing timing.cc)
add_executable(ring-finder ring-finder.cc)
add_executable(offset_calibrator offset_calibrator.cc)
add_executable(validate_timing_estimator timing/test/validate_timing_estimator.cc)
target_include_directories(validate_timing_estimator PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/timing/include
)
target_link_libraries(validate_timing_estimator PRIVATE timing_estimator)

foreach(target decoder dcalib calibrator coordinator cleaner checker sorter after-pulse-suppressor merger trigger timing ring-finder offset_calibrator)
  target_compile_options(${target} PRIVATE ${ROOT_CFLAGS_LIST})
  target_link_libraries(${target} PRIVATE ${ROOT_LIBS_LIST} Boost::program_options)
endforeach()
target_include_directories(timing PRIVATE ${CMAKE_CURRENT_LIST_DIR}/timing/include)
target_link_libraries(timing PRIVATE timing_estimator)

install(TARGETS
  timing_estimator
  decoder
  dcalib
  calibrator
  coordinator
  cleaner
  checker
  sorter
  after-pulse-suppressor
  merger
  trigger
  timing
  ring-finder
  offset_calibrator
  validate_timing_estimator
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
)

install(FILES
  timing/include/TimingEstimator.h
  DESTINATION include
)
