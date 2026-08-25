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

include(CheckLanguage)
check_language(CUDA)
set(RING_HOUGH_CUDA_ENABLED OFF)
if(CMAKE_CUDA_COMPILER)
  enable_language(CUDA)
  if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    # The supported deployment GPU is the Turing-generation T550 (sm_75).
    # Users can override this cache variable for another CUDA device.
    set(CMAKE_CUDA_ARCHITECTURES 75 CACHE STRING "CUDA architectures")
  endif()
  find_package(CUDAToolkit QUIET)
  if(CUDAToolkit_FOUND)
    set(RING_HOUGH_CUDA_ENABLED ON)
    add_library(ring_hough_cuda STATIC ring-finder-hough-cuda.cu)
    set_target_properties(ring_hough_cuda PROPERTIES
      CUDA_STANDARD 17
      CUDA_STANDARD_REQUIRED ON
    )
    target_include_directories(ring_hough_cuda PUBLIC ${CMAKE_CURRENT_LIST_DIR})
    target_link_libraries(ring_hough_cuda PUBLIC CUDA::cudart)
  endif()
endif()

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
add_executable(ring-finder-ransac ring-finder-ransac.cc)
add_executable(ring-finder-hough ring-finder-hough.cc)
add_executable(validate_timing_estimator timing/test/validate_timing_estimator.cc)
target_include_directories(validate_timing_estimator PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/timing/include
)
target_link_libraries(validate_timing_estimator PRIVATE timing_estimator)

foreach(target decoder dcalib calibrator coordinator cleaner checker sorter after-pulse-suppressor merger trigger timing ring-finder-ransac ring-finder-hough)
  target_compile_options(${target} PRIVATE ${ROOT_CFLAGS_LIST})
  target_link_libraries(${target} PRIVATE ${ROOT_LIBS_LIST} Boost::program_options)
endforeach()
target_include_directories(timing PRIVATE ${CMAKE_CURRENT_LIST_DIR}/timing/include)
target_link_libraries(timing PRIVATE timing_estimator)
if(RING_HOUGH_CUDA_ENABLED)
  target_compile_definitions(ring-finder-hough PRIVATE RING_HOUGH_HAS_CUDA=1)
  target_link_libraries(ring-finder-hough PRIVATE ring_hough_cuda)
endif()

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
  ring-finder-ransac
  ring-finder-hough
  validate_timing_estimator
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
)

install(FILES
  timing/include/TimingEstimator.h
  DESTINATION include
)
