cmake_minimum_required(VERSION 3.16)
project(sorting LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

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

add_executable(calibrator calibrator.cc)
add_executable(sorter sorter.cc)
add_executable(after-pulse-suppressor after-pulse-suppressor.cc)
add_executable(merger merger.cc)
add_executable(trigger trigger.cc)

foreach(target calibrator sorter after-pulse-suppressor merger trigger)
  target_compile_options(${target} PRIVATE ${ROOT_CFLAGS_LIST})
  target_link_libraries(${target} PRIVATE ${ROOT_LIBS_LIST} Boost::program_options)
endforeach()

install(TARGETS
  calibrator
  sorter
  after-pulse-suppressor
  merger
  trigger
  RUNTIME DESTINATION bin
)
