cmake_minimum_required(VERSION 3.16)
project(SeqPositionProbs CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_compile_options(-Wall -march=native -ffast-math)
add_compile_options(
        $<$<CONFIG:Debug>:-fsanitize=address,undefined>
        $<$<CONFIG:Debug>:-fno-omit-frame-pointer>
)
add_link_options(
        $<$<CONFIG:Debug>:-fsanitize=address,undefined>
)

find_package(Git QUIET)
if(GIT_FOUND)
    execute_process(
        COMMAND ${GIT_EXECUTABLE} describe --dirty
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_VERSION
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
else()
    set(GIT_VERSION "unknown")
endif()

file(WRITE "${CMAKE_BINARY_DIR}/version.hh.tmp" "\"${GIT_VERSION}\"\n")
execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_BINARY_DIR}/version.hh.tmp" "${CMAKE_SOURCE_DIR}/version.hh")

set(Kokkos_ARCH_NATIVE ON CACHE BOOL "" FORCE)
add_subdirectory(kokkos-5.1.1)

add_executable(dummer dummer.cc dummer-util.hh tantan-wrapper.hh version.hh can_i_haz_simd.hh)
target_compile_definitions(dummer PRIVATE DOUBLE)
target_include_directories(dummer PRIVATE "${CMAKE_SOURCE_DIR}/include")
target_link_libraries(dummer Kokkos::kokkos)

add_executable(dummerl dummer.cc dummer-util.hh tantan-wrapper.hh version.hh can_i_haz_simd.hh)
target_include_directories(dummerl PRIVATE "${CMAKE_SOURCE_DIR}/include")
target_link_libraries(dummerl Kokkos::kokkos)

add_executable(dummer-build dummer-build.cc dummer-util.hh version.hh priors/gap-priors.hh priors/wheeler4.hh priors/blocks9.hh)

set_target_properties(dummer dummerl dummer-build simd_hello PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin"
)
add_executable(test_vectorization test_vectorization.cc)
target_link_libraries(test_vectorization Kokkos::kokkos)
