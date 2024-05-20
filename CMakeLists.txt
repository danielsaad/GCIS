CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
project("GCIS" CXX C)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


if (CMAKE_VERSION VERSION_LESS "3.1")
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      set (CMAKE_CXX_FLAGS "--std=c++23 -fmodules-ts -fno-permissive ${CMAKE_CXX_FLAGS}")
    endif ()
else ()
    set (CMAKE_CXX_STANDARD 23)
endif ()

# Check if Profile mode is enabled

SET(CMAKE_VERBOSE_MAKEFILE on)
if(CMAKE_BUILD_TYPE MATCHES DEBUG)
	MESSAGE(STATUS "Debug mode")
    SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ")
elseif(CMAKE_BUILD_TYPE MATCHES RELEASE)
    MESSAGE(STATUS "Release mode")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native -Wall")
elseif(${PROFILE})
    MESSAGE(STATUS "Gprof Profile  mode")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg -no-pie")
else()  
    MESSAGE(STATUS "Default mode")
endif(CMAKE_BUILD_TYPE MATCHES DEBUG)

include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/external/mem_monitor)
include_directories(${CMAKE_SOURCE_DIR}/external/malloc_count)
link_directories(${CMAKE_SOURCE_DIR}/lib)

add_subdirectory(${CMAKE_SOURCE_DIR}/src)
# add_subdirectory(${CMAKE_SOURCE_DIR}/test)
add_subdirectory(${CMAKE_SOURCE_DIR}/lib)

#add_subdirectory(${CMAKE_SOURCE_DIR}/benchmark)
