include(ExternalProject)

### Git clone google test library and add to our project
ExternalProject_Add(
  celeroproject
  GIT_REPOSITORY https://github.com/DigitalInBlue/Celero.git
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  CMAKE_ARGS "-DCELERO_COMPILE_DYNAMIC_LIBRARIES=OFF"
  LOG_DOWNLOAD ON
  LOG_CONFIGURE ON
  LOG_BUILD ON)

ExternalProject_Get_Property(celeroproject source_dir)
set(CELERO_INCLUDE_DIRS ${source_dir}/include)


ExternalProject_Get_Property(celeroproject binary_dir)
set(CELERO_LIBRARY_PATH ${binary_dir}/${CMAKE_FIND_LIBRARY_PREFIXES}celero.a)
set(CELERO_LIBRARY celero)
add_library(${CELERO_LIBRARY} UNKNOWN IMPORTED)
set_target_properties(${CELERO_LIBRARY} PROPERTIES
  IMPORTED_LOCATION ${CELERO_LIBRARY_PATH}
)
add_dependencies(${CELERO_LIBRARY} celeroproject)

# Including header files in our project
include_directories(${CELERO_INCLUDE_DIRS})

add_executable(benchmark_main
    benchmark_main.cpp)
target_link_libraries(benchmark_main celero gc-is sdsl)
install(TARGETS benchmark_main RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/bin)
