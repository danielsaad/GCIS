include(ExternalProject)

### Git clone google test library and add to our project
ExternalProject_Add(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  LOG_DOWNLOAD ON
  LOG_CONFIGURE ON
  LOG_BUILD ON)

ExternalProject_Get_Property(googletest source_dir binary_dir)

set(GTEST_INCLUDE_DIRS ${source_dir}/googletest/include)
set(GMOCK_INCLUDE_DIRS ${source_dir}/googlemock/include)



ExternalProject_Get_Property(googletest binary_dir)
set(GTEST_LIBRARY_PATH ${binary_dir}/lib/${CMAKE_FIND_LIBRARY_PREFIXES}gtest.a)
set(GTEST_LIBRARY gtest)
add_library(${GTEST_LIBRARY} UNKNOWN IMPORTED)
set_target_properties(${GTEST_LIBRARY} PROPERTIES
  IMPORTED_LOCATION ${GTEST_LIBRARY_PATH}
)
add_dependencies(${GTEST_LIBRARY} googletest)

set(GTEST_MAIN_LIBRARY_PATH ${binary_dir}/lib/${CMAKE_FIND_LIBRARY_PREFIXES}gtest_main.a)
set(GTEST_MAIN_LIBRARY gtest_main)
add_library(${GTEST_MAIN_LIBRARY} UNKNOWN IMPORTED)
set_target_properties(${GTEST_MAIN_LIBRARY} PROPERTIES
  IMPORTED_LOCATION ${GTEST_MAIN_LIBRARY_PATH}
)
add_dependencies(${GTEST_MAIN_LIBRARY} googletest)

set(GMOCK_LIBRARY_PATH ${binary_dir}/lib/${CMAKE_FIND_LIBRARY_PREFIXES}gmock.a)
set(GMOCK_LIBRARY gmock)
add_library(${GMOCK_LIBRARY} UNKNOWN IMPORTED)
set_target_properties(${GMOCK_LIBRARY} PROPERTIES
  IMPORTED_LOCATION ${GMOCK_LIBRARY_PATH}
)
add_dependencies(${GMOCK_LIBRARY} googletest)

set(GMOCK_MAIN_LIBRARY_PATH ${binary_dir}/lib/${CMAKE_FIND_LIBRARY_PREFIXES}gmock_main.a)
set(GMOCK_MAIN_LIBRARY gmock_main)
add_library(${GMOCK_MAIN_LIBRARY} UNKNOWN IMPORTED)
set_target_properties(${GMOCK_MAIN_LIBRARY} PROPERTIES
  IMPORTED_LOCATION ${GMOCK_MAIN_LIBRARY_PATH}
)
add_dependencies(${GMOCK_MAIN_LIBRARY} ${GTEST_LIBRARY})

# Including header files in our project
include_directories(${GTEST_INCLUDE_DIRS})
include_directories(${GMOCK_INCLUDE_DIRS})

add_executable(test_main test_main.cpp test_simple8b.cpp ../lib/gcis_s8b_codec.cpp test_eliasfano.cpp)
add_executable(test_gcis_index_grid test_gcis_index_grid.cpp)
add_executable(test_gcis_index_bs test_gcis_index_bs.cpp)
add_executable(test_gcis_index test_gcis_index.cpp)

link_directories(${CMAKE_SOURCE_DIR}/lib)
target_link_libraries(test_main gtest gc-is sdsl gmock pthread)
target_link_libraries(test_gcis_index_grid gtest gtest_main gc-is sdsl divsufsort divsufsort64 gmock pthread)
target_link_libraries(test_gcis_index gtest gtest_main gc-is sdsl divsufsort divsufsort64 gmock pthread)
target_link_libraries(test_gcis_index_bs gtest gtest_main gc-is sdsl divsufsort divsufsort64  gmock pthread)

# install(TARGETS test_main RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/bin)
