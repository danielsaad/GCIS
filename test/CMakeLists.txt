include(FetchContent)
FetchContent_Declare(
  googletest
  # Specify the commit you depend on and update it regularly>.
  URL https://github.com/google/googletest/releases/download/v1.15.2/googletest-1.15.2.tar.gz)
 
link_directories(${CMAKE_SOURCE_DIR}/lib)
add_executable(test_main test_main.cpp test_simple8b.cpp ../lib/gcis_s8b_codec.cpp test_eliasfano.cpp)
target_link_libraries(test_main gtest gc-is sdsl pthread)
install(TARGETS test_main RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/bin)