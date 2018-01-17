add_executable(gc-is-codec gc-is-codec.cpp)
target_link_libraries(gc-is-codec gc-is sdsl pthread)
add_executable(gc-is-codec-memory gc-is-codec.cpp)
target_link_libraries(gc-is-codec-memory gc-is-statistics sdsl pthread)
target_compile_definitions(gc-is-codec-memory PRIVATE MEM_MONITOR REPORT )

install(TARGETS gc-is-codec gc-is-codec-memory RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/bin)
