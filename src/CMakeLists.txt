add_executable(gc-is-codec gc-is-codec.cpp)
target_link_libraries(gc-is-codec gc-is sdsl pthread sais)
add_executable(gc-is-codec-memory gc-is-codec.cpp ${CMAKE_SOURCE_DIR}/external/malloc_count/malloc_count.c)
target_link_libraries(gc-is-codec-memory gc-is-statistics sdsl pthread dl sais)
target_compile_definitions(gc-is-codec-memory PRIVATE MEM_MONITOR REPORT )

add_executable(sais-yuta sais-yuta.cpp)
target_link_libraries(sais-yuta gc-is sdsl pthread sais)

add_executable(sais-nong sais-nong.cpp)

install(TARGETS sais-nong sais-yuta gc-is-codec gc-is-codec-memory RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/bin)
