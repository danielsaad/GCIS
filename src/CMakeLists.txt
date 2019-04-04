add_executable(gc-is-codec gc-is-codec.cpp)
target_link_libraries(gc-is-codec gc-is sdsl pthread sais)
add_executable(gc-is-codec-memory gc-is-codec.cpp ${CMAKE_SOURCE_DIR}/external/malloc_count/malloc_count.c)
target_link_libraries(gc-is-codec-memory gc-is-statistics sdsl pthread dl sais)
target_compile_definitions(gc-is-codec-memory PRIVATE MEM_MONITOR REPORT )

add_executable(sais-yuta sais-yuta.cpp)
add_executable(sais-nong sais-nong.cpp)
add_executable(decode-sais-nong decode-sais-nong.cpp)
add_executable(decode-sais-yuta decode-sais-yuta.cpp)
add_executable(sais-divsufsort divsufsort.cpp)
add_executable(decode-sais-divsufsort decode-divsufsort.cpp)

target_link_libraries(sais-yuta gc-is sdsl pthread sais)
target_link_libraries(sais-nong gc-is sdsl )
target_link_libraries(decode-sais-nong gc-is sdsl )
target_link_libraries(decode-sais-yuta gc-is sdsl sais)
target_link_libraries(sais-divsufsort divsufsort)
target_link_libraries(decode-sais-divsufsort gc-is sdsl divsufsort)


install(TARGETS sais-nong sais-yuta decode-sais-nong decode-sais-yuta gc-is-codec gc-is-codec-memory sais-divsufsort decode-sais-divsufsort RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/bin)
