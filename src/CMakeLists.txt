add_executable(gcis gc-is-codec.cpp)
add_executable(gcis-memory gc-is-codec.cpp ${CMAKE_SOURCE_DIR}/external/malloc_count/malloc_count.c)
add_executable(gcis-64 gc-is-codec.cpp)
add_executable(gcis-64-memory gc-is-codec.cpp ${CMAKE_SOURCE_DIR}/external/malloc_count/malloc_count.c)

target_compile_definitions(gcis-64 PRIVATE m64 )
target_compile_definitions(gcis-memory PRIVATE MEM_MONITOR REPORT )
target_compile_definitions(gcis-64-memory PRIVATE MEM_MONITOR REPORT m64 )

target_link_libraries(gcis gc-is sdsl pthread sais )
target_link_libraries(gcis-memory gc-is-statistics sdsl pthread dl sais)
target_link_libraries(gcis-64 gc-is-64 sdsl pthread sais )
target_link_libraries(gcis-64-memory  gc-is-statistics-64 sdsl pthread dl sais)



add_executable(sais-yuta sais-yuta.cpp)
add_executable(sais-lcp-yuta sais-lcp-yuta.cpp)
add_executable(sais-nong sais-nong.cpp)
add_executable(decode-sais-nong decode-sais-nong.cpp)
add_executable(decode-sais-yuta decode-sais-yuta.cpp)
add_executable(decode-sais-lcp-yuta decode-sais-lcp-yuta.cpp)
add_executable(sais-divsufsort divsufsort.cpp)
add_executable(sais-divsufsort-lcp divsufsort-lcp.cpp)
add_executable(decode-sais-divsufsort decode-divsufsort.cpp)
add_executable(decode-sais-divsufsort-lcp decode-divsufsort-lcp.cpp)

target_link_libraries(sais-yuta  sais sdsl)
target_link_libraries(sais-lcp-yuta saislcp sdsl)
target_link_libraries(sais-nong gc-is sdsl)
target_link_libraries(decode-sais-nong gc-is sdsl )
target_link_libraries(decode-sais-yuta gc-is sdsl sais)
target_link_libraries(decode-sais-lcp-yuta gc-is sdsl saislcp)
target_link_libraries(sais-divsufsort divsufsort)
target_link_libraries(sais-divsufsort-lcp divsufsort-lcp)
target_link_libraries(sais-divsufsort-lcp divsufsort-lcp)
target_link_libraries(decode-sais-divsufsort gc-is sdsl divsufsort)
target_link_libraries(decode-sais-divsufsort-lcp gc-is sdsl divsufsort-lcp)


install(TARGETS sais-nong sais-yuta sais-lcp-yuta decode-sais-nong decode-sais-yuta decode-sais-lcp-yuta gcis gcis-memory gcis-64 gcis-64-memory sais-divsufsort sais-divsufsort-lcp decode-sais-divsufsort decode-sais-divsufsort-lcp  RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/bin)
