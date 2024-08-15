ADD_LIBRARY(gc-is STATIC
        util.cpp
        simple8b.cpp
        gcis_s8b_codec.cpp
        gcis_unary_codec.cpp
        gcis_eliasfano_codec.cpp
        gcis_eliasfano_codec_no_lcp.cpp
        gcis_gap_codec.cpp
        eliasfano.cpp
        sais_nong.cpp
        ../include/grammar_builder.hpp
        ../include/eliasfano.hpp
        ../include/gcis_unary.hpp
        ../include/gcis_s8b.hpp
        ../include/gcis_eliasfano.hpp
        ../include/gcis_eliasfano_no_lcp.hpp
        ../include/gcis_gap_codec.hpp
        ../include/gcis_gap.hpp
        ../include/sais_nong.hpp
)


ADD_LIBRARY(gc-is-64 STATIC
        util.cpp
        simple8b.cpp
        gcis_s8b_codec.cpp
        gcis_unary_codec.cpp
        gcis_eliasfano_codec.cpp
        gcis_eliasfano_codec_no_lcp.cpp
        gcis_gap_codec.cpp
        eliasfano.cpp
        sais_nong.cpp
        ../include/grammar_builder.hpp
        ../include/eliasfano.hpp
        ../include/gcis_unary.hpp
        ../include/gcis_s8b.hpp
        ../include/gcis_eliasfano.hpp
        ../include/gcis_eliasfano_no_lcp.hpp
        ../include/gcis_gap_codec.hpp
        ../include/gcis_gap.hpp
        ../include/sais_nong.hpp
)

ADD_LIBRARY(gc-is-statistics STATIC
        util.cpp
        simple8b.cpp
        gcis_s8b_codec.cpp
        gcis_unary_codec.cpp
        gcis_eliasfano_codec.cpp
        gcis_eliasfano_codec_no_lcp.cpp
        gcis_gap_codec.cpp
        eliasfano.cpp
        sais_nong.cpp
        ../include/grammar_builder.hpp
        ../include/eliasfano.hpp
        ../include/gcis_unary.hpp
        ../include/gcis_s8b.hpp
        ../include/gcis_eliasfano.hpp
        ../include/gcis_eliasfano_no_lcp.hpp
        ../include/gcis_gap_codec.hpp
        ../include/gcis_gap.hpp
        ../include/sais_nong.hpp
)

ADD_LIBRARY(gc-is-statistics-64 STATIC
        util.cpp
        simple8b.cpp
        gcis_s8b_codec.cpp
        gcis_unary_codec.cpp
        gcis_eliasfano_codec.cpp
        gcis_eliasfano_codec_no_lcp.cpp
        gcis_gap_codec.cpp
        eliasfano.cpp
        sais_nong.cpp
        ../include/grammar_builder.hpp
        ../include/eliasfano.hpp
        ../include/gcis_unary.hpp
        ../include/gcis_s8b.hpp
        ../include/gcis_eliasfano.hpp
        ../include/gcis_eliasfano_no_lcp.hpp
        ../include/gcis_gap_codec.hpp
        ../include/gcis_gap.hpp
        ../include/sais_nong.hpp
)

target_compile_definitions(gc-is-64 PRIVATE m64)
target_compile_definitions(gc-is-statistics PRIVATE  MEM_MONITOR REPORT )
target_compile_definitions(gc-is-statistics-64 PRIVATE m64 MEM_MONITOR REPORT )
install (TARGETS gc-is ARCHIVE DESTINATION ${CMAKE_SOURCE_DIR}/lib)
install (TARGETS gc-is-statistics ARCHIVE DESTINATION ${CMAKE_SOURCE_DIR}/lib)
install (TARGETS gc-is-64 ARCHIVE DESTINATION ${CMAKE_SOURCE_DIR}/lib)
install (TARGETS gc-is-statistics-64 ARCHIVE DESTINATION ${CMAKE_SOURCE_DIR}/lib)