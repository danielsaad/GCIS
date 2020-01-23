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
        ../include/gcis_eliasfano_index.hpp
        ../include/grammar_builder.hpp
        ../include/index_builder.hpp
        ../include/eliasfano.hpp
        ../include/gcis_unary.hpp
        ../include/gcis_s8b.hpp
        ../include/gcis_eliasfano.hpp
        ../include/gcis_eliasfano_no_lcp.hpp
        ../include/gcis_gap_codec.hpp
        ../include/gcis_gap.hpp
        ../include/sais_nong.hpp
        ../include/gcis_index.hpp
        ../include/gcis_index_bs.hpp
        ../include/gcis_index_grid.hpp
        ../include/gcis_index_dfuds.hpp)

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
        ../include/gcis_eliasfano_index.hpp
        ../include/grammar_builder.hpp
        ../include/index_builder.hpp
        ../include/eliasfano.hpp
        ../include/gcis_unary.hpp
        ../include/gcis_s8b.hpp
        ../include/gcis_eliasfano.hpp
        ../include/gcis_eliasfano_no_lcp.hpp
        ../include/gcis_gap_codec.hpp
        ../include/gcis_gap.hpp
        ../include/sais_nong.hpp
        ../include/gcis_index.hpp
        ../include/gcis_index_bs.hpp
        ../include/gcis_index_grid.hpp
        ../include/gcis_index_dfuds.hpp
        )

target_compile_definitions(gc-is-statistics PRIVATE  MEM_MONITOR REPORT )
install (TARGETS gc-is ARCHIVE DESTINATION ${CMAKE_SOURCE_DIR}/lib)
install (TARGETS gc-is-statistics ARCHIVE DESTINATION ${CMAKE_SOURCE_DIR}/lib)