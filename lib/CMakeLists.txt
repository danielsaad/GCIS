ADD_LIBRARY(gc-is STATIC util.cpp simple8b.cpp gcis_s8b_codec.cpp gcis_unary_codec.cpp gcis_eliasfano_codec.cpp ../include/eliasfano.hpp eliasfano.cpp ../src/gcis-extract.cpp)
ADD_LIBRARY(gc-is-statistics STATIC util.cpp simple8b.cpp gcis_s8b_codec.cpp gcis_unary_codec.cpp gcis_eliasfano_codec.cpp eliasfano.cpp)
target_compile_definitions(gc-is-statistics PRIVATE  MEM_MONITOR REPORT )
install (TARGETS gc-is ARCHIVE DESTINATION ${CMAKE_SOURCE_DIR}/lib)
