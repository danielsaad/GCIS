ADD_LIBRARY(gc-is STATIC util.cpp simple8b.cpp gcis_s8b_codec.cpp gcis_unary_codec.cpp gcis_eliasfano_codec.cpp ../include/eliasfano.hpp eliasfano.cpp ../include/gcis_unary.hpp  ../include/gcis_s8b.hpp ../include/gcis_eliasfano.hpp)
ADD_LIBRARY(gc-is-statistics STATIC util.cpp simple8b.cpp gcis_s8b_codec.cpp gcis_unary_codec.cpp gcis_eliasfano_codec.cpp eliasfano.cpp  ../include/gcis_unary.hpp ../include/gcis_s8b.hpp ../include/gcis_eliasfano.hpp)
target_compile_definitions(gc-is-statistics PRIVATE  MEM_MONITOR REPORT )
install (TARGETS gc-is ARCHIVE DESTINATION ${CMAKE_SOURCE_DIR}/lib)
