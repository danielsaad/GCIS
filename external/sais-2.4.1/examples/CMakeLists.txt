## Targets ##
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../include"
                    "${CMAKE_CURRENT_BINARY_DIR}/../include")
link_directories("${CMAKE_CURRENT_BINARY_DIR}/../lib")
foreach(src suftest bwt unbwt)
  add_executable(${src} ${src}.c)
  target_link_libraries(${src} sais)
endforeach(src)

if(BUILD_SAIS64)
  foreach(src suftest64 bwt64 unbwt64)
    add_executable(${src} ${src}.c)
    target_link_libraries(${src} sais64)
  endforeach(src)
endif(BUILD_SAIS64)
