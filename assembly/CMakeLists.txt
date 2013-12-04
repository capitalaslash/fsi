set(HDRS
  interface.hpp
  mat.hpp
)
install(FILES ${HDRS}
  DESTINATION ${CMAKE_INSTALL_PREFIX}/include/assembly
)

set(SRCS
  interface.cpp
  mat.cpp
)

add_library(fsiassembly ${SRCS})

install(TARGETS fsiassembly
  DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
)
