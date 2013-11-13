set(HDRS
  extvtkio.hpp
  init.hpp
)
install(FILES ${HDRS}
  DESTINATION ${CMAKE_INSTALL_PREFIX}/include
)

set(SRCS
  extvtkio.cpp
  init.cpp
)

add_library(fsiutil ${SRCS})

install(TARGETS fsiutil
  DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
)
