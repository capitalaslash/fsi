set(HDRS
  extvtkio.hpp
)
install(FILES ${HDRS}
  DESTINATION ${CMAKE_INSTALL_PREFIX}/include
)

set(SRCS
  extvtkio.cpp
)

add_library(fsiutil ${SRCS} ${LIBS})

install(TARGETS fsiutil
  DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
)

