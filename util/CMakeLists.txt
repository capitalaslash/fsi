set(HDRS
  extvtkio.hpp
  init.hpp
  biquad.hpp
  stress.hpp
)
install(FILES ${HDRS}
  DESTINATION ${CMAKE_INSTALL_PREFIX}/include/util
)

set(SRCS
  extvtkio.cpp
  init.cpp
  biquad.cpp
  stress.cpp
)

add_library(fsiutil ${SRCS})

install(TARGETS fsiutil
  DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
)
