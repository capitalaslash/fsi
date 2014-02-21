set(HDRS
  ext_pressure.hpp
  pressureramp.hpp
  damp.hpp
)
install(FILES ${HDRS}
  DESTINATION ${CMAKE_INSTALL_PREFIX}/include/bc
)

