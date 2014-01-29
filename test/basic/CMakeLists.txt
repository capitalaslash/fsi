
add_executable(basic_test basic_test.cpp)
target_link_libraries(basic_test ${LIBS})
add_test(test1 basic_test)

add_executable(laplacian_test laplacian_test.cpp)
target_link_libraries(laplacian_test ${LIBS})

add_executable(monomial_2d monomial_2d.cpp)
target_link_libraries(monomial_2d ${LIBS})

add_executable(bidomain bidomain.cpp)
target_link_libraries(bidomain ${LIBS})
