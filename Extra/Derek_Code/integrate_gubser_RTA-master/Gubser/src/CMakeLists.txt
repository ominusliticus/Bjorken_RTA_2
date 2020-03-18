set (SOURCES
    integrate_gubser_RTA.cpp
    )

add_executable(integrate_gubser_RTA ${SOURCES})
target_link_libraries (integrate_gubser_RTA ${LIBS})

install(TARGETS integrate_gubser_RTA DESTINATION ${CMAKE_HOME_DIRECTORY})
