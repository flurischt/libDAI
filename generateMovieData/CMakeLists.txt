ADD_EXECUTABLE(generate_inputfiles generate.cpp)
TARGET_LINK_LIBRARIES(generate_inputfiles ${EXTRA_LIBS})
INSTALL(TARGETS generate_inputfiles DESTINATION share/dai/examples)
INSTALL(FILES ${dai_EXAMPLE_RESOURCES} DESTINATION share/dai/examples)
