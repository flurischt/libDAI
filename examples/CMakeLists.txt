IF(WITH_ADVANCED)
    SET(dai_EXAMPLES
        example
        example_bipgraph
        example_varset
        example_permute
        example_sprinkler
       )


    LIST(APPEND dai_EXAMPLES example_sprinkler_em)

    IF(WITH_GIBBS)
        SET(dai_EXAMPLES ${dai_EXAMPLES} example_sprinkler_gibbs)
    ENDIF()

    FOREACH(EXAMPLE ${dai_EXAMPLES})
        ADD_EXECUTABLE(${EXAMPLE} ${EXAMPLE}.cpp)
        TARGET_LINK_LIBRARIES(${EXAMPLE} ${EXTRA_LIBS})
        INSTALL(TARGETS ${EXAMPLE} DESTINATION share/dai/examples)
    ENDFOREACH(EXAMPLE)

    ADD_EXECUTABLE(uai2010-aie-solver uai2010-aie-solver.cpp)
    TARGET_LINK_LIBRARIES(uai2010-aie-solver ${EXTRA_LIBS})
    INSTALL(TARGETS uai2010-aie-solver DESTINATION share/dai/examples)

    SET(dai_EXAMPLE_RESOURCES
        example_sprinkler.dot
        example_sprinkler.png
        sprinkler.em
        example_bipgraph.out
        example_permute.out
        example_varset.out
       )
ENDIF(WITH_ADVANCED)


IF(WITH_CIMG)
    FIND_PACKAGE(CImg)
    IF(NEW_CIMG)
        ADD_DEFINITIONS("-DNEW_CIMG")
    ENDIF()

    INCLUDE_DIRECTORIES("${CIMG_INCLUDE_DIRS}")
    SET(EXTRA_LIBS ${EXTRA_LIBS} ${CIMG_LIBRARIES})

    SET(dai_EXAMPLE_RESOURCES ${dai_EXAMPLE_RESOURCES}
        example_img_in1.jpg
        example_img_in2.jpg
    )

    ADD_EXECUTABLE(example_imagesegmentation example_imagesegmentation.cpp)
    TARGET_LINK_LIBRARIES(example_imagesegmentation ${EXTRA_LIBS})
    INSTALL(TARGETS example_imagesegmentation DESTINATION share/dai/examples)

    # For convenience, copy default input images into binary directory
    ADD_CUSTOM_COMMAND(TARGET example_imagesegmentation POST_BUILD
                       COMMAND cp "${CMAKE_SOURCE_DIR}/examples/example_img_in*.jpg"
                                  "${CMAKE_CURRENT_BINARY_DIR}/"
                      )

ENDIF()

ADD_EXECUTABLE(example_recommendation example_recommendation.cpp)
TARGET_LINK_LIBRARIES(example_recommendation ${EXTRA_LIBS})
INSTALL(TARGETS example_recommendation DESTINATION share/dai/examples)
ADD_CUSTOM_COMMAND(TARGET example_recommendation POST_BUILD
                   COMMAND cp "${CMAKE_SOURCE_DIR}/examples/ml-100k/*"
                              "${CMAKE_CURRENT_BINARY_DIR}/"
                  )
ADD_CUSTOM_COMMAND(TARGET example_recommendation POST_BUILD
                   COMMAND cp "${CMAKE_SOURCE_DIR}/examples/recommendation_ref_files/*"
                              "${CMAKE_CURRENT_BINARY_DIR}/"
                  )
#set_target_properties(example_recommendation PROPERTIES COMPILE_FLAGS "-save-temps")
#set_target_properties(example_recommendation PROPERTIES COMPILE_FLAGS "-S")

INSTALL(FILES ${dai_EXAMPLE_RESOURCES} DESTINATION share/dai/examples)
