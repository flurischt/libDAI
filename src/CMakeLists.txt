SET (dai_srcs
    graph.cpp
    dag.cpp
    bipgraph.cpp
    varset.cpp
    daialg.cpp
    alldai.cpp
    clustergraph.cpp
    factor.cpp
    factorgraph.cpp
    properties.cpp
    regiongraph.cpp
    util.cpp
    weightedgraph.cpp
    exceptions.cpp
    exactinf.cpp
    evidence.cpp
    emalg.cpp
    io.cpp
)

IF(WITH_BP)
    SET(dai_srcs ${dai_srcs} bp.cpp)
ENDIF(WITH_BP)

IF(WITH_FBP)
    SET(dai_srcs ${dai_srcs} fbp.cpp)
ENDIF(WITH_FBP)

IF(WITH_TRWBP)
    SET(dai_srcs ${dai_srcs} trwbp.cpp)
ENDIF(WITH_TRWBP)

IF(WITH_MF)
    SET(dai_srcs ${dai_srcs} mf.cpp)
ENDIF(WITH_MF)

IF(WITH_HAK)
    SET(dai_srcs ${dai_srcs} hak.cpp)
ENDIF(WITH_HAK)

IF(WITH_LC)
    SET(dai_srcs ${dai_srcs} lc.cpp)
ENDIF(WITH_LC)

IF(WITH_TREEEP)
    SET(dai_srcs ${dai_srcs} treeep.cpp)
ENDIF(WITH_TREEEP)

IF(WITH_JTREE)
    SET(dai_srcs ${dai_srcs} jtree.cpp)
ENDIF(WITH_JTREE)

IF(WITH_MR)
    SET(dai_srcs ${dai_srcs} mr.cpp)
ENDIF(WITH_MR)

IF(WITH_GIBBS)
    SET(dai_srcs ${dai_srcs} gibbs.cpp)
ENDIF(WITH_GIBBS)

IF(WITH_CBP)
    SET(dai_srcs ${dai_srcs} bbp.cpp cbp.cpp bp_dual.cpp)
ENDIF(WITH_CBP)

IF(WITH_DECMAP)
    SET(dai_srcs ${dai_srcs} decmap.cpp)
ENDIF(WITH_DECMAP)


ADD_LIBRARY (dai ${dai_srcs})
TARGET_LINK_LIBRARIES (dai ${EXTRA_LIBS})

INSTALL (TARGETS dai
    RUNTIME DESTINATION bin COMPONENT main
    LIBRARY DESTINATION lib COMPONENT main
    ARCHIVE DESTINATION lib COMPONENT main
)
