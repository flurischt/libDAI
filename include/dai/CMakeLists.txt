SET (dai_hdrs
    alldai.h
    bipgraph.h
    clustergraph.h
    dag.h
    daialg.h
    doc.h
    emalg.h
    enum.h
    evidence.h
    exactinf.h
    exceptions.h
    factor.h
    factorgraph.h
    graph.h
    index.h
    io.h
    prob.h
    properties.h
    regiongraph.h
    smallset.h
    util.h
    var.h
    varset.h
    weightedgraph.h
)

IF(WITH_BP)
    SET(dai_hdrs ${dai_hdrs} bp.h)
ENDIF(WITH_BP)

IF(WITH_FBP)
    SET(dai_hdrs ${dai_hdrs} fbp.h)
ENDIF(WITH_FBP)

IF(WITH_TRWBP)
    SET(dai_hdrs ${dai_hdrs} trwbp.h)
ENDIF(WITH_TRWBP)

IF(WITH_MF)
    SET(dai_hdrs ${dai_hdrs} mf.h)
ENDIF(WITH_MF)

IF(WITH_HAK)
    SET(dai_hdrs ${dai_hdrs} hak.h)
ENDIF(WITH_HAK)

IF(WITH_LC)
    SET(dai_hdrs ${dai_hdrs} lc.h)
ENDIF(WITH_LC)

IF(WITH_TREEEP)
    SET(dai_hdrs ${dai_hdrs} treeep.h)
ENDIF(WITH_TREEEP)

IF(WITH_JTREE)
    SET(dai_hdrs ${dai_hdrs} jtree.h)
ENDIF(WITH_JTREE)

IF(WITH_MR)
    SET(dai_hdrs ${dai_hdrs} mr.h)
ENDIF(WITH_MR)

IF(WITH_GIBBS)
    SET(dai_hdrs ${dai_hdrs} gibbs.h)
ENDIF(WITH_GIBBS)

IF(WITH_CBP)
    SET(dai_hdrs ${dai_hdrs} bbp.h cbp.h bp_dual.h)
ENDIF(WITH_CBP)

IF(WITH_DECMAP)
    SET(dai_hdrs ${dai_hdrs} decmap.h)
ENDIF(WITH_DECMAP)

INSTALL (FILES ${dai_hdrs} DESTINATION include/dai COMPONENT main)
