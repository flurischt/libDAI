SET(dai_srcs
    graph.cpp
    bipgraph.cpp
    varset.cpp
    daialg.cpp
    alldai.cpp
    factor.cpp
    factorgraph.cpp
    properties.cpp
    util.cpp
    exceptions.cpp
    exactinf.cpp
)

SET(dai_include_dir "../include/dai")
SET(dai_hdrs
    "${dai_include_dir}/alldai.h"
    "${dai_include_dir}/bipgraph.h"
    "${dai_include_dir}/daialg.h"
    "${dai_include_dir}/dai_config.h"
    "${dai_include_dir}/doc.h"
    "${dai_include_dir}/enum.h"
    "${dai_include_dir}/exactinf.h"
    "${dai_include_dir}/exceptions.h"
    "${dai_include_dir}/factor.h"
    "${dai_include_dir}/factorgraph.h"
    "${dai_include_dir}/graph.h"
    "${dai_include_dir}/index.h"
    "${dai_include_dir}/prob.h"
    "${dai_include_dir}/properties.h"
    "${dai_include_dir}/smallset.h"
    "${dai_include_dir}/util.h"
    "${dai_include_dir}/var.h"
    "${dai_include_dir}/varset.h"
    "${dai_include_dir}/utils/timer.h"
)

# Files used only for inference methods other than BP, FBP or TRWBP
# TODO: fix this properly
IF(WITH_ADVANCED)
    ADD_DEFINITIONS("-DDAI_WITH_ADVANCED")
    LIST(APPEND dai_hdrs
         "${dai_include_dir}/clustergraph.h"
         "${dai_include_dir}/dag.h"
         "${dai_include_dir}/emalg.h"
         "${dai_include_dir}/evidence.h"
         "${dai_include_dir}/io.h"
         "${dai_include_dir}/regiongraph.h"
    )
    LIST(APPEND dai_srcs
         clustergraph.cpp
         dag.cpp
         emalg.cpp
         evidence.cpp
         io.cpp
         regiongraph.cpp
    )
ENDIF(WITH_ADVANCED)

IF(WITH_BP)
    SET(dai_hdrs ${dai_hdrs} "${dai_include_dir}/bp.h")
    SET(dai_srcs ${dai_srcs} bp.cpp)
    SET(dai_srcs ${dai_srcs} bp_ext.cpp)
ENDIF(WITH_BP)

IF(WITH_FBP)
    SET(dai_hdrs ${dai_hdrs} "${dai_include_dir}/fbp.h")
    SET(dai_srcs ${dai_srcs} fbp.cpp)
ENDIF(WITH_FBP)

IF(WITH_TRWBP)
    SET(dai_hdrs ${dai_hdrs}
        "${dai_include_dir}/trwbp.h"
        "${dai_include_dir}/weightedgraph.h"
    )
    SET(dai_srcs ${dai_srcs}
        trwbp.cpp weightedgraph.cpp)
ENDIF(WITH_TRWBP)

IF(WITH_MF)
    SET(dai_hdrs ${dai_hdrs} "${dai_include_dir}/mf.h")
    SET(dai_srcs ${dai_srcs} mf.cpp)
ENDIF(WITH_MF)

IF(WITH_HAK)
    SET(dai_hdrs ${dai_hdrs} "${dai_include_dir}/hak.h")
    SET(dai_srcs ${dai_srcs} hak.cpp)
ENDIF(WITH_HAK)

IF(WITH_LC)
    SET(dai_hdrs ${dai_hdrs} "${dai_include_dir}/lc.h")
    SET(dai_srcs ${dai_srcs} lc.cpp)
ENDIF(WITH_LC)

IF(WITH_TREEEP)
    SET(dai_hdrs ${dai_hdrs} "${dai_include_dir}/treeep.h")
    SET(dai_srcs ${dai_srcs} treeep.cpp)
ENDIF(WITH_TREEEP)

IF(WITH_JTREE)
    SET(dai_hdrs ${dai_hdrs} "${dai_include_dir}/jtree.h")
    SET(dai_srcs ${dai_srcs} jtree.cpp)
ENDIF(WITH_JTREE)

IF(WITH_MR)
    SET(dai_hdrs ${dai_hdrs} "${dai_include_dir}/mr.h")
    SET(dai_srcs ${dai_srcs} mr.cpp)
ENDIF(WITH_MR)

IF(WITH_GIBBS)
    SET(dai_hdrs ${dai_hdrs} "${dai_include_dir}/gibbs.h")
    SET(dai_srcs ${dai_srcs} gibbs.cpp)
ENDIF(WITH_GIBBS)

IF(WITH_CBP)
    SET(dai_hdrs ${dai_hdrs}
        "${dai_include_dir}/bbp.h"
        "${dai_include_dir}/cbp.h"
        "${dai_include_dir}/bp_dual.h")
    SET(dai_srcs ${dai_srcs} bbp.cpp cbp.cpp bp_dual.cpp)
ENDIF(WITH_CBP)

IF(WITH_DECMAP)
    SET(dai_hdrs ${dai_hdrs} "${dai_include_dir}/decmap.h")
    SET(dai_srcs ${dai_srcs} decmap.cpp)
ENDIF(WITH_DECMAP)

IF(WITH_GLC)
    SET(dai_hdrs ${dai_hdrs}
        "${dai_include_dir}/glc.h"
        "${dai_include_dir}/cobwebgraph.h")
    SET(dai_srcs ${dai_srcs}
        glc.cpp
        cobwebgraph.cpp)
ENDIF(WITH_GLC)

ADD_LIBRARY(dai ${dai_hdrs} ${dai_srcs})
TARGET_LINK_LIBRARIES(dai ${EXTRA_LIBS})

INSTALL(TARGETS dai
    RUNTIME DESTINATION bin COMPONENT main
    LIBRARY DESTINATION lib COMPONENT main
    ARCHIVE DESTINATION lib COMPONENT main
)

#set_target_properties(dai PROPERTIES COMPILE_FLAGS "-S")
#set_target_properties(dai PROPERTIES COMPILE_FLAGS "-save-temps")

INSTALL(FILES ${dai_hdrs} DESTINATION include/dai COMPONENT main)
