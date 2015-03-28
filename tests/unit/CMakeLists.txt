ADD_DEFINITIONS (-DBOOST_TEST_DYN_LINK)

SET (EXTRA_LIBS ${EXTRA_LIBS} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

SET (dai_UNITTESTS
var_test
smallset_test
varset_test
graph_test
dag_test
bipgraph_test
weightedgraph_test
enum_test
util_test
exceptions_test
properties_test
index_test
prob_test
factor_test
factorgraph_test
clustergraph_test
regiongraph_test
daialg_test
alldai_test
)

FOREACH (UNITTEST ${dai_UNITTESTS})
    ADD_EXECUTABLE (${UNITTEST} ${UNITTEST}.cpp)
    TARGET_LINK_LIBRARIES (${UNITTEST} ${EXTRA_LIBS})
    ADD_TEST (${UNITTEST}_run ${UNITTEST})
ENDFOREACH (UNITTEST)
