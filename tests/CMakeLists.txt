ADD_EXECUTABLE (testdai testdai.cpp)
TARGET_LINK_LIBRARIES (testdai ${EXTRA_LIBS} ${Boost_PROGRAM_OPTIONS_LIBRARY})
INSTALL (TARGETS testdai DESTINATION share/dai/tests)

ADD_EXECUTABLE (testem testem/testem.cpp)
TARGET_LINK_LIBRARIES (testem ${EXTRA_LIBS} ${Boost_PROGRAM_OPTIONS_LIBRARY})
INSTALL (TARGETS testem DESTINATION share/dai/tests/testem)

IF(WITH_CBP)
    ADD_EXECUTABLE (testbbp testbbp.cpp)
    TARGET_LINK_LIBRARIES (testbbp ${EXTRA_LIBS})
    INSTALL (TARGETS testbbp DESTINATION share/dai/tests)
ENDIF(WITH_CBP)

SET (dai_TESTS_RESOURCES
alarm.fg
aliases.conf
hoi1.fg
hoi2.fg
hoi3.fg
hoi4.fg
jtreemapbug.fg
maxprodbug.fg
maxprodbug2.fg
maxprodbug3.fg
testall
testall.bat
testfast.fg
testfast.out
testregression
testregression.bat
twofactors.fg
zeroes1.fg
)

SET (dai_TESTS_TESTEM_RESOURCES
testem/2var.em
testem/2var.fg
testem/2var_data.tab
testem/3var.em
testem/3var.fg
testem/hoi1_data.tab
testem/hoi1_infer_f2.em
testem/hoi1_share_f0_f1_f2.em
testem/hoi1_share_f0_f2.em
testem/runtests
testem/runtests.bat
testem/testem.out
)

INSTALL (FILES ${dai_TESTS_RESOURCES} DESTINATION share/dai/tests)
INSTALL (FILES ${dai_TESTS_TESTEM_RESOURCES} DESTINATION share/dai/tests/testem)

# unittests
ADD_SUBDIRECTORY (unit)
