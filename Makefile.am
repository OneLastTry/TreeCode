ACLOCAL_AMFLAGS = -I m4
AM_CPPFLAGS = $(BOOST_CPPFLAGS) -I$(srcdir)/TreeCode2

bin_PROGRAMS = simulator generator analyser
simulator_SOURCES = TreeCode2/Simulator.cpp
simulator_LDADD   = $(BOOST_PROGRAM_OPTIONS_LIBS)
simulator_LDFLAGS = $(BOOST_PROGRAM_OPTIONS_LDFLAGS)

generator_SOURCES = TreeCode2/ParticleGenerator.cpp
generator_LDADD   = $(BOOST_PROGRAM_OPTIONS_LIBS)
generator_LDFLAGS = $(BOOST_PROGRAM_OPTIONS_LDFLAGS)

analyser_SOURCES = TreeCode2/Analyser.cpp
analyser_LDADD   = $(BOOST_PROGRAM_OPTIONS_LIBS)
analyser_LDFLAGS = $(BOOST_PROGRAM_OPTIONS_LDFLAGS)

check_PROGRAMS = treecode_tests
treecode_tests_LDFLAGS = $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS) 
treecode_tests_LDADD   = $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)
treecode_tests_SOURCES = TreeCode2/tests/custom_asserts.cpp \
						 TreeCode2/tests/ParticleTest.cpp \
						 TreeCode2/tests/TreeTest.cpp \
						 TreeCode2/tests/IOTest.cpp \
						 TreeCode2/tests/DistributionTest.cpp
