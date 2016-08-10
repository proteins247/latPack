
TESTINCLUDES = \
	LatticeProteinUtil.hh

TESTSOURCES = \
	testUtilString.cc \
	testOptionParser.cc \
	testLimitedHash.cc \
	testRNG_LCG.cc \
	testQuasiRandomNumberGenerator.cc \
	testRNG_ARS4x32.cc \
	testVirtualList.cc \
	testMatrix.cc \
	testMatrixSparse.cc \
	testSquareMatrix.cc \
	testAlphabet.cc \
	testAllowedBasePairs.cc \
	testRNAStructure_TB.cc \
	testPoint.cc \
	testNeighborVector.cc \
	testLatticeNeighborhood.cc \
	testLatticeDescriptors.cc \
	testLatticeDescriptorSQR.cc \
	testLatticeDescriptorCUB.cc \
	testLatticeModel.cc \
	testLatticeFrame.cc \
	testUtilLattice.cc \
	testPivotMoveSet.cc \
	testPullMoveSet.cc \
	testContactEnergyFunction.cc \
	testIntervalEnergyFunction.cc \
	testLatticeProtein_Ipnt.cc \
	testLatticeProteinUtil.cc \
	testGraph_UD.cc 


EXTRA_DIST = $(TESTSOURCES) $(TESTINCLUDES) $(TESTSOURCES:.cc=.verified)

# set test specific compiler and linker options
TESTCXXFLAGS = -I$(top_srcdir)/src 
TESTLDFLAGS  = 
TESTLIBS     = $(top_builddir)/src/libbiu.a


#################################################################
# AUXILIARY VARIABLES FOR TEST BUILD AND EXECUTION
#################################################################

# final compiler and linker flag adds
AM_CXXFLAGS = $(TESTCXXFLAGS) @AM_CXXFLAGS@ -g -O0 -DDEBUG=1
AM_LDFLAGS = $(TESTLDFLAGS) @AM_LDFLAGS@
LIBS = $(TESTLIBS) @LIBS@

DIFFS=$(TESTSOURCES:.cc=.diff)

TESTCXXCOMPILE = $(CXX) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) \
	$(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CXXFLAGS) $(CXXFLAGS)
TESTCXXLINK = $(CXX) $(AM_CXXFLAGS) $(CXXFLAGS) $(AM_LDFLAGS) $(LDFLAGS) \
	-o $@


#################################################################
# MAKE TARGETS FOR TEST BUILD AND EXECUTION
#################################################################

all:

tests: $(DIFFS)


.cc.diff:
	@name=$*; name=$${name/test/}; \
	$(PRINTF) " Building test %-38s" $$name;
	    @$(TESTCXXCOMPILE) -c -o $*.o $<
	    @$(TESTCXXLINK) $*.o $(LIBS) -o $*.test
	@$(PRINTF) " --> execution :";
	    @./$*.test > $*.lastout
	@name=$*; name=$${name/test/}; \
	if [ -e $*.verified ] ; then \
	    if $(DIFF) -w $*.verified $*.lastout >$@ ; then \
	        $(PRINTF) " OK\n"; \
	    else \
	        $(PRINTF) " FAILED\n";\
	        $(PRINTF) "------------------------------------------------------\n" ;\
	        $(PRINTF) " Diff for $$name (verified '<', current '>')\n" ;\
	        $(PRINTF) "------------------------------------------------------\n" ;\
	        $(CAT) $@ ;\
	        $(PRINTF) "------------------------------------------------------\n" ;\
	        $(PRINTF) " [ If new output is correct: mv $*.lastout $*.verified ]\n" ;\
	        $(PRINTF) "------------------------------------------------------\n" ;\
	        $(PRINTF) "\n" ;\
	    fi ;\
	    $(RM) $@ ;\
	else \
	    $(PRINTF) " DONE\n"; \
	    $(PRINTF) "------------------------------------------------------\n" ;\
	    $(PRINTF) " Output of $$name :\n" ;\
	    $(PRINTF) "------------------------------------------------------\n" ;\
	    $(CAT) $*.lastout ;\
	    $(PRINTF) "------------------------------------------------------\n" ;\
	    $(PRINTF) " WARNING: verified output of $$name is missing. Please generate\n" ;\
	    $(PRINTF) " a file $*.verified that contains the verfied output of $$name.\n"; \
	    $(PRINTF) " [ If output is correct: mv $*.lastout $*.verified ]\n" ;\
	    $(PRINTF) "------------------------------------------------------\n" ;\
	    $(PRINTF) "\n" ;\
	fi

clean-local:
	$(RM) -rf $(TESTSOURCES:.cc=.test) $(TESTSOURCES:.cc=.o) \
	    $(TESTSOURCES:.cc=.lastout) $(TESTSOURCES:.cc=.diff) \
	    $(TESTSOURCES:.cc=.exe)

#################################################################
