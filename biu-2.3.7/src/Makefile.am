
################################################################
# the library header to distribute
################################################################
					
lib_biu_HEADERS = 	biu/OptionParser.hh \
					biu/VirtualList.hh \
					biu/HashMap.hh \
					biu/LimitedHash.hh \
					biu/LimitedHash.icc \
					biu/Matrix.hh \
					biu/Matrix.icc \
					biu/MatrixSparse.hh \
					biu/MatrixSparse.icc \
					biu/SquareMatrix.hh \
					biu/Timer.hh \
					biu/assertbiu.hh \
					biu/Alphabet.hh \
					biu/Alphabet.icc \
					biu/BioMolecule.hh \
					biu/BackboneStructure3D.hh \
					biu/Point.hh \
					biu/NeighborVector.hh \
					biu/LatticeDescriptor.hh \
					biu/LatticeDescriptor.icc \
					biu/LatticeDescriptorSQR.hh \
					biu/LatticeDescriptorCUB.hh \
					biu/LatticeDescriptorFCC.hh \
					biu/LatticeNeighborhood.hh \
					biu/LatticeNeighborhood.icc \
					biu/LatticeModel.hh \
					biu/LatticeModel.icc \
					biu/LatticeFrame.hh \
					biu/DistanceEnergyFunction.hh \
					biu/LatticeProtein.hh \
					biu/LatticeProtein_I.hh \
					biu/LatticeProtein_Ipnt.hh \
					biu/OffLatticeProtein.hh \
					biu/AllowedBasePairs.hh \
					biu/AllowedBasePairs.icc \
					biu/RNAStructure.hh \
					biu/RNAStructure_TB.hh \
					biu/RNAStructure_TB.icc \
					biu/Graph_UD.hh \
					biu/LatticeMoveSet.hh \
					biu/LatticeMoveSet.icc \
					biu/RandomNumberGenerator.hh \
					biu/RandomNumberFactory.hh \
					biu/QuasiRandomNumberGenerator.hh \
					biu/PullMoveSet.hh \
					biu/PivotMoveSet.hh

#					biu/ContactEnergyFunction.hh \
#
#

lib_biudir = $(includedir)/biu


lib_biu_util_HEADERS = 	biu/util/point.hh \
						biu/util/point.icc \
						biu/util/biu-string.h \
						biu/util/Util.h \
						biu/util/Util_String.h
						
lib_biu_utildir = $(includedir)/biu/util


################################################################
# headers and sources for quasi randum number generators
################################################################

QRNG_H = \
	biu/qrng/gsl_errno.h \
	biu/qrng/gsl_types.h \
	biu/qrng/gsl_qrng.h 

# header install directory via dummy lib
lib_biu_qrng_HEADERS = $(QRNG_H)
lib_biu_qrngdir =  $(includedir)/biu/qrng

QRNG_SOURCES = \
	biu/qrng/qrng.c \
	biu/qrng/niederreiter-2.c \
	biu/qrng/sobol.c


################################################################
# headers for Random123 random number generators
################################################################

Random123_H = \
	biu/Random123/array.h \
	biu/Random123/gsl_microrng.h \
	biu/Random123/threefry.h \
	biu/Random123/MicroURNG.hpp \
	biu/Random123/aes.h \
	biu/Random123/ars.h \
	biu/Random123/philox.h \
	biu/Random123/ReinterpretCtr.hpp
	
Random123_conventional_H = \
	biu/Random123/conventional/Engine.hpp \
	biu/Random123/conventional/gsl_cbrng.h
	
Random123_features_H = \
	biu/Random123/features/clangfeatures.h \
	biu/Random123/features/gccfeatures.h \
	biu/Random123/features/msvcfeatures.h \
	biu/Random123/features/open64features.h \
	biu/Random123/features/pgccfeatures.h \
	biu/Random123/features/sunprofeatures.h \
	biu/Random123/features/compilerfeatures.h \
	biu/Random123/features/iccfeatures.h \
	biu/Random123/features/nvccfeatures.h \
	biu/Random123/features/openclfeatures.h \
	biu/Random123/features/sse.h \
	biu/Random123/features/xlcfeatures.h


# header install directory via dummy lib
lib_biu_Random123_HEADERS = $(Random123_H)
lib_biu_Random123dir =  $(includedir)/biu/Random123
lib_biu_Random123_conventional_HEADERS = $(Random123_conventional_H)
lib_biu_Random123_conventionaldir =  $(includedir)/biu/Random123/conventional
lib_biu_Random123_features_HEADERS = $(Random123_features_H)
lib_biu_Random123_featuresdir =  $(includedir)/biu/Random123/features



################################################################
# the library objects
################################################################

lib_LIBRARIES = libbiu.a

libbiu_a_SOURCES =	$(lib_biu_HEADERS) \
					$(QRNG_H) \
					$(QRNG_SOURCES) \
					$(Random123_H) \
					biu/OptionParser.cc \
					biu/LatticeDescriptor.cc \
					biu/LatticeDescriptorSQR.cc \
					biu/LatticeDescriptorCUB.cc \
					biu/LatticeDescriptorFCC.cc \
					biu/LatticeModel.cc \
					biu/LatticeFrame.cc \
					biu/DistanceEnergyFunction.cc \
					biu/OffLatticeProtein.cc \
					biu/RNAStructure.cc \
					biu/RNAStructure_TB.cc \
					biu/LatticeProtein_I.cc \
					biu/LatticeProtein_Ipnt.cc \
					$(lib_biu_util_HEADERS) \
					biu/util/biu-string.cpp \
					biu/util/Util.cpp \
					biu/util/Util_String.cpp \
					biu/Graph_UD.cc \
					biu/RandomNumberGenerator.cc \
					biu/RandomNumberFactory.cc \
					biu/QuasiRandomNumberGenerator.cc \
					biu/PullMoveSet.cc \
					biu/PivotMoveSet.cc

#					biu/ContactEnergyFunction.cc \
#
#


################################################################
# further subdirs
################################################################
					
SUBDIRS = 
