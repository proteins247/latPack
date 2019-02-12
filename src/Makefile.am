
####################################################################

# the program to build
if with_latfold
bin_PROGRAMS = test latFit latVec latFoldVec latConv latMap latMapTraj latFold latNum # latFit-old
else
bin_PROGRAMS = latFit latVec latConv latMap latMapTraj latNum
endif

AM_CPPFLAGS = -I$(top_srcdir)/src

#######################   LatFit PROGRAM   #########################

latFit_SOURCES =	\
					biu/LatticeDescriptorCKW.hh \
					biu/LatticeDescriptorCKW.cc \
					biu/Rotator3D.hh \
					biu/Rotator3D.cc \
					biu/LatticeProteinUtil.hh \
					biu/LatticeProteinUtil.cc \
					biu/SuperPos_Kabsch.hh \
					biu/SuperPos_Kabsch.cc \
					latFit.cc

#latFit_old_SOURCES =	\
#					biu/LatticeDescriptorCKW.hh \
#					biu/LatticeDescriptorCKW.cc \
#					biu/Rotator3D.hh \
#					biu/Rotator3D.cc \
#					biu/LatticeProteinUtil.hh \
#					biu/LatticeProteinUtil.cc \
#					biu/SuperPos_Kabsch.hh \
#					biu/SuperPos_Kabsch.cc \
#					latFit-old.cc
#
test_SOURCES =	\
				kabsch.hh \
				kabsch.cc \
				biu/LatticeProteinUtil.hh \
				biu/LatticeProteinUtil.cc \
				biu/SuperPos_Kabsch.hh \
				biu/SuperPos_Kabsch.cc \
				testSuperPos.cc


#######################   LatFit PROGRAM   #########################

latConv_SOURCES =	\
					Exception.hh \
					Exception.cc \
					PDBsupport.hh \
					latConv.hh \
					latConv.cc


#######################   LatVec PROGRAM   #########################

latVec_SOURCES =	\
					PDBsupport.hh \
					energyFileSupport.hh \
					latVec.cc


#######################   LatMap PROGRAM   #########################

latNum_SOURCES =	\
					biu/LatticeDescriptorCKW.hh \
					biu/LatticeDescriptorCKW.cc \
					biu/LatticeProteinUtil.hh \
					biu/LatticeProteinUtil.cc \
					latNum.cc

#######################   LatMap PROGRAM   #########################

latMap_SOURCES =	\
					biu/LatticeDescriptorCKW.hh \
					biu/LatticeDescriptorCKW.cc \
					biu/LatticeProteinUtil.hh \
					biu/LatticeProteinUtil.cc \
					biu/SuperPos_Kabsch.hh \
					biu/SuperPos_Kabsch.cc \
					latMap.cc

#######################   LatMap PROGRAM   #########################

latMapTraj_SOURCES =	\
					biu/LatticeDescriptorCKW.hh \
					biu/LatticeDescriptorCKW.cc \
					biu/LatticeProteinUtil.hh \
					biu/LatticeProteinUtil.cc \
					biu/SuperPos_Kabsch.hh \
					biu/SuperPos_Kabsch.cc \
					HDF5_support.hh \
					HDF5_support.cc \
					latMapTraj.cc

#######################   LatFold PROGRAM   ########################

if with_latfold
latFold_SOURCES =	\
					SC_CountTarget.hh \
					SC_CountTarget.cc \
					SC_MinE.hh \
					SC_MinE.cc \
					SC_OutAbs.hh \
					SC_OutAbs.cc \
					SC_OutEnergy.hh \
					SC_OutEnergy.cc \
					energyFileSupport.hh \
					HDF5_support.hh \
					HDF5_support.cc \
					latFold.cc 
endif

#######################   LatFoldVec PROGRAM   ########################

if with_latfold
latFoldVec_SOURCES =	\
					SC_CountTarget.hh \
					SC_CountTarget.cc \
					SC_MinE.hh \
					SC_MinE.cc \
					SC_OutAbs.hh \
					SC_OutAbs.cc \
					SC_OutEnergy.hh \
					SC_OutEnergy.cc \
					energyFileSupport.hh \
					HDF5_support.hh \
					HDF5_support.cc \
					latFoldVec.cc 
endif

####################################################################

