#include "biu/PivotMoveSet.hh"
#include "biu/LatticeDescriptorCUB.hh"
#include "biu/LatticeDescriptorSQR.hh"
#include "biu/DistanceEnergyFunction.hh"

using namespace biu;
	
void doFolding( std::string seqStr, std::string absMoveStr, LatticeDescriptor * ld) {
	LatticeModel lm(ld);
	Alphabet alph("HP",1);
	EnergyMatrix E = EnergyMatrix(	alph.getAlphabetSize(),
									alph.getAlphabetSize(),
									0.0 );
	E[alph.getIndex("H")][alph.getIndex("H")] = -1.0;
	ContactEnergyFunction f(&alph, &E, &lm);
	PivotMoveSet pms(&lm);
	
		/*
	 *		Test 1:
	 */ 
	std::cout	<<"=> "<<seqStr <<" = seqStr" <<std::endl;
	bool seqShared = false;
	bool isAbsoluteMove = true;
	std::cout	<<"=> "<<absMoveStr <<"  = absMoveStr" <<std::endl;
	std::cout	<<"=> "<<seqShared <<" = seqShared" <<std::endl;
	
	std::cout	<<"=> creating biu::LatticeProtein_Ipnt(lm,f,seqStr,seqShared"
				<<",absMoveStr,isAbsoluteMove)"
				<<" = lp" 
				<<std::endl;
	biu::Sequence seq = alph.getSequence(seqStr);
	LatticeProtein_Ipnt* lp = new LatticeProtein_Ipnt
							(&lm,&f,&seq,seqShared,absMoveStr,isAbsoluteMove);
	std::cout	<< "->getMoveStrRel = " 
				<< lp->getMoveStrRel() << std::endl;
	std::cout	<< "pms.getMoveNumber(lp) = "
				<< pms.getMoveNumber(lp) << std::endl;
	
	// applying moves
	for (size_t i=0; i<pms.getMoveNumber(lp); i++)
	{
		std::cout << "Move: " << i << std::endl;
		LatticeProtein* neigh = lp->clone();
		neigh = pms.applyMoveInPlace(neigh,i);
		if (neigh!=NULL) {
			std::cout 	<< "pms.applyMoveInPlace(lp," << i << ") = "
			<< neigh->getMoveStrRel()
			<< " " << (neigh->isSelfAvoiding() ? "selfavoiding" : "NOT selfavoiding!")
			<< std::endl;
			std::cout 	<< "pms.undoLastMove(lp)->getMoveStrRel() = "
			<< pms.undoLastMove(neigh)->getMoveStrRel()
			<< std::endl;
			delete neigh;
		}
	}

	delete lp;
}

int main() {
	std::cout	<<"========== BEGIN OF TEST ==================\n"
	<<" test function : PivotMoveSet\n"
	<<"===========================================" <<std::endl;

	LatticeDescriptorSQR ld_SQR;
	LatticeDescriptorCUB ld_CUB;

	doFolding("HHHH","FFF",&ld_CUB);
	doFolding("HPHPPHP","FRRBLB",&ld_SQR);
	
	std::cout	<<"============= END OF TEST =================\n" <<std::endl;
}
