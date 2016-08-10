

#include <iostream>
#include <biu/LatticeDescriptorCUB.hh>
#include <biu/LatticeDescriptorFCC.hh>
#include <biu/LatticeProtein_Ipnt.hh>
//#include <biu/PivotMoveSet.hh>

void testLatticeProtein_Ipnt(bool seqShared) {
	
	using namespace biu;
	
		// erzeugen	
	std::cout	<<"=> creating biu::LatticeDescriptorCUB() = ld" <<std::endl;
		LatticeDescriptorCUB ld ;
	std::cout	<<"=> creating biu::LatticeModel(ld) = lm" <<std::endl;
		LatticeModel lm(&ld);
	std::cout	<<"=> creating biu::Alphabet(HP,1) = alph" <<std::endl;
		Alphabet alph("HP",1);
	std::cout	<<"=> creating biu::EnergyMatrix(0.0,2,2) = E" <<std::endl;
		EnergyMatrix E = EnergyMatrix(	alph.getAlphabetSize(),
										alph.getAlphabetSize(),
										0.0 );
		E[alph.getIndex("H")][alph.getIndex("H")] = -1.0;
	std::cout	<<E <<std::endl;
	std::cout	<<"=> creating biu::ContactEnergyFunction( alph, E, lm) = f" 
				<<std::endl;
		ContactEnergyFunction f(&alph, &E, &lm);
		std::string seqStr = "PHPHHPH";
	std::cout	<<"=> "<<seqStr <<" = seqStr" <<std::endl;
		std::string relMoveStr = "FRUULU";
	std::cout	<<"=> "<<relMoveStr <<"  = relMoveStr" <<std::endl;
		std::string absMoveStr = lm.getString(
									lm.relMovesToAbsMoves(
										lm.getDescriptor()->getSequence(
											relMoveStr)));
	std::cout	<<"=> "<<absMoveStr<<"  = absMoveStr" <<std::endl;

	std::cout	<<"=> creating biu::LatticeProtein_Ipnt(lm,f,"
				<<"alph.getSequence(seqStr),seqShared,relMoveStr,false)"
				<<" = lp" 
				<<std::endl;
	biu::Sequence seq = alph.getSequence(seqStr);
		LatticeProtein_Ipnt lp(&lm,&f,&seq,seqShared,relMoveStr,false);

		// getLatticeModel()
	std::cout	<<"->(*(lp.getLatticeModel()) == lm)  = "
				<<(*(lp.getLatticeModel()) == lm?"true":"false")
				<<std::endl;
				
		// getEnergyFunction()
	std::cout	<<"->(*(getEnergyFunction()) == f)  = "
				<<(*(lp.getEnergyFunction()) == f?"true":"false")
				<<std::endl;
				
		// getSequence()
	std::cout	<<"->(alph.getString(lp.getSequence()) == seqStr)  = "
				<<(alph.getString(lp.getSequence()) == seqStr ? "true":"false")
				<<std::endl;
				
		// getStructure()
	std::cout	<<"->(lm.getString(lp.getStructure()) == relMoveStr)  = "
				<<(lm.getString(lp.getStructure()) == relMoveStr ?
				   "true":"false")
				<<std::endl;

		// getStringRepresentation
	std::cout	<<"->getStringRepresentation()  = "
				<<lp.getStringRepresentation()
				<<std::endl;

		// getEnergy()
	std::cout	<<"->getEnergy()  = "
				<<lp.getEnergy()
				<<std::endl;

		// getAbsoluteMoveString()
	std::cout	<<"->getMoveStrAbs()  = "
				<<lp.getMoveStrAbs()
				<<std::endl;
	std::cout	<<"->(lp.getMoveStrAbs() == absMoveStr)  = "
				<<(lp.getMoveStrAbs() == absMoveStr?"true":"false")
				<<std::endl;
				
		// getAbsoluteMoves()
	std::cout	<<"->lm.getString(lp.getMoveSeqAbs()) == absMoveStr  = "
				<<( lm.getString(lp.getMoveSeqAbs()) == absMoveStr ? 
						"true":"false")
				<<std::endl;
				
		// getRelativeMoveString()
	std::cout	<<"->getMoveStrRel()  = "
				<<lp.getMoveStrRel()
				<<std::endl;
	std::cout	<<"->(lm.getString(lp.getMoveSeqRel())"
				<<" == lp.getMoveStrRel()  = "
				<<(lm.getString(lp.getMoveSeqRel()) 
					== lp.getMoveStrRel() ? "true":"false")
				<<std::endl;
				
		// getRelativeMoves()
	std::cout <<"->(lm.getString(lm.relMovesToAbsMoves(lp.getMoveSeqRel()))"
				<<" == absMoveStr)  = "
				<<(lm.getString(lm.relMovesToAbsMoves(lp.getMoveSeqRel())) ==
					absMoveStr ? "true":"false")
				<<std::endl;
				
		// get3Ddata()
	std::cout	<<"->get3Ddata().size()  = "
				<<lp.get3Ddata().size()
				<<std::endl;
				
		// getDRMSD()
	std::cout	<<"->getDRMSD(*this)  = "
				<<((int)(lp.getDRMSD(lp)*1000))
				<<std::endl;
				
		// isSelfAvoiding()
	std::cout	<<"->isSelfAvoiding()  = "
				<<(lp.isSelfAvoiding()?"true":"false")
				<<std::endl;

		// isConnected()
	std::cout	<<"->isConnected()  = "
				<<(lp.isConnected()?"true":"false")
				<<std::endl;

		// isValid()
	std::cout	<<"->isValid()  = "
				<<(lp.isValid()?"true":"false")
				<<std::endl;
				
		// assignement
	std::cout	<<"->lp2 = lp"
				<<std::endl;
	LatticeProtein_Ipnt lp2 = lp;
	
		// comparison
	std::cout	<<"-> lp2 == lp ? "
				<<( lp2 == lp ? "true" : "false")
				<<std::endl;
	std::cout	<<"-> lp2 != lp ? "
				<<( lp2 != lp ? "true" : "false")
				<<std::endl;
	
		// setMoveStrAbs
	std::cout	<< "lp2.getStringRepresentation() = " 
				<< lp2.getStringRepresentation() << std::endl;
	std::string moveString = "FFFFFF";
	std::cout	<< "lp2.setMoveStrAbs(moveString)" << std::endl;
	lp2.setMoveStrAbs(moveString);
	std::cout	<< "lp2.getStringRepresentation() = " 
					<< lp2.getStringRepresentation() << std::endl;
	
	// test getLenght
	std::cout	<< "lp2.getLength() = " << lp2.getLength() << std::endl;
	
	
	// test operator= (LatticeProtein_Ipnt)
	std::cout	<<"LatticeProtein_Ipnt lp_ipnt(&lm,&f,&seq,seqShared,\"FFLFFF\",false);"
				<< std::endl;
	LatticeProtein_Ipnt lp_ipnt(&lm,&f,&seq,seqShared,"FFLFFF",false);
	std::cout	<<"lp_ipnt = lp;" << std::endl;
	lp_ipnt = lp;
	std::cout	<<"-> lp_ipnt == lp ? "
				<<( lp_ipnt == lp ? "true" : "false")
				<<std::endl;
	std::cout	<<"-> lp_ipnt != lp ? "
				<<( lp_ipnt != lp ? "true" : "false")
				<<std::endl;
	
	// test operator= (LatticeProtein_I)
	std::cout	<<"LatticeProtein_Ipnt lp_ipnt2(&lm,&f,&seq,seqShared,relMoveStr,false);"
				<< std::endl;
	LatticeProtein_Ipnt lp_ipnt2(&lm,&f,&seq,seqShared,relMoveStr,false);
	std::cout	<<"lp_ipnt2 = *(dynamic_cast<LatticeProtein_I*>(&lp2));"
				<< std::endl;
	lp_ipnt2 = *(dynamic_cast<LatticeProtein_I*>(&lp2));
	std::cout	<<"-> lp_ipnt2 == lp2 ? "
				<<( lp_ipnt2 == lp2 ? "true" : "false")
				<<std::endl;
	std::cout	<<"-> lp_ipnt2 != lp2 ? "
				<<( lp_ipnt2 != lp2 ? "true" : "false")
				<<std::endl;
	
	// test fromString
	std::cout	<< "LatticeProtein_Ipnt* lp_fromString = "
					"lp.fromString(\"RRRRRR(PPPPPPP)\")"
				<< std::endl;
	std::string stringRep = "RRRRRR(PPPPPPP)";
	LatticeProtein_I* lp_fromString = lp.fromString(stringRep);
	std::cout 	<< "lp_fromString->getStringRepresentation() = "
				<< lp_fromString->getStringRepresentation()
				<< std::endl;

}

void testLatticeProtein_Ipnt_FCC(bool seqShared) {
	
	using namespace biu;
	
		// erzeugen	
	std::cout	<<"=> creating biu::LatticeDescriptorFCC() = ld" <<std::endl;
		LatticeDescriptorFCC ld ;
	std::cout	<<"=> creating biu::LatticeModel(ld) = lm" <<std::endl;
		LatticeModel lm(&ld);
	std::cout	<<"=> creating biu::Alphabet(HP,1) = alph" <<std::endl;
		Alphabet alph("HP",1);
	std::cout	<<"=> creating biu::EnergyMatrix(0.0,2,2) = E" <<std::endl;
		EnergyMatrix E = EnergyMatrix(	alph.getAlphabetSize(),
										alph.getAlphabetSize(),
										0.0 );
		E[alph.getIndex("H")][alph.getIndex("H")] = -1.0;
	std::cout	<<E <<std::endl;
	std::cout	<<"=> creating biu::ContactEnergyFunction( alph, E, lm) = f" 
				<<std::endl;
		ContactEnergyFunction f(&alph, &E, &lm);
		std::string seqStr = "PHPHHPH";
	std::cout	<<"=> "<<seqStr <<" = seqStr" <<std::endl;
	std::string absMoveStr = "LULULULULULU";
	std::cout	<<"=> "<<absMoveStr <<"  = absMoveStr" <<std::endl;
	
	std::cout	<<"=> creating biu::LatticeProtein_Ipnt(lm,f,"
				<<"alph.getSequence(seqStr),seqShared,relMoveStr,true)"
				<<" = lp" 
				<<std::endl;
	biu::Sequence seq = alph.getSequence(seqStr);
		LatticeProtein_Ipnt lp(&lm,&f,&seq,seqShared,absMoveStr,true);

		// getLatticeModel()
	std::cout	<<"->(*(lp.getLatticeModel()) == lm)  = "
				<<(*(lp.getLatticeModel()) == lm?"true":"false")
				<<std::endl;
				
		// getEnergyFunction()
	std::cout	<<"->(*(getEnergyFunction()) == f)  = "
				<<(*(lp.getEnergyFunction()) == f?"true":"false")
				<<std::endl;
				
		// getSequence()
	std::cout	<<"->(alph.getString(lp.getSequence()) == seqStr)  = "
				<<(alph.getString(lp.getSequence()) == seqStr ? "true":"false")
				<<std::endl;
	
	std::cout<<"->printing coordinates:" << std::endl;
	DPointVec vec = lp.get3Ddata();
	for (DPointVec::const_iterator it = vec.begin();  
	it != vec.end(); it++) 
	{
		std::cout << "  "<<  *it << std::endl;
	}	

		// getStructure()
	std::cout	<<"->(lm.getString(lp.getStructure()) == absMoveStr)  = "
				<<(lm.getString(lp.getStructure()) == absMoveStr ?
				   "true":"false")
				<<std::endl;

		// getStringRepresentation
	std::cout	<<"->getStringRepresentation()  = "
				<<lp.getStringRepresentation()
				<<std::endl;

		// getEnergy()
	std::cout	<<"->getEnergy()  = "
				<<lp.getEnergy()
				<<std::endl;

		// getAbsoluteMoveString()
	std::cout	<<"->getMoveStrAbs()  = "
				<<lp.getMoveStrAbs()
				<<std::endl;
	std::cout	<<"->(lp.getMoveStrAbs() == absMoveStr)  = "
				<<(lp.getMoveStrAbs() == absMoveStr?"true":"false")
				<<std::endl;
				
		// getAbsoluteMoves()
	std::cout	<<"->lm.getString(lp.getMoveSeqAbs()) == absMoveStr  = "
				<<( lm.getString(lp.getMoveSeqAbs()) == absMoveStr ? 
						"true":"false")
				<<std::endl;
				
		// getRelativeMoveString()
	std::cout	<<"->getMoveStrRel()  = "
				<<lp.getMoveStrRel()
				<<std::endl;
	std::cout	<<"->(lm.getString(lp.getMoveSeqRel())"
				<<" == lp.getMoveStrRel()  = "
				<<(lm.getString(lp.getMoveSeqRel()) 
					== lp.getMoveStrRel() ? "true":"false")
				<<std::endl;
			
	
	std::cout 	<< "lm.getString(lp.getMoveSeqAbs()) = " 
				<< lm.getString(lp.getMoveSeqAbs()) << std::endl;
	std::cout 	<< "lm.getString(lp.getMoveSeqRel()) = " 
				<< lm.getString(lp.getMoveSeqRel()) << std::endl;
	
		// getRelativeMoves()
	std::cout << "lm.getString(lm.relMovesToAbsMoves(lp.getMoveSeqRel()) = "
			<< lm.getString(lm.relMovesToAbsMoves(lp.getMoveSeqRel()))
			<< std::endl;
	std::cout << "absMoveStr = " << absMoveStr << std::endl;
	std::cout <<"->(lm.getString(lm.relMovesToAbsMoves(lp.getMoveSeqRel()))"
				<<" == absMoveStr)  = "
				<<(lm.getString(lm.relMovesToAbsMoves(lp.getMoveSeqRel())) ==
					absMoveStr ? "true":"false")
				<<std::endl;
				
		// get3Ddata()
	std::cout	<<"->get3Ddata().size()  = "
				<<lp.get3Ddata().size()
				<<std::endl;
				
		// getDRMSD()
	std::cout	<<"->getDRMSD(*this)  = "
				<<((int)(lp.getDRMSD(lp)*1000))
				<<std::endl;
				
		// isSelfAvoiding()
	std::cout	<<"->isSelfAvoiding()  = "
				<<(lp.isSelfAvoiding()?"true":"false")
				<<std::endl;

		// isConnected()
	std::cout	<<"->isConnected()  = "
				<<(lp.isConnected()?"true":"false")
				<<std::endl;

		// isValid()
	std::cout	<<"->isValid()  = "
				<<(lp.isValid()?"true":"false")
				<<std::endl;
				
		// assignement
	std::cout	<<"->lp2 = lp"
				<<std::endl;
	LatticeProtein_Ipnt lp2 = lp;
	
		// comparison
	std::cout	<<"-> lp2 == lp ? "
				<<( lp2 == lp ? "true" : "false")
				<<std::endl;
	std::cout	<<"-> lp2 != lp ? "
				<<( lp2 != lp ? "true" : "false")
				<<std::endl;
	
		// setMoveStrAbs
	std::cout	<< "lp2.getStringRepresentation() = " 
				<< lp2.getStringRepresentation() << std::endl;
	std::string moveString = "FLFLFLFLFLFL";
	std::cout	<< "lp2.setMoveStrAbs(moveString)" << std::endl;
	lp2.setMoveStrAbs(moveString);
	std::cout	<< "lp2.getStringRepresentation() = " 
					<< lp2.getStringRepresentation() << std::endl;
	
	// test getLenght
	std::cout	<< "lp2.getLength() = " << lp2.getLength() << std::endl;
	
	
	// test operator= (LatticeProtein_Ipnt)
	std::cout	<<"LatticeProtein_Ipnt lp_ipnt(&lm,&f,&seq,seqShared,\"FLFLFLFLFLFL\",false);"
				<< std::endl;
	LatticeProtein_Ipnt lp_ipnt(&lm,&f,&seq,seqShared,"FLFLFLFLFLFL",false);
	std::cout	<<"lp_ipnt = lp;" << std::endl;
	lp_ipnt = lp;
	std::cout	<<"-> lp_ipnt == lp ? "
				<<( lp_ipnt == lp ? "true" : "false")
				<<std::endl;
	std::cout	<<"-> lp_ipnt != lp ? "
				<<( lp_ipnt != lp ? "true" : "false")
				<<std::endl;
	
	// test operator= (LatticeProtein_I)
	std::cout	<<"LatticeProtein_Ipnt lp_ipnt2(&lm,&f,&seq,seqShared,relMoveStr,false);"
				<< std::endl;
	LatticeProtein_Ipnt lp_ipnt2(&lm,&f,&seq,seqShared,absMoveStr,true);
	std::cout	<<"lp_ipnt2 = *(dynamic_cast<LatticeProtein_I*>(&lp2));"
				<< std::endl;
	lp_ipnt2 = *(dynamic_cast<LatticeProtein_I*>(&lp2));
	std::cout	<<"-> lp_ipnt2 == lp2 ? "
				<<( lp_ipnt2 == lp2 ? "true" : "false")
				<<std::endl;
	std::cout	<<"-> lp_ipnt2 != lp2 ? "
				<<( lp_ipnt2 != lp2 ? "true" : "false")
				<<std::endl;

	// test fromString
	std::cout	<< "LatticeProtein_Ipnt* lp_fromString = "
					"lp.fromString(\"RRRRRR(PPPPPPP)\")"
				<< std::endl;
	std::string stringRep = "FRFRFRFRFRFR(PPPPPPP)";
	LatticeProtein_I* lp_fromString = lp.fromString(stringRep);
	std::cout 	<< "lp_fromString->getStringRepresentation() = "
				<< lp_fromString->getStringRepresentation()
				<< std::endl;
}

void testFCCEnergy() {
	using namespace biu;
	
		// erzeugen	
	std::cout	<<"=> creating biu::LatticeDescriptorFCC() = ld" <<std::endl;
		LatticeDescriptorFCC ld ;
	std::cout	<<"=> creating biu::LatticeModel(ld) = lm" <<std::endl;
		LatticeModel lm(&ld);
	std::cout	<<"=> creating biu::Alphabet(HP,1) = alph" <<std::endl;
		Alphabet alph("HP",1);
	std::cout	<<"=> creating biu::EnergyMatrix(0.0,2,2) = E" <<std::endl;
		EnergyMatrix E = EnergyMatrix(	alph.getAlphabetSize(),
										alph.getAlphabetSize(),
										0.0 );
		E[alph.getIndex("H")][alph.getIndex("H")] = -1.0;
	std::cout	<<E <<std::endl;
	std::cout	<<"=> creating biu::ContactEnergyFunction( alph, E, lm) = f" 
				<<std::endl;
		ContactEnergyFunction f(&alph, &E, &lm);
		std::string seqStr = "HHHH";
	std::cout	<<"=> "<<seqStr <<" = seqStr" <<std::endl;
	std::string absMoveStr = "FLFRBR";
	std::cout	<<"=> "<<absMoveStr <<"  = absMoveStr" <<std::endl;
	
	std::cout	<<"=> creating biu::LatticeProtein_Ipnt(lm,f,"
				<<"alph.getSequence(seqStr),seqShared,absoveStr,true)"
				<<" = lp" 
				<<std::endl;
	bool seqShared = false;
	biu::Sequence seq = alph.getSequence(seqStr);
	LatticeProtein_Ipnt lp(&lm,&f,&seq,seqShared,absMoveStr,true);
	std::cout << "lp.getEnergy() = " << lp.getEnergy() << std::endl;
}

int main(int argc, char **argv) {


	std::cout	<<"========== BEGIN OF TEST =================\n"
	<<" test class : LatticeProtein_Ipnt\n"
	<<"==========================================" <<std::endl;
	
	std::cout	<<"seqShared = false" << std::endl;
	testLatticeProtein_Ipnt(false);
	std::cout	<< std::endl;
	std::cout	<<"seqShared = true" << std::endl;
	testLatticeProtein_Ipnt(true);
	
	std::cout	<< std::endl << "using FCC:" << std::endl;
	std::cout	<<"seqShared = false" << std::endl;
	testLatticeProtein_Ipnt_FCC(false);
	std::cout	<< std::endl;
	std::cout	<<"seqShared = true" << std::endl;
	testLatticeProtein_Ipnt_FCC(true);
	std::cout	<<"testing getEnergy with FCC:" << std::endl;
	testFCCEnergy();
	
	std::cout	<<"============= END OF TEST ================\n" <<std::endl;
	
	return 0;
}
