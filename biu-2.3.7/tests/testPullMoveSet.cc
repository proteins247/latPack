#include "biu/PullMoveSet.hh"
#include "biu/LatticeDescriptorSQR.hh"
#include "biu/LatticeDescriptorCUB.hh"
#include "biu/LatticeDescriptorFCC.hh"
#include "biu/LatticeModel.hh"
#include "biu/LatticeProtein.hh"
#include "biu/LatticeProtein_Ipnt.hh"
#include "biu/DistanceEnergyFunction.hh"
#include <iostream>

#include "LatticeProteinUtil.hh"

using namespace biu;

void print3DData(LatticeProtein* lp)
{
	if (lp!=NULL) {
		std::cout<<"->printing coordinates:" << std::endl;
		DPointVec vec = lp->get3Ddata();
		for (DPointVec::const_iterator it = vec.begin();
		it != vec.end(); it++)
		{
			std::cout << "  "<<  *it << std::endl;
		}
	}
}

void testPullMoveDecoderNeighborhood(LatticeDescriptor * ld)
{
	std::cout	<<"========== BEGIN OF TEST ==================\n"
				<<" test function : PullMoveSet::PullMoveDecoder\n"
				<<"===========================================" <<std::endl;

	LatticeModel * lattice = new LatticeModel(ld);
	std::cout 	<< "=> creating PullMoveSet::PullMoveDecoder(new Lattice "
				<< ld->getName() <<")"
				<< std::endl;
	PullMoveSet::PullMoveDecoder * decoder = new PullMoveSet::PullMoveDecoder(
			lattice);

	const LatticeNeighborhood neighbors = lattice->getNeighborhood();
	
	std::cout << " # get all possible pull moves for each neighboring vector :" << std::endl;
	// for every possivle direction vector
	for (LatticeNeighborhood::const_iterator it = neighbors.begin();
			it != neighbors.end(); it++) {
		std::cout << "direction vector: " << (*it) << std::endl;
		// get every set of C and L
		for (size_t i=0; i<decoder->getPullMoveNumber(); i++) {
			PullMoveSet::IPointPair clPair = decoder->lookupStdMove(IntPoint(0,0,0), (*it), i);
			std::cout 	<< "  Pull-Move: " << i << std::endl;
			std::cout 	<< "    C: " << clPair.first << std::endl;
			std::cout 	<< "    L: " << clPair.second << std::endl;
			// check validity
			std::cout 	<< "    C and L neighbored: "
						<< (lattice->areNeighbored(clPair.first, clPair.second) ? "true" : "false")
						<< std::endl;
			std::cout 	<< "    C is neighbored to (0,0,0)? "
						<< (lattice->areNeighbored(IntPoint(0,0,0), clPair.first) ? "true" : "false")
						<< std::endl;
			std::cout 	<< "    L is neighbored to fixed direction? "
						<< (lattice->areNeighbored((*it), clPair.second) ? "true" : "false")
						<< std::endl;
		}
	}
	
	std::cout << " # get all possible end moves for (0,0,0) :" << std::endl;
	for (size_t i=0; i<decoder->getEndMoveNumber(); i++) {
		PullMoveSet::IPointPair clPair = decoder->lookupEndMove(IntPoint(0,0,0), i);
		std::cout 	<< "  End-Move: " << i << std::endl;
		std::cout 	<< "    1: " << clPair.first << std::endl;
		std::cout 	<< "    2: " << clPair.second << std::endl;
		// check validity
		std::cout 	<< "    1 is neighbored to (0,0,0)? "
					<< (lattice->areNeighbored(IntPoint(0,0,0), clPair.first) ? "true" : "false")
					<< std::endl;
		std::cout 	<< "    1 and 2 neighbored: "
					<< (lattice->areNeighbored(clPair.first, clPair.second) ? "true" : "false")
					<< std::endl;
	}

	delete decoder;
	delete lattice;

	std::cout	<<"============= END OF TEST =================\n" <<std::endl;
}

void testPullMoveDecoderLookup(LatticeDescriptor* ld, std::string seqStr, std::string absMoveStr)
{
	std::cout	<<"========== BEGIN OF TEST ==================\n"
				<<" test PullMoveDecoder lookup\n"
				<<"===========================================" <<std::endl;
	std::cout << "=> creating helper objects" << std::endl;
	const LatticeModel lm(ld);
	const Alphabet alph("HP",1);
	EnergyMatrix E = EnergyMatrix(	alph.getAlphabetSize(),
									alph.getAlphabetSize(),
									0.0 );
	E[alph.getIndex("H")][alph.getIndex("H")] = -1.0;
	const ContactEnergyFunction f(&alph, &E, &lm);
	bool seqShared = false;

	std::cout 	<< "=> creating PullMoveSet(new LatticeModel(new LatticeDescriptor "<<ld->getName()<<" ))"
				<< std::endl;
	PullMoveSet pms(&lm);

	std::cout	<<"=> "<<seqStr <<" = seqStr" <<std::endl;
	const biu::Sequence seq = alph.getSequence(seqStr);
	std::cout	<<"=> "<<absMoveStr <<"  = absMoveStr" <<std::endl;
	std::cout	<<"=> creating biu::LatticeProtein_Ipnt(lm,f,seq,absMoveStr,true)"
	<<" = lp"
	<<std::endl;
	LatticeProtein_Ipnt* lp = new LatticeProtein_Ipnt
								(&lm,&f,&seq,seqShared,absMoveStr,true);
	std::cout	<< "->getMoveStrAbs = "
				<< lp->getMoveStrAbs() << std::endl;
	std::cout	<< "->getMoveStrRel = "
				<< lp->getMoveStrRel() << std::endl;
	std::cout	<< "pms.getMoveNumber(lp) = "
				<< pms.getMoveNumber(lp) << std::endl;

	bool pullFront;
	bool stdPull;
	size_t moveIndex;
	size_t movePosition;
	size_t size = lp->getSequence().size();

	for (size_t i=0; i<pms.getMoveNumber(lp); i++)
	{
		moveIndex = i;
		std::cout 	<< "pms->decoder.lookupMove(lp.getSequence().size(), "
		<< i
		<< ",pullFront,stdPull,movePosition)"
		<< std::endl;
		pms.getDecoder()->lookupMove(size, moveIndex,pullFront,stdPull,movePosition);

		std::cout	<< "  pullFront = " << (pullFront ? "true" : "false")
		<< std::endl;
		std::cout	<< "  stdPull = " << (stdPull ? "true" : "false")
		<< std::endl;
		if (!stdPull)
		{
			PullMoveSet::IPointPair pullMoves = pms.getDecoder()->lookupEndMove(IntPoint(0,0,0), moveIndex);
			std::cout 	<< "  endPullPos0 = "
						<< pullMoves.first
						<< std::endl;
			std::cout 	<< "  endPullPos1 = "
						<< pullMoves.second
						<< std::endl;

		}
		std::cout	<< "  moveIndex = " << moveIndex
		<< std::endl;
		std::cout	<< "  movePosition = " << movePosition
		<< std::endl;
	}

	delete lp;

  	std::cout	<<"============= END OF TEST =================\n" <<std::endl;
}

void testPullMoveSet(LatticeDescriptor* ld, std::string seqStr, std::string relMoveStr, bool isAbsoluteMove )
{
	std::cout	<<"========== BEGIN OF TEST ==================\n"
	<<" test function : PullMoveSet\n"
	<<"===========================================" <<std::endl;

	std::cout << "=> creating helper objects" << std::endl;
	LatticeModel lm(ld);
	Alphabet alph("HP",1);
	EnergyMatrix E = EnergyMatrix(	alph.getAlphabetSize(),
									alph.getAlphabetSize(),
									0.0 );
	E[alph.getIndex("H")][alph.getIndex("H")] = -1.0;
	ContactEnergyFunction f(&alph, &E, &lm);

	std::cout 	<< "=> creating PullMoveSet(new LatticeModel"
				<< "(new LatticeDescriptor "<<ld->getName()<<"))"
				<< std::endl;
	PullMoveSet pms(&lm);

	/*
	 *		Test 1:
	 */
	std::cout	<<"=> "<<seqStr <<" = seqStr" <<std::endl;
	bool seqShared = false;
	std::cout	<<"=> "<<relMoveStr <<"  = relMoveStr" <<std::endl;
	std::cout	<<"=> "<<seqShared <<" = seqShared" <<std::endl;

	std::cout	<<"=> creating biu::LatticeProtein_Ipnt(lm,f,seqStr,seqShared"
				<<",relMoveStr,isAbsoluteMove)"
				<<" = lp"
				<<std::endl;
	biu::Sequence seq = alph.getSequence(seqStr);
	LatticeProtein_Ipnt* lp = new LatticeProtein_Ipnt
							(&lm,&f,&seq,seqShared,relMoveStr,isAbsoluteMove);
	std::cout	<< "->getMoveStrRel = "
				<< lp->getMoveStrRel() << std::endl;
	std::cout	<< "->getMoveStrAbs = "
				<< lp->getMoveStrAbs() << std::endl;
	std::cout	<< "pms.getMoveNumber(lp) = "
				<< pms.getMoveNumber(lp) << std::endl;
	
	if (ld->getName() == "sqr") {
		std::cout <<" initial structure : \n" <<toASCII(lp->getMoveSeqAbs()) <<std::endl;
	}

	// applying moves
	for (size_t i=0; i<pms.getMoveNumber(lp); i++)
	{
		LatticeProtein* neigh = lp->clone();
		neigh = pms.applyMoveInPlace(neigh,i);
		if (neigh!=NULL) {
//std::cout <<"\n orig  : ";
//for (size_t j=0; j<lp->getPointsRef()->size(); j++) {
//	std::cout <<" ("<<lp->getPointsRef()->at(j) <<"),";
//}
//std::cout <<"\n";
//std::cout <<" neigh : ";
//for (size_t j=0; j<static_cast<LatticeProtein_Ipnt*>(neigh)->getPointsRef()->size(); j++) {
//	std::cout <<" ("<<static_cast<LatticeProtein_Ipnt*>(neigh)->getPointsRef()->at(j) <<"),";
//}
//std::cout <<"\n";
			std::cout 	<< "pms.applyMoveInPlace(lp," << i << ") = rel "
			<< neigh->getMoveStrRel() <<" / abs "
			<< neigh->getMoveStrAbs()
			<< " " << (neigh->isSelfAvoiding() ? "selfavoiding" : "NOT selfavoiding!")
			<< std::endl;
			if (ld->getName() == "sqr") {
				std::cout <<" structure : \n" <<toASCII(neigh->getMoveSeqAbs()) <<std::endl;
			}
			std::cout 	<< "pms.undoLastMove(lp)";
			neigh = pms.undoLastMove(neigh);
			std::cout <<" ->getMoveStr..() = rel "
					<< neigh->getMoveStrRel() <<" / abs "
					<< neigh->getMoveStrAbs()
			<< std::endl;
			std::cout <<" undo equals origin : " 
					<<(lp->getMoveStrAbs()==neigh->getMoveStrAbs()?"true":"false")
					<<std::endl;
//std::cout <<" done  : ";
//for (size_t j=0; j<static_cast<LatticeProtein_Ipnt*>(neigh)->getPointsRef()->size(); j++) {
//	std::cout <<" ("<<static_cast<LatticeProtein_Ipnt*>(neigh)->getPointsRef()->at(j) <<"),";
//}
//std::cout <<"\n";
			delete neigh;
		}
	}

	delete lp;

	std::cout	<<"============= END OF TEST =================\n" <<std::endl;
}

void testClone() {
	std::cout	<<"========== BEGIN OF TEST ==================\n"
	<<" test cloning PullMoveSet\n"
	<<"===========================================" <<std::endl;
	std::cout << "=> creating helper objects" << std::endl;
	LatticeDescriptorCUB ld;
	LatticeModel lm(&ld);
	Alphabet alph("HP",1);
	EnergyMatrix E = EnergyMatrix(	alph.getAlphabetSize(),
									alph.getAlphabetSize(),
									0.0 );
	E[alph.getIndex("H")][alph.getIndex("H")] = -1.0;
	ContactEnergyFunction f(&alph, &E, &lm);

	std::cout 	<< "=> creating PullMoveSet(new LatticeModel(new LatticeDescriptorCUB()))"
	<< std::endl;
	PullMoveSet* pms = new PullMoveSet(&lm);

	std::string seqStr2 = "PHPHHH";
	std::cout	<<"=> "<<seqStr2 <<" = seqStr2" <<std::endl;
	std::string relMoveStr = "FLLRF";
	std::cout	<<"=> "<<relMoveStr <<"  = relMoveStr" <<std::endl;
	bool seqShared = false;
	std::cout	<<"=> "<<seqShared <<"  = seqShared" <<std::endl;
	bool isAbsoluteMoveStr = false;
	std::cout	<<"=> "<<isAbsoluteMoveStr <<"  = isAbsoluteMoveStr" <<std::endl;

	std::cout	<<"=> creating biu::LatticeProtein_Ipnt(lm,f,seqStr2,seqShared,"
				<<"relMoveStr,isAbsoluteMoveStr)"
				<<" = lp"
				<<std::endl;
	biu::Sequence seq = alph.getSequence(seqStr2);
	LatticeProtein_Ipnt* lp = new LatticeProtein_Ipnt(&lm,&f,&seq,seqShared,
									relMoveStr,isAbsoluteMoveStr);
	std::cout<< "->getMoveStrRel = " << lp->getMoveStrRel() << std::endl;
	std::cout<< "->isSelfAvoiding() = " << lp->isSelfAvoiding() << std::endl;

	// clone MoveSet!
	LatticeMoveSet* msClone = pms->clone();
	delete pms;

	size_t i=5;

	msClone->applyMoveInPlace(lp,i);

	if (lp!=NULL) {
		std::cout 	<< "msClone->applyMoveInPlace(lp," << i << ") = "
		<< lp->getMoveStrRel()
		<< std::endl;
		std::cout 	<< "msClone->undoLastMove(lp)->getMoveStrRel() = "
		<< msClone->undoLastMove(lp)->getMoveStrRel()
		<< std::endl;

		delete lp;
	}

	delete msClone;

	std::cout	<<"============= END OF TEST =================\n" <<std::endl;
}

void testReversibility() {
	std::cout	<<"========== BEGIN OF TEST ==================\n"
	<<" test PullMove reversibility\n"
	<<"===========================================" <<std::endl;
	std::cout << "=> creating helper objects" << std::endl;
	LatticeDescriptorCUB ld;
	LatticeModel lm(&ld);
	Alphabet alph("HP",1);
	EnergyMatrix E = EnergyMatrix(	alph.getAlphabetSize(),
									alph.getAlphabetSize(),
									0.0 );
	E[alph.getIndex("H")][alph.getIndex("H")] = -1.0;
	ContactEnergyFunction f(&alph, &E, &lm);

	std::cout 	<< "=> creating PullMoveSet(new LatticeModel(new LatticeDescriptorCUB()))"
	<< std::endl;
	PullMoveSet* pms = new PullMoveSet(&lm);

	std::string seqStr2 = "PHPHHHHHH";
	std::cout	<<"=> "<<seqStr2 <<" = seqStr2" <<std::endl;
	std::string relMoveStr = "FFRLRRFF";
	std::cout	<<"=> "<<relMoveStr <<"  = relMoveStr" <<std::endl;
	bool seqShared = false;
	std::cout	<<"=> "<<seqShared <<"  = seqShared" <<std::endl;
	bool isAbsoluteMoveStr = false;
	std::cout	<<"=> "<<isAbsoluteMoveStr <<"  = isAbsoluteMoveStr" <<std::endl;

	std::cout	<<"=> creating biu::LatticeProtein_Ipnt(lm,f,seqStr2,seqShared,"
				<<"relMoveStr,isAbsoluteMoveStr)"
				<<" = lp"
				<<std::endl;
	biu::Sequence seq = alph.getSequence(seqStr2);
	LatticeProtein_Ipnt* lp = new LatticeProtein_Ipnt(&lm,&f,&seq,seqShared,
									relMoveStr,isAbsoluteMoveStr);
	print3DData(lp);
	std::cout<< "->getMoveStrRel = " << lp->getMoveStrRel() << std::endl;
	std::cout<< "->isSelfAvoiding() = " << lp->isSelfAvoiding() << std::endl;

	//
	// Test Reversibility
	//
	// for each neighbor to lp \ straight Pulls
	for (size_t i=51; i<pms->getMoveNumber(lp); i++) {
		LatticeProtein_Ipnt* neigh =
			dynamic_cast<LatticeProtein_Ipnt*>(pms->applyMove(lp, i));
		if (neigh!=NULL) {
			bool isReversed = false;
			for (size_t j=0; j<pms->getMoveNumber(neigh); j++) {
				// create all neighbors to neigh and check if move was reversed
				LatticeProtein_Ipnt* neigh2 =
					dynamic_cast<LatticeProtein_Ipnt*>(pms->applyMove(neigh, j));

				if (neigh2!=NULL && *neigh2 == *lp) isReversed = true;
				delete neigh2;
			}
			if (isReversed==false) {
				std::cout<< "sequence derived by move " << i << ": "
				<<neigh->getMoveStrRel() << " can't be reversed to "
				<< lp->getMoveStrRel() << std::endl;
			}
		delete neigh;
		}
	}

	delete lp;
	delete pms;

	std::cout	<<"============= END OF TEST =================\n" <<std::endl;
}

void
testCopyConstructor() {
	std::cout	<<"========== BEGIN OF TEST ==================\n"
				<<" test copy constructor\n"
				<<"===========================================" <<std::endl;
	std::cout << "=> creating helper objects" << std::endl;
	LatticeDescriptorCUB ld;
	LatticeModel lm(&ld);
	Alphabet alph("HP",1);
	EnergyMatrix E = EnergyMatrix(	alph.getAlphabetSize(),
									alph.getAlphabetSize(),
									0.0 );
	E[alph.getIndex("H")][alph.getIndex("H")] = -1.0;
	ContactEnergyFunction f(&alph, &E, &lm);

	std::cout 	<< "=> creating PullMoveSet(new LatticeModel(new LatticeDescriptorCUB()))"
	<< std::endl;
	PullMoveSet* pms = new PullMoveSet(&lm);

	std::string seqStr2 = "PHPHHH";
	std::cout	<<"=> "<<seqStr2 <<" = seqStr2" <<std::endl;
	std::string relMoveStr = "FLLRF";
	std::cout	<<"=> "<<relMoveStr <<"  = relMoveStr" <<std::endl;
	bool seqShared = false;
	std::cout	<<"=> "<<seqShared <<"  = seqShared" <<std::endl;
	bool isAbsoluteMoveStr = false;
	std::cout	<<"=> "<<isAbsoluteMoveStr <<"  = isAbsoluteMoveStr" <<std::endl;

	std::cout	<<"=> creating biu::LatticeProtein_Ipnt(lm,f,seqStr2,seqShared,"
				<<"relMoveStr,isAbsoluteMoveStr)"
				<<" = lp"
				<<std::endl;
	biu::Sequence seq = alph.getSequence(seqStr2);
	LatticeProtein_Ipnt* lp = new LatticeProtein_Ipnt(&lm,&f,&seq,seqShared,
									relMoveStr,isAbsoluteMoveStr);
	std::cout<< "->getMoveStrRel = " << lp->getMoveStrRel() << std::endl;
	std::cout<< "->isSelfAvoiding() = " << lp->isSelfAvoiding() << std::endl;

	// test copy constructor
	PullMoveSet* moveSet2 = new PullMoveSet(*pms);
	delete pms;

	size_t i=5;

	moveSet2->applyMoveInPlace(lp,i);

	if (lp!=NULL) {
		std::cout 	<< "moveSet2->applyMoveInPlace(lp," << i << ") = "
		<< lp->getMoveStrRel()
		<< std::endl;
		std::cout 	<< "moveSet2->undoLastMove(lp)->getMoveStrRel() = "
		<< moveSet2->undoLastMove(lp)->getMoveStrRel()
		<< std::endl;
	}


	// test operator=
	std::cout	<< "PullMoveSet* moveSet3 = new PullMoveSet(&lm);"
				<< std::endl;
	PullMoveSet* moveSet3 = new PullMoveSet(&lm);
	std::cout	<< "moveSet3 = moveSet2;" << std::endl;
	moveSet3->operator=(*moveSet2);

	delete moveSet2;

	std::cout 	<< "moveSet3->undoLastMove(lp)->getMoveStrRel() = "
				<< moveSet3->undoLastMove(lp)->getMoveStrRel()
				<< std::endl;

	delete moveSet3;
	delete lp;

	std::cout	<<"============= END OF TEST =================\n" <<std::endl;
}

void testSharedDecoder() {
	std::cout	<<"========== BEGIN OF TEST ==================\n"
	<<" test function : PullMoveSet with shared pull move decoder\n"
	<<"===========================================" <<std::endl;
	std::cout << "=> creating helper objects" << std::endl;
	LatticeDescriptorCUB ld;
	LatticeModel lm(&ld);
	Alphabet alph("HP",1);
	EnergyMatrix E = EnergyMatrix(	alph.getAlphabetSize(),
									alph.getAlphabetSize(),
									0.0 );
	E[alph.getIndex("H")][alph.getIndex("H")] = -1.0;
	ContactEnergyFunction f(&alph, &E, &lm);

	std::cout 	<< "=> creating PullMoveSet(new LatticeModel(new LatticeDescriptorCUB()))"
	<< std::endl;
	PullMoveSet::PullMoveDecoder pmd(&lm);
	bool decoderShared = true;
	PullMoveSet pms(&pmd, decoderShared);

	//	/*
	//	 *		Test 3:
	//	 */
	std::string seqStr2 = "PHPHHHHHH";
	std::cout	<<"=> "<<seqStr2 <<" = seqStr2" <<std::endl;
	std::string relMoveStr = "FFRLRRFF";
	std::cout	<<"=> "<<relMoveStr <<"  = relMoveStr" <<std::endl;
	bool seqShared = false;
	std::cout	<<"=> "<<seqShared <<"  = seqShared" <<std::endl;
	bool isAbsoluteMoveStr = false;
	std::cout	<<"=> "<<isAbsoluteMoveStr <<"  = isAbsoluteMoveStr" <<std::endl;

	std::cout	<<"=> creating biu::LatticeProtein_Ipnt(lm,f,seqStr2,seqShared,"
				<<"relMoveStr,isAbsoluteMoveStr)"
				<<" = lp"
				<<std::endl;
	biu::Sequence seq = alph.getSequence(seqStr2);
	LatticeProtein_Ipnt* lp = new LatticeProtein_Ipnt(&lm,&f,&seq,seqShared,
									relMoveStr,isAbsoluteMoveStr);
	std::cout<< "->getMoveStrRel = " << lp->getMoveStrRel() << std::endl;
	std::cout<< "->isSelfAvoiding() = " << lp->isSelfAvoiding() << std::endl;

	// applying moves
	for (size_t i=0; i<pms.getMoveNumber(lp); i++)
	{
		LatticeProtein* neigh = lp->clone();
		pms.applyMoveInPlace(neigh,i);
		if (lp!=NULL) {
			std::cout 	<< "pms.applyMoveInPlace(lp," << i << ") = "
			<< neigh->getMoveStrRel()
			<< std::endl;
			std::cout 	<< "pms.undoLastMove(lp)->getMoveStrRel() = "
			<< pms.undoLastMove(neigh)->getMoveStrRel()
			<< std::endl;
			delete neigh;
		}
	}

	delete lp;

	std::cout	<<"============= END OF TEST =================\n" <<std::endl;
}

void testPullMoveSet_selfavoiding()
{
	std::cout	<<"========== BEGIN OF TEST ==================\n"
	<<" test function : PullMoveSet\n"
	<<"===========================================" <<std::endl;
	std::cout << "=> creating helper objects" << std::endl;
	LatticeDescriptorCUB ld;
	LatticeModel lm(&ld);
	Alphabet alph("HP",1);
	EnergyMatrix E = EnergyMatrix(	alph.getAlphabetSize(),
									alph.getAlphabetSize(),
									0.0 );
	E[alph.getIndex("H")][alph.getIndex("H")] = -1.0;
	ContactEnergyFunction f(&alph, &E, &lm);

	std::cout 	<< "=> creating PullMoveSet(new LatticeModel(new LatticeDescriptorCUB()))"
	<< std::endl;
	PullMoveSet pms(&lm);

	//	/*
	//	 *		Test 3:
	//	 */
	std::string seqStr2 = "HHHHHHHPHHPHHPHHPHPHPHPHHPH";
	std::cout	<<"=> "<<seqStr2 <<" = seqStr2" <<std::endl;
	std::string absMoveStr = "BDFFULDFRFFDFURFLFLUBURDFR";
	std::cout	<<"=> "<<absMoveStr <<"  = absMoveStr" <<std::endl;
	bool seqShared = false;
	std::cout	<<"=> "<<seqShared <<"  = seqShared" <<std::endl;
	bool isAbsoluteMoveStr = true;
	std::cout	<<"=> "<<isAbsoluteMoveStr <<"  = isAbsoluteMoveStr" <<std::endl;

	std::cout	<<"=> creating biu::LatticeProtein_Ipnt(lm,f,seqStr2,seqShared,"
				<<"absMoveStr,isAbsoluteMoveStr)"
				<<" = lp"
				<<std::endl;
	biu::Sequence seq = alph.getSequence(seqStr2);
	LatticeProtein_Ipnt* lp = new LatticeProtein_Ipnt(&lm,&f,&seq,seqShared,
									absMoveStr,isAbsoluteMoveStr);
	std::cout<< "->getMoveStrRel = " << lp->getMoveStrRel() << std::endl;
	std::cout<< "->isSelfAvoiding() = " << lp->isSelfAvoiding() << std::endl;

	// applying moves
	for (size_t i=0; i<pms.getMoveNumber(lp); i++)
	{
		LatticeProtein* neigh = lp->clone();
		pms.applyMoveInPlace(neigh,i);
		if (lp!=NULL) {
			std::cout 	<< "pms.applyMoveInPlace(lp," << i << ") = "
			<< neigh->getMoveStrRel()
			<< " " << (neigh->isSelfAvoiding() ? "selfavoiding" : "NOT selfavoiding!")
			<< std::endl;
			if (!(pms.undoLastMove(neigh)->getMoveStrAbs() == absMoveStr))
				std::cout << "Error: undo failed" << std::endl;
			delete neigh;
		}
	}

	delete lp;

	std::cout	<<"============= END OF TEST =================\n" <<std::endl;
}

int main(int argc, char **argv)
{
	
	std::cout	<< "Testing CUB:" << std::endl;
	LatticeDescriptorCUB ld_CUB;
	testPullMoveDecoderNeighborhood(&ld_CUB);
	testPullMoveDecoderLookup(&ld_CUB, "PHPHHHH", "FFFFFF");

	testPullMoveSet(&ld_CUB, "PHPHHHH", "FFFFFF", false);
	testPullMoveSet(&ld_CUB, "PHPHHH", "FLLRF", false);
	testPullMoveSet(&ld_CUB, "PHPHHHHHH", "FFRLRRFF", false);

	testPullMoveSet_selfavoiding();
	testClone();
	testCopyConstructor();
	testSharedDecoder();

	testReversibility();

	std::cout	<< std::endl << "Testing FCC:" << std::endl;
	LatticeDescriptorFCC ld_FCC;
	testPullMoveDecoderNeighborhood(&ld_FCC);
	testPullMoveDecoderLookup(&ld_FCC, "PHPHHHH", "FLFLFLFLFLFL");
	testPullMoveSet(&ld_FCC, "PHPHHHHHH", "FLFLFLFLFLFLFLFL", false);

	std::cout	<< std::endl << "Testing SQR:" << std::endl;
	LatticeDescriptorSQR ld_SQR;
	testPullMoveDecoderNeighborhood(&ld_SQR);
	testPullMoveDecoderLookup(&ld_SQR, "HPHPPHP", "FRRBLB");
	testPullMoveSet(&ld_SQR, "HPHPPHP","FRRBLB",true);
	
	return 0;
}

