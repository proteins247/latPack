
#include <iostream>
#include <biu/LatticeDescriptorCUB.hh>
#include <biu/LatticeModel.hh>


int
main() {
	using namespace biu;


	std::cout	<<"========== BEGIN OF TEST =================\n"
				<<" test class : LatticeModel\n"
				<<"==========================================" <<std::endl;
				
		// erzeugen	
	std::cout <<"=> creating biu::LatticeDescriptorCUB() = ld" <<std::endl;
		LatticeDescriptorCUB ld ;
	std::cout <<"=> creating biu::LatticeModel(ld)" <<std::endl;
		LatticeModel lm(&ld);
		
		// getDescriptor()
	std::cout	<<"->getDescriptor()->getName()  = " 
				<<lm.getDescriptor()->getName()
				<<std::endl;

		// areNeighbored()
	std::cout	<<"->areNeighbored(IntPoint(0,0,0), IntPoint(2,2,2))  = " 
				<<(lm.areNeighbored(IntPoint(0,0,0), IntPoint(2,2,2))?
					"true":"false")
				<<std::endl;
	std::cout	<<"->areNeighbored(IntPoint(0,0,0), IntPoint(0,0,1))  = " 
				<<(lm.areNeighbored(IntPoint(0,0,0), IntPoint(0,0,1))?
					"true":"false")
				<<std::endl;

		// getString()
		// parseMoveString()
	std::string moveStr = "FFLUURBLDD";
	std::cout	<<"-> " <<moveStr <<" = moveStr" <<std::endl;
	std::cout	<<"->(getString(lm.parseMoveString(moveStr)) == moveStr)  = "
				<<(lm.getString(lm.parseMoveString(moveStr)) == moveStr?
					"true":"false")
				<<std::endl;
				
	IPointVec points;
		points.push_back(biu::IntPoint(0,0,0));	
		points.push_back(biu::IntPoint(0,1,0));	
		points.push_back(biu::IntPoint(0,2,0));	
		points.push_back(biu::IntPoint(1,2,0));	
		points.push_back(biu::IntPoint(2,2,0));
	std::cout	<<"-> IPointVec(...) = points" <<std::endl;
	
		// absMovesToPoints()
		// pointsToAbsMoves()
	std::cout	<<"->(absMovesToPoints(lm.pointsToAbsMoves(points))"
				<<" == points)  = " 
				<< (lm.absMovesToPoints(lm.pointsToAbsMoves(points)) == points?
						"true":"false")
				<<std::endl;
				
		// relMovesToPoints()
		// pointsToRelMoves()
	std::cout	<<"->(relMovesToPoints(lm.pointsToRelMoves(points))"
				<<" == points)  = " 
				<< (lm.relMovesToPoints(lm.pointsToRelMoves(points)) == points?
						"true":"false")
				<<std::endl;

	IPointVec toFill;
	std::cout	<<"->(relMovesToPoints(lm.pointsToRelMoves(points), toFill)"
				<<std::endl;
	lm.relMovesToPoints(lm.pointsToRelMoves(points), toFill);
	IPointVec::iterator it;
	std::cout	<< "toFill =" << std::endl;
	for (it = toFill.begin(); it != toFill.end(); it++) {
		std::cout << *it << std::endl;;
	}
		
				
		// getAbsMove()
	std::cout	<<"->getString(MoveSequence(1,getAbsMove(IntPoint(0,0,0), IntPoint(0,0,1))))  = " 
				<<lm.getString(MoveSequence(1,lm.getAbsMove(IntPoint(0,0,0), IntPoint(0,0,1))))
				<<std::endl;

		// getAllNeighPoints()
	std::cout	<<"->getAllNeighPoints(IntPoint(0,0,0)).size()  = " 
				<<lm.getAllNeighPoints(IntPoint(0,0,0)).size()
				<<std::endl;
	
	IPointSet neighbors = lm.getAllNeighPoints(IntPoint(0,0,0));
	bool retVal = true;
	for (IPointSet::const_iterator it = neighbors.begin(); 
			it != neighbors.end(); it++) 
	{
		retVal &= lm.getDescriptor()->getNeighborhood().isElement(*it);
	}
	std::cout	<<"-> test of return values of getAllNeighPoints  = "
				<<(retVal?"ok":"failure")
				<< std::endl;
	
		// applyAbsMove
	std::cout	<<"->applyAbsMove(IntPoint(2,2,2),"
				<<" lm.parseMoveString(U)[0])  = (" 
				<<lm.applyAbsMove(IntPoint(2,2,2), lm.parseMoveString("U")[0])
				<<")"
				<<std::endl;

	// ----------------------------------------
	// test conversion rel <-> abs moves
	std::cout	<<"->(lm.getString( lm.absMovesToRelMoves( "
				<<"lm.relMovesToAbsMoves( lm.parseMoveString( moveStr))))"
				<<" == moveStr)  = "
				<<(lm.getString( 
					lm.absMovesToRelMoves(
						lm.relMovesToAbsMoves(lm.parseMoveString(moveStr)) ) )
					== moveStr ? "true":"false")
				<< std::endl;
	
	// test copy constructor
	std::cout	<<"LatticeModel lm2 = LatticeModel(lm);" << std::endl;
	LatticeModel lm2 = LatticeModel(lm);
	std::cout	<<"->applyAbsMove(IntPoint(2,2,2),"
					<<" lm2.parseMoveString(U)[0])  = (" 
					<<lm2.applyAbsMove(IntPoint(2,2,2), lm2.parseMoveString("U")[0])
					<<")"
					<<std::endl;
	
	// test operator==
	std::cout	<< "(lm==lm2) = " 
				<<((lm==lm2) == true ? "true":"false") 
				<< std::endl;

	std::cout	<<"============= END OF TEST ================\n" <<std::endl;
	
	return 0;
}
