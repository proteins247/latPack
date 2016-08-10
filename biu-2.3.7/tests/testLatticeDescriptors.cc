


#include <iostream>
#include <biu/LatticeDescriptor.hh>
#include <biu/LatticeDescriptorSQR.hh>
#include <biu/LatticeDescriptorCUB.hh>
#include <biu/LatticeDescriptorFCC.hh>


void test_classLatticeDescriptor(const biu::LatticeDescriptor& ld,
		const std::string& moveStr ) {

	using namespace biu;

		// getName()
	std::cout	<<"->getName()  = " << ld.getName()
				<<std::endl;

	std::cout	<<"-> moveStr = " <<moveStr
				<<"\n-> moveStr.size = "  <<moveStr.size() 
				<<std::endl;
								
		// getSequence()
	std::cout	<<"->getSequence(moveStr).size()  = " 
				<< ld.getSequence(moveStr).size()
				<<std::endl;

		// getString(getSequence(X))
	std::cout	<<"->getString(ld.getSequence(moveStr)) == moveStr = "
				<< (ld.getString(ld.getSequence(moveStr)) == moveStr?
					"true":"false")
				<<std::endl;

		// getBase()
	std::cout	<<"->getBase().size() = " <<ld.getBase().size()
				<<std::endl;

		// getNeighborhood()
	std::cout	<<"->getNeighborhood().size()  = " 
				<< ld.getNeighborhood().size()
				<<std::endl;
	
		// getAutomorphisms()
	std::cout	<<"->getAutomorphisms().size()  = " 
				<< ld.getAutomorphisms().size()
				<<std::endl;
	
		// normalizeSequence()
	std::cout	<<"->normalizeSequence(" <<moveStr <<") = "
				<<ld.getString(ld.normalizeSequence(ld.getSequence(moveStr)))
				<<std::endl;
	
		// getAllSymmetricSequences()
	std::cout	<<"->getAllSymmetricSequences(" <<moveStr <<") + normalize : "
				<<std::endl;
	typedef std::set<biu::MoveSequence> MSset;
	MSset allSymm = ld.getAllSymmetricSequences(ld.getSequence(moveStr));
	for (MSset::const_iterator it=allSymm.begin(); it!=allSymm.end(); it++) {
		std::cout	<<"  '" <<ld.getString(*it) <<"' == '"
					<<ld.getString(ld.normalizeSequence(*it)) <<"'"
					<<std::endl;
	}
	
		// isLatticeNode(...)
	biu::IPointVec x;
	x.push_back( IntPoint(0,0,0) );
	x.push_back( IntPoint(1,0,0) );
	x.push_back( IntPoint(0,1,0) );
	x.push_back( IntPoint(0,0,1) );
	x.push_back( IntPoint(1,1,0) );
	x.push_back( IntPoint(0,1,1) );
	x.push_back( IntPoint(1,0,1) );
	x.push_back( IntPoint(1,1,1) );
	for (size_t i=0; i<x.size(); i++) {
		std::cout <<"->isLatticeNode(" <<x.at(i) <<") = "
				<<(ld.isLatticeNode(x.at(i))?"true":"false")
				<<std::endl;
		std::cout <<"->isLatticeNode(" <<(x.at(i)*2) <<") = "
				<<(ld.isLatticeNode((x.at(i)*2))?"true":"false")
				<<std::endl;
	}
	
	
	for (size_t i=0; i<10; i++) {
		std::cout <<"->isPossibleRing(" <<i<<") = "
				<<(ld.isPossibleRing(i)?"true":"false")
				<<std::endl;
	}

}


int main(int argc, char **argv) {

	using namespace biu;

	std::cout	<<"========== BEGIN OF TEST =================\n"
				<<" test class : LatticeDescriptor\n"
				<<"==========================================" <<std::endl;
				
		// erzeugen	

	using namespace biu;

		//cubic
	std::cout <<"\n=> creating biu::LatticeDescriptorSQR()" <<std::endl;
		test_classLatticeDescriptor(LatticeDescriptorSQR(), "FLBR");

		// square
	std::cout <<"\n=> creating biu::LatticeDescriptorCUB()" <<std::endl;
		test_classLatticeDescriptor(LatticeDescriptorCUB(), "FLUBRD");

		// fcc
	std::cout <<"\n=> creating biu::LatticeDescriptorFCC()" <<std::endl;
		test_classLatticeDescriptor(LatticeDescriptorFCC(), 
			"FLFRFUFDBLBRBUBDLULDRDRD");

		// == and !=
	std::cout	<<"\n-> (LatticeDescriptorCUB() == LatticeDescriptorCUB())  = "
				<<(LatticeDescriptorCUB() == LatticeDescriptorCUB()?
					"true":"false")
				<<std::endl;
	std::cout	<<"-> (LatticeDescriptorCUB() != LatticeDescriptorSQR())  = "
				<<(LatticeDescriptorCUB() != LatticeDescriptorSQR()?
					"true":"false")
				<<std::endl;
		
	std::cout	<<"============= END OF TEST ================\n" <<std::endl;
	
	return 0;
}
