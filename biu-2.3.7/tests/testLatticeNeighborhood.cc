

#include <iostream>
#include <biu/LatticeNeighborhood.hh>

int main(int argc, char **argv) {

	using namespace biu;

	std::cout	<<"========== BEGIN OF TEST =================\n"
				<<" test class : LatticeNeighborhood\n"
				<<"==========================================" <<std::endl;

		// erzeugen
	std::cout <<"=> creating biu::MoveAlphabet(FBLR,1) = alph" <<std::endl;
		MoveAlphabet* alph = new MoveAlphabet("FBLR",1);	// Square lattice
	std::cout <<"=> creating corresponding NeighSet = neighbors" <<std::endl;
		NeighSet neighbors;
		// nachbarschaftsvektoren
		// AUTOMORPHISMEN NUR EXEMPLARISCH !! NICHT UNBEDINGT RICHTIG !!!

		int m1[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
		int m2[3][3] = {{-1,0,0},{0,-1,0},{0,0,1}};
		int m3[3][3] = {{0,1,0},{-1,0,0},{0,0,1}};
		int m4[3][3] = {{0,-1,0},{1,0,0},{0,0,1}};

		neighbors.insert(NeighborVector(+1,0,0,alph->getElement("F"),
			Automorphism(m1),	// rel2abs
			Automorphism(m1)));	// abs2rel
		neighbors.insert(NeighborVector(-1,0,0,alph->getElement("B"),
			Automorphism(m2),	// rel2abs
			Automorphism(m2)));	// abs2rel
		neighbors.insert(NeighborVector(0,+1,0,alph->getElement("L"),
			Automorphism(m3),	// rel2abs
			Automorphism(m3)));	// abs2rel
		neighbors.insert(NeighborVector(0,-1,0,alph->getElement("R"),
			Automorphism(m4),	// rel2abs
			Automorphism(m4)));	// abs2rel
	std::cout <<"=> creating biu::LatticeNeighborhood(alph,neighbors)" <<std::endl;
		LatticeNeighborhood nh(alph, neighbors);

		// size()
	std::cout	<<"->size()  = " << nh.size()
				<<std::endl;

		// isElement()
	std::cout	<<"->isElement(IntPoint(3,3,3))  = "
				<<(nh.isElement(IntPoint(3,3,3))?"true":"false")
				<<std::endl;

		// test neighbor in nh
		// isElement()
		// getElement(IntPoint)
		// getElement(Move)
	std::cout	<<"-> test elements of neighbor in LatticeNeighborhood:"
				<<std::endl;
	bool isElem = true, getElem1 = true, getElem2 = true;
	for (NeighSet::const_iterator it = neighbors.begin();
			it != neighbors.end(); it++)
	{
		isElem &= nh.isElement(*it);
		getElem1 &= nh.getElement(*it) == *it;
		getElem2 &= nh.getElement(it->getMove()) == *it;
	}
	std::cout	<<"->isElement(..) = " <<(isElem?"all ok":"failed") <<std::endl;
	std::cout	<<"->getElement(IntPoint(..)) = " <<(getElem1?"all ok":"failed")
				<<std::endl;
	std::cout	<<"->getElement(Move(..)) = " <<(getElem2?"all ok":"failed")
				<<std::endl;

	bool retVal = true;
	for (LatticeNeighborhood::const_iterator it = nh.begin();
			it != nh.end(); it++) {
		retVal &= neighbors.find(*it) != neighbors.end();	// iterator testen
	}
	std::cout	<<"-> all iterator elements of LatticeNeighborhood are"
				<<" inside neighbors  = "
				<<(retVal?"true":"false")
				<<std::endl;


	delete alph;

	std::cout	<<"============= END OF TEST ================\n" <<std::endl;

	return 0;
}
