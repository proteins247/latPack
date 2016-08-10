

#include <iostream>
#include <biu/NeighborVector.hh>


int main(int argc, char **argv) {

	using namespace biu;

	std::cout	<<"========== BEGIN OF TEST =================\n"
				<<" test class : NeighborVector\n"
				<<"==========================================" <<std::endl;

		// erzeugen
	std::cout <<"=> creating biu::MoveAlphabet(FBLR,1)" <<std::endl;
	MoveAlphabet alph = MoveAlphabet("FBLR",1);
	std::cout <<"=> creating biu::Automorphism({{1,0,0},{0,1,0},{0,0,1}})"
				<<" = a2r" <<std::endl;
	int a2array[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
	Automorphism a2r = Automorphism(a2array);
	std::cout <<"=> creating a2r*2 = r2a" <<std::endl;
	Automorphism r2a = a2r * 2;
	std::cout <<"=> creating biu::NeighborVector(1,0,0, alph.getElement(B)"
				<<", r2a, a2r)" <<std::endl;
	NeighborVector nv = NeighborVector(-1,0,0, alph.getElement("B"), r2a, a2r);

		// getX()
	std::cout	<<"->getX()  = "  << nv.getX() <<std::endl;

		// getY()
	std::cout	<<"->getY()  = " << nv.getY() <<std::endl;

		// getZ()
	std::cout	<<"->getZ()  = "  << nv.getZ() <<std::endl;

		// getMove()
	std::cout	<<"->alph.getString(getMove())  = "
				<< alph.getString(nv.getMove())
				<<std::endl;

		// getAbs2RelRotation()
	std::cout	<<"->getAbs2RelRotation()  = \n"  << nv.getAbs2RelRotation();
	std::cout.flush();

		// getRel2AbsRotation()
	std::cout	<<"->getRel2AbsRotation()  = \n"  << nv.getRel2AbsRotation();
	std::cout.flush();

		// operator ==
	std::cout	<<"->this == IntPoint(-1,0,0)  = "
				<< (nv == IntPoint(-1,0,0)?"true":"false")
				<<std::endl;



	std::cout	<<"============= END OF TEST ================\n" <<std::endl;

	return 0;
}
