#include <biu/LatticeDescriptorCUB.hh>

int
main() {
	
	biu::LatticeDescriptorCUB ld;

	std::cout <<ld.getName()<<std::endl;

	biu::AutomorphismVec as = ld.getAutomorphisms();
	
	std::cout <<as.size()<<std::endl;
	std::cout << std::endl;

	for (unsigned int i=0; i<as.size(); i++) {
		biu::Automorphism &a = as[i];
		
		std::cout << a << std::endl;
	}
	
	std::vector< std::string > absMoves;
	absMoves.push_back("RDFLUUFRDDFRLLLFR");
	absMoves.push_back("RBULFFURBBURLLLUR");
	absMoves.push_back("RFULBBURFFURLLLUR");
	absMoves.push_back("RUFLDDFRUUFRLLLFR");
	absMoves.push_back("DRFULLFDRRFDUUUFD");
	absMoves.push_back("DBLUFFLDBBLDUUULD");
	absMoves.push_back("DFLUBBLDFFLDUUULD");
	absMoves.push_back("DLFURRFDLLFDUUUFD");
	absMoves.push_back("BRUFLLUBRRUBFFFUB");
	absMoves.push_back("BDLFUULBDDLBFFFLB");
	absMoves.push_back("BULFDDLBUULBFFFLB");
	absMoves.push_back("BLUFRRUBLLUBFFFUB");
	absMoves.push_back("FRUBLLUFRRUFBBBUF");
	absMoves.push_back("FDLBUULFDDLFBBBLF");
	absMoves.push_back("FULBDDLFUULFBBBLF");
	absMoves.push_back("FLUBRRUFLLUFBBBUF");
	absMoves.push_back("URFDLLFURRFUDDDFU");
	absMoves.push_back("UBLDFFLUBBLUDDDLU");
	absMoves.push_back("UFLDBBLUFFLUDDDLU");
	absMoves.push_back("ULFDRRFULLFUDDDFU");
	absMoves.push_back("LDFRUUFLDDFLRRRFL");
	absMoves.push_back("LBURFFULBBULRRRUL");
	absMoves.push_back("LFURBBULFFULRRRUL");
	absMoves.push_back("LUFRDDFLUUFLRRRFL");
	std::cout <<"normalizing : \n";
	for (size_t i=0; i<absMoves.size(); i++) {
		std::cout <<"   " <<absMoves[i] <<" --> "
				<<ld.getString(ld.normalizeSequence(ld.getSequence(absMoves[i]))) 
				<<"\n";
	}
	std::cout <<std::endl;
}
