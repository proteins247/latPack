

#include <iostream>
#include "biu/AllowedBasePairs.hh"


int main(int argc, char **argv) {

	using namespace biu;

	std::cout	<<"========== BEGIN OF TEST =================\n"
				<<" test class : AllowedBasePairs\n"
				<<"==========================================" <<std::endl;
				
		// erzeugen	
	std::cout <<"=> creating biu::Alphabet(ACGU,1) = alph" <<std::endl;
	Alphabet alph = Alphabet("ACGU",1);
	std::cout	<<"=> creating string with basepairs = {A-U, C-G, G-U} "
				<<std::endl;
	std::string bpStr = "AU,GU,CG";
	std::cout	<<" bpStr = '" <<bpStr <<"'" <<std::endl;
	std::cout <<"=> creating biu::AllowedBasePairs( alph, bpStr)" <<std::endl;
	AllowedBasePairs bp(&alph, bpStr);

		// getAlphabet()	
	std::cout	<<"->getAlphabet().getSize() = " 
				<<bp.getAlphabet()->getAlphabetSize() 
				<<std::endl;

		// allowedBasePair()
	std::cout	<<"->allowedBasePair(C,A) = " 
				<< (bp.allowedBasePair(	bp.getAlphabet()->getElement("C"),
										bp.getAlphabet()->getElement("A"))
					?"true":"false")
				<<std::endl;
	std::cout	<<"->allowedBasePair(C,G) = " 
				<< (bp.allowedBasePair(	bp.getAlphabet()->getElement("C"),
										bp.getAlphabet()->getElement("G"))
					?"true":"false")
				<<std::endl;
		
	std::cout	<<"============= END OF TEST ================\n" <<std::endl;
	
	return 0;
}
