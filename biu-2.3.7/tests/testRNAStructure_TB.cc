
#include <iostream>

#include "biu/RNAStructure_TB.hh"


int main(int arg, char** argv) {
	

	using namespace biu;

	std::cout	<<"========== BEGIN OF TEST =================\n"
				<<" test class : RNAStructure_TB \n"
				<<"==========================================" <<std::endl;
	
		// create structural alphabet and allowed base pairs 	
	std::cout <<"=> creating biu::Alphabet(ACGU,1) = alph" <<std::endl;
	Alphabet alph = Alphabet("ACGU",1);
	std::cout	<<"=> creating string with basepairs = {A-U, C-G, G-U} "
				<<std::endl;
	std::string bpStr = "AU,GU,CG";
	std::cout	<<" bpStr = '" <<bpStr <<"'" <<std::endl;
	std::cout <<"=> creating biu::AllowedBasePairs( alph, bpStr)" <<std::endl;
	AllowedBasePairs bp(&alph, bpStr);

	{
		std::string rnaSeq = "ACGUGCUGUAAGUCCCGUUU";
		std::string rnaStr = "..(...)..(..(...).).";
		
		std::cout <<"-->RNAStructure_TB rna(" <<rnaSeq <<"," <<rnaStr<<")" <<std::endl;
		RNAStructure_TB rna(rnaSeq,rnaStr,&bp);
		
		size_t i = RNAStructure_TB::INVALID_INDEX;
		size_t j = RNAStructure_TB::INVALID_INDEX;
		
		std::cout <<"\n-->rna.getNextSingleMove(...) : " <<std::endl;
		while (rna.getNextSingleMove(i,j)) {
			std::cout <<" next = (" <<i<<","<<j<<")" <<std::endl;
		}
		std::cout <<std::endl;
	}
	
	{
		std::string rnaSeq = "AGGCUUGUAGCUCAGGUGGUUAGAGCGCACCCCUGAUAAGGGUGAGGUCGGUGGUUCAAGUCCACUCAGGCCUACCA";
		std::string rnaStr = "(((..(((.((((.........)))))))..)))......(((.((((((((((.......)))))..)))))))).";
		
		std::cout <<"-->RNAStructure_TB rna(" <<rnaSeq <<"," <<rnaStr<<")" <<std::endl;
		RNAStructure_TB rna(rnaSeq,rnaStr,&bp);
		
		size_t i = RNAStructure_TB::INVALID_INDEX;
		size_t j = RNAStructure_TB::INVALID_INDEX;
		
		std::cout <<"\n-->rna.getNextSingleMove(...) : " <<std::endl;
		while (rna.getNextSingleMove(i,j)) {
			std::cout <<" next = (" <<i<<","<<j<<")" <<std::endl;
		}
		std::cout <<std::endl;
	}

	

	
	
	
	std::cout	<<"============= END OF TEST ================\n" <<std::endl;

	return 0;
}
