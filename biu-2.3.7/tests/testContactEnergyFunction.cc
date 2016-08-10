#include <iostream>
//#include <biu/ContactEnergyFunction.hh>
#include <biu/DistanceEnergyFunction.hh>

#include <biu/LatticeDescriptorCUB.hh>


int main(int argc, char **argv) {

	using namespace biu;

	std::cout	<<"========== BEGIN OF TEST =================\n"
				<<" test class : ContactEnergyFunction\n"
				<<"==========================================" <<std::endl;
				
		// erzeugen	
	std::cout <<"=> creating biu::Alphabet(HPNX,1) = alph" <<std::endl;
	Alphabet alph("HPNX",1);
	std::cout <<"=> creating biu::EnergyMatrix(4,4,0.0) = E" <<std::endl;
	EnergyMatrix E = EnergyMatrix(	alph.getAlphabetSize(),
									alph.getAlphabetSize(),
									0.0 );
		E[alph.getIndex("H")][alph.getIndex("H")] = -4.0;
		E[alph.getIndex("P")][alph.getIndex("P")] =  1.0;
		E[alph.getIndex("N")][alph.getIndex("N")] =  1.0;
		E[alph.getIndex("P")][alph.getIndex("N")] = -1.0;
		E[alph.getIndex("N")][alph.getIndex("P")] = -1.0;
	std::cout <<E <<std::endl;
	std::cout <<"=> creating biu::LatticeModel( CUB )" <<std::endl;
	biu::LatticeDescriptorCUB latDesc;
	biu::LatticeModel lattice(&latDesc);
	std::cout <<" lattice = " <<lattice.getDescriptor()->getName() <<std::endl;
	std::cout <<"=> creating biu::ContactEnergyFunction( alph, E, lattice)" <<std::endl;
	ContactEnergyFunction f(&alph, &E, &lattice);
	
		// getAlphabet()
	std::cout	<<"->getAlphabet()->getAlphabetSize() = " 
				<<f.getAlphabet()->getAlphabetSize() 
				<<std::endl;
	
		// getContactEnergy()
	std::cout <<"->getContactEnergy( 0..3, 0..3)" <<std::endl;
	EnergyMatrix T = EnergyMatrix( E.numRows(), E.numColumns(), -1.0);
	for (unsigned int i = 0; i<alph.getAlphabetSize(); i++) {
		for (unsigned int j = 0; j<alph.getAlphabetSize(); j++) {
			std::cout <<f.getContactEnergy(	alph.getElement(i), 
											alph.getElement(j)) <<" ";
		}
		std::cout <<std::endl;
	}

		// ==
	std::cout	<<"-> (*this == *this) = " 
				<<(f == f?"true":"false")
				<<std::endl;

		// !=
	std::cout	<<"-> (*this != *this) = " 
				<<(f != f?"true":"false")
				<<std::endl;

	std::cout	<<"============= END OF TEST ================\n" <<std::endl;
	
	return 0;
}
