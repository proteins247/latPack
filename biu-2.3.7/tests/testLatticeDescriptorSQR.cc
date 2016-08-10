#include <biu/LatticeDescriptorSQR.hh>

int
main() {
	
	biu::LatticeDescriptorSQR ldsqr;

	std::cout <<ldsqr.getName()<<std::endl;

	biu::AutomorphismVec as = ldsqr.getAutomorphisms();
	
	std::cout <<as.size()<<std::endl;
	std::cout << std::endl;

	for (unsigned int i=0; i<as.size(); i++) {
		biu::Automorphism &a = as[i];
		
		std::cout << a << std::endl;
	}
}
