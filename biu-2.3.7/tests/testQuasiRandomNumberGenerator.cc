#include "biu/QuasiRandomNumberGenerator.hh"
#include <iostream>


using namespace biu;

int main(int argc, char **argv) {

	std::cout	<<"========== BEGIN OF TEST =================\n"
	<<" test class : QuasiRandomNumberGenerator\n"
	<<"==========================================" <<std::endl;

	std::cout	<<"=> creating QuasiRandomNumberGenerator" << std::endl;
	QuasiRandomNumberGenerator* qrng = new QuasiRandomNumberGenerator();
	for (int i=0; i<10; i++) {
		std::cout << "->qrng->getQuasiRN() = " << qrng->getQuasiRN() << std::endl;
	}
	
	delete qrng;
	
	std::cout	<<"* Check if two QRNG's can be used independently:" << std::endl;
	std::cout	<<"=> creating QuasiRandomNumberGenerators qrng1 & qrng2" << std::endl;
	QuasiRandomNumberGenerator* qrng1 = new QuasiRandomNumberGenerator();
	QuasiRandomNumberGenerator* qrng2 = new QuasiRandomNumberGenerator();
	for (int i=0; i<20; i++) {
		std::cout << "qrng1->getQuasiRN() = " << qrng1->getQuasiRN() << std::endl;
		std::cout << "qrng2->getQuasiRN() = " << qrng2->getQuasiRN() << std::endl;
	}
	
	delete qrng1;
	delete qrng2;
	
	std::cout	<<"============= END OF TEST ================\n" <<std::endl;

	return 0;
}
