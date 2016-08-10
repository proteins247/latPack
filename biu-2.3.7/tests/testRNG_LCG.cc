#include "biu/RandomNumberGenerator.hh"
#include <iostream>
#include <string>

using namespace biu;

int main(int argc, char **argv) {
	
	std::cout	<<"========== BEGIN OF TEST ==================\n"
				<<" test function : RNG_LCG\n"
				<<"===========================================" <<std::endl;
	
	RandomNumberGenerator* rng = new RNG_LCG();
	
	unsigned int seed = 0;
	
	std::cout	<<"=> creating RNG_LCG()" << std::endl;
	rng = new RNG_LCG();
	std::cout	<<"-> getMaxRN() = " << rng->getMaxRN() << std::endl;
	std::cout	<<"=> testing generator" << std::endl;
	for (int i=0; i<5; i++) {
		std::cout	<<"-> getRN() = " << rng->getRN() << std::endl;
	}
	
	std::cout	<< "==> testing rng & seed" << std::endl;
		std::cout	<< "-> setSeed(1)" << std::endl;
		seed = 1;
		rng->setSeed(seed);
		for (int i=0; i<5; i++) {
			std::cout	<<"-> getRN() = " << rng->getRN() << std::endl;
	}
	
	std::cout	<< "==> testing rng & seed" << std::endl;
	std::cout	<< "-> setSeed(50)" << std::endl;
	seed = 50;
	rng->setSeed(seed);
	for (int i=0; i<5; i++) {
		std::cout	<<"-> getRN() = " << rng->getRN() << std::endl;
	}
	
	std::cout	<< "==> testing rng & seed" << std::endl;
	std::cout	<< "-> setSeed(50)" << std::endl;
	seed = 50;
	rng->setSeed(seed);
	for (int i=0; i<5; i++) {
		std::cout	<<"-> getRN() = " << rng->getRN() << std::endl;
	}

	std::cout	<< "==> testing rng & seed" << std::endl;
	std::cout	<< "-> setSeed(123456)" << std::endl;
	seed = 123456;
	rng->setSeed(seed);
	for (int i=0; i<5; i++) {
		std::cout	<<"-> getRN() = " << rng->getRN() << std::endl;
	}
	
	delete rng;
  	
  	std::cout	<<"============= END OF TEST =================\n" <<std::endl;
  	
	return 0; 
}
