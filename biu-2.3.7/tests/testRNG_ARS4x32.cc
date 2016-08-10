#include "biu/RandomNumberGenerator.hh"
#include <iostream>
#include <string>
#include <vector>

using namespace biu;

void runRNG( RandomNumberGenerator & rng, const unsigned int seed, const unsigned int testSize ) {

	std::cout	<< "==> testing rng & seed" << std::endl;
	std::cout	<< "-> setSeed("<<seed<<")" << std::endl;

	std::vector<unsigned int> randNumbers(testSize);

	// fill master set
	rng.setSeed(seed);
	for (int i=0; i<5; i++) {
		randNumbers[i] = rng.getRN();
	}
	// test master set
	rng.setSeed(seed);
	for (int i=0; i<5; i++) {
		if (randNumbers.at(i) != rng.getRN()) {
			std::cout <<"-> test failed since the "<<i<<". number of "<<testSize<<" is different when called with identic seed!!" <<std::endl;
			return;
		}
	}
	std::cout <<"-> test passed for " <<testSize<<" random numbers" <<std::endl;

}

int main(int argc, char **argv) {
	
	std::cout	<<"========== BEGIN OF TEST ==================\n"
				<<" test function : RNG_ARS4x32\n"
				<<"===========================================" <<std::endl;
	
	RandomNumberGenerator* rng = new RNG_ARS4x32();
	
	std::cout	<<"=> creating RNG_ARS4x32()" << std::endl;
	rng = new RNG_ARS4x32();
	std::cout	<<"-> getMaxRN() = " << rng->getMaxRN() << std::endl;

	// should return same as above
	runRNG( *rng, 1, 10 );
	
	// should return the same
	runRNG( *rng, 50, 10 );
	
	// shoul return something else
	runRNG( *rng, 123456, 10 );

	delete rng;
  	
  	std::cout	<<"============= END OF TEST =================\n" <<std::endl;
  	
	return 0; 
}
