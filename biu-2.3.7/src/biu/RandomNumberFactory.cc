// $Id: RandomNumberFactory.cc,v 1.2 2016/08/08 12:41:57 mmann Exp $


#include "biu/RandomNumberFactory.hh"
#include <biu/assertbiu.hh>

namespace biu
{

	// default random number generator
	RandomNumberGenerator* RandomNumberFactory::rng = new RNG_ISO();
	
	// pointer object to random generator function for RNF usage with STL
	// algorithms like stl::random_shuffle
	unsigned int (*RandomNumberFactory::pt_getRN)(unsigned int) = &RandomNumberFactory::getRN;


	RandomNumberFactory::RandomNumberFactory()
	{
	}
	
	RandomNumberFactory::~RandomNumberFactory()
	{
	}
	
	RandomNumberGenerator&
	RandomNumberFactory::getRNG() {
		return *rng;
	}
	
	void
	RandomNumberFactory::setRNG(RandomNumberGenerator& rng_) {

		delete rng; // delete old RNG
		
		rng = rng_.copy(); // set new one
	}
	
	void
	RandomNumberFactory::setRNG(RandomNumberGenerator* rng_) {

		assertbiu( rng_ != NULL, "given RNG is not available (NULL)");
		delete rng; // delete old RNG
		
		rng = rng_->copy(); // set new one
	}
	
	unsigned int
	RandomNumberFactory::getRN(void) {
		return rng->getRN();
	}
	
        unsigned int
        RandomNumberFactory::getRN(unsigned int max) {
                assertbiu(max != 0, "maximal value == 0 would cause division by zero");
                // The distribution of random numbers is slightly non-uniform
                // if max is not a divisor of maxRN, so we ignore numbers greater
                // than (maxRN/max)*max
                // Note: we might have to ignore up to halve the RNs drawn if max is
                // larger than getMaxRN/2
                unsigned int maxAllowedRn = max * (unsigned int)(getMaxRN() / max);
                // get an rn that is within allowed range
                unsigned int rn = getRN();
                while ( rn > maxAllowedRn )
                        rn = getRN();
                // compute final random number within [0,max)
                return rn%max;
        }
	
	unsigned int
	RandomNumberFactory::getMaxRN(void) {
		return rng->getMaxRN();
	}
	
} // namespace 
