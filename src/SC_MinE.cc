#include "SC_MinE.hh"

#include <algorithm>
#include <limits.h>

namespace ell
{

	SC_MinE::SC_MinE()
	 :	SC_CountTarget(), minE((double)INT_MAX)
	{
	}

	SC_MinE::SC_MinE(size_t previousCount)
	 :	SC_CountTarget(previousCount), minE((double)INT_MAX)
	{
	}

	SC_MinE::~SC_MinE()
	{
	}

	  // This function is used to track all added intermediate States.
	  // @param s the added State
	void
	SC_MinE::add(const State& s) {
		  //  call handler of superclass
		SC_CountTarget::add(s);

		  // check for minimal energy
		minE = std::min( s.getEnergy(), minE );
	}

	double
	SC_MinE::getMinE() const
	{
		return minE;
	}

	void
	SC_MinE::outputLast()
	{
	}


}
