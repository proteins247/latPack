#include "SC_OutEnergy.hh"

#include <iomanip>

namespace ell
{

	SC_OutEnergy::SC_OutEnergy(	std::ostream& out_)
	 :	SC_MinE(), out(out_), outFreq(1)
	{
	}
	
        SC_OutEnergy::SC_OutEnergy(	std::ostream& out_, size_t outFreq_, size_t previousCount)
	 :	SC_MinE(previousCount), out(out_), outFreq(outFreq_)
	{
	}
	
	SC_OutEnergy::~SC_OutEnergy()
	{
	}

	  // This function is used to track all added intermediate States.
	  // @param s the added State
	void 
	SC_OutEnergy::add(const State& s) {
		  //  call handler of superclass
		SC_MinE::add(s);

		  // print to stream
		  // stateCount starts at 1 (original structure)
		if  ( !((stateCount-1) % outFreq) ) 
			out << std::setw(10) << totalCount - 1 << " "
			    << std::setw(6) << std::setprecision(2) << s.getEnergy() << std::endl;
	}

}
