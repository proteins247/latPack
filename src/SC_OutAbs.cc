#include "SC_OutAbs.hh"

#include <biu/assertbiu.hh>

namespace ell
{

	SC_OutAbs::SC_OutAbs(	std::ostream& out_, size_t cutoff_, size_t outFreq_)
	  :	SC_MinE(), out(out_), cutoff(cutoff_), outFreq(outFreq_)
	{
	}
	
	SC_OutAbs::~SC_OutAbs()
	{
	}

	  // This function is used to track all added intermediate States.
	  // @param s the added State
	void 
	SC_OutAbs::add(const State& s) {
		  //  call handler of superclass
		SC_MinE::add(s);
		
		
		  // print to stream
		if ( !(stateCount % outFreq) ) {
			// extract structure information
			std::string abs = s.toString();
			abs = abs.substr(0, cutoff);
			out << s.getEnergy() << " " << abs << std::endl;
		}
	}

}
