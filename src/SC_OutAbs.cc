#include "SC_OutAbs.hh"

#include <biu/assertbiu.hh>
#include <iomanip>

#include <biu/Point.hh>		// added
#include <biu/LatticeProtein_I.hh> // added
#include <ell/protein/S_LP.hh>	// added

namespace ell
{

	SC_OutAbs::SC_OutAbs(	std::ostream& out_, size_t cutoff_)
	  :	SC_MinE(), out(out_), cutoff(cutoff_), outFreq(1)
	{
	}
	
        SC_OutAbs::SC_OutAbs(	std::ostream& out_, size_t cutoff_, size_t outFreq_, size_t previousCount)
	  :	SC_MinE(previousCount), out(out_), cutoff(cutoff_), outFreq(outFreq_)
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
		  // stateCount starts at 1 (original structure)
		if ( !((stateCount-1) % outFreq) ) {
			// extract structure information
			std::string abs = s.toString();
			abs = abs.substr(0, cutoff);
			out << std::setw(10) << totalCount - 1 << " "
			    << std::setw(6) << std::setprecision(2) << s.getEnergy() << " "
			    << abs << std::endl;
                        // below is added for debugging
                        // const S_LP* slp = dynamic_cast<const S_LP*>(&s);
                        // biu::LatticeProtein_I* latProt = slp->getProtein();
                        // biu::IPointVec points = latProt->getPoints();
                        // for (biu::IPointVec::iterator it = points.begin(); it!= points.end(); ++it) {
                        //          out << *it << std::endl;
                        // }

		}
	}

}
