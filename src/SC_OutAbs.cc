#include "SC_OutAbs.hh"

#include <biu/assertbiu.hh>
#include <iomanip>

// the below three are for debugging by printing the points occupied
//   by the lattice protein
#include <biu/Point.hh>		// added(VZ)
#include <biu/LatticeProtein_I.hh> // added(VZ)
#include <ell/protein/S_LP.hh>	// added(VZ)

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
	
        // have to initialize out (it's type std::ostream&), but we will not use it
        //   when outputting to writer
        SC_OutAbs::SC_OutAbs(	HDF5TrajWriter* writer_, size_t cutoff_, size_t outFreq_, size_t previousCount)
	  :	SC_MinE(previousCount), out(std::cout), cutoff(cutoff_), outFreq(outFreq_), writer(writer_)
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
		
		// note, stateCount, totalCount starts at 1 (original structure, step 0)
		//   so subtract 1 to get step number
		
		if (writer && !((stateCount-1) % outFreq)) {
 			// extract structure information
		        std::string abs = (s.toString()).substr(0, cutoff);
			writer->write_buffered_traj(totalCount - 1, s.getEnergy(), &abs);
		}
		// else print to stream
		else if ( !((stateCount-1) % outFreq) ) {
			// extract structure information
			std::string abs = s.toString();
			abs = abs.substr(0, cutoff);
			out << std::setw(10) << totalCount - 1 << " "
			    << std::setw(6) << std::fixed << std::setprecision(2) << s.getEnergy()
			    << " " << abs << std::endl;
                        // below is added for debugging
                        // const S_LP* slp = dynamic_cast<const S_LP*>(&s);
                        // biu::LatticeProtein_I* latProt = slp->getProtein();
                        // biu::IPointVec points = latProt->getPoints();
                        // for (biu::IPointVec::iterator it = points.begin(); it!= points.end(); ++it) {
                        //          out << *it << std::endl;
                        // }

		}
	}

        void
	SC_OutAbs::outputLast() {
	        const State& s = *getLastAdded();

		if (writer) {
 			// extract structure information
		        std::string abs = (s.toString()).substr(0, cutoff);
			writer->write_buffered_traj(totalCount - 1, s.getEnergy(), &abs);
		} else {
			// extract structure information
			std::string abs = s.toString();
			abs = abs.substr(0, cutoff);
			out << std::setw(10) << totalCount - 1 << " "
			    << std::setw(6) << std::fixed << std::setprecision(2) << s.getEnergy()
			    << " " << abs << std::endl;
		}

	}

}
