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
	
        // have to initialize out (it's type std::ostream&), but we will not use it
        //   when outputting to writer
        SC_OutEnergy::SC_OutEnergy(	HDF5TrajWriter* writer_, size_t outFreq_, size_t previousCount)
	  :	SC_MinE(previousCount), out(std::cout), outFreq(outFreq_), writer(writer_)
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

		// note, stateCount, totalCount starts at 1 (original structure)
		//   subtract 1 to get step number
		
		if (writer && !((stateCount-1) % outFreq)) {
		        writer->write_buffered_traj(totalCount - 1, s.getEnergy(), nullptr);
		}
		// else print to stream
		else if  ( !((stateCount-1) % outFreq) ) {
			out << std::setw(10) << totalCount - 1 << " "
			    << std::setw(6) << std::fixed << std::setprecision(2) << s.getEnergy() << std::endl;
		}
	}

}
