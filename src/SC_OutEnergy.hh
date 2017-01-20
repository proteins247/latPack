#ifndef SC_OUTSTREAMENERGY_HH_
#define SC_OUTSTREAMENERGY_HH_

#include <iostream>
#include <string>

#include "SC_MinE.hh"

namespace ell
{

	/*! A WalkCollector that prints the energy of each reported state to the 
	 * given outstream.
	 * 
	 * @author Daniel Maticzka
	 */
	class SC_OutEnergy : public SC_MinE
	{
	protected:
		  //! the stream to write the states to
		std::ostream& out;
                const size_t outFreq;
		
	public:
	
		SC_OutEnergy(	std::ostream& out);
       	        SC_OutEnergy(	std::ostream& out, size_t outFreq, size_t previousCount);
		
		virtual ~SC_OutEnergy();
	
		  //! This function is used to track all accepted intermediate States.
		  //! @param s the accepted State
		virtual void add(const State& s);
			
	};

}

#endif /*SC_OUTSTREAMENERGY_HH_*/
