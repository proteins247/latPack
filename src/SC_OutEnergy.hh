#ifndef SC_OUTSTREAMENERGY_HH_
#define SC_OUTSTREAMENERGY_HH_

#include <iostream>
#include <string>

#include "SC_MinE.hh"
#include "HDF5_support.hh" // added(vz)

namespace ell
{

	/*! A WalkCollector that prints the energy of each reported state to the 
	 * given outstream.
	 * 
	 * This class was subsequently updated with facilities for writing
	 * to HDF5 file via HDF5_support.hh's HDF5TrajWriter class. Depending
	 * on constructor, the output format will differ. -VZ
	 * 
	 * @author Daniel Maticzka
	 */
	class SC_OutEnergy : public SC_MinE
	{
	protected:
		  //! the stream to write the states to
		std::ostream& out;
                const size_t outFreq;
		//! alternatively the handler class for HDF5 output
		//! can initialize here with c++11
		HDF5TrajWriter* writer = nullptr;

	public:
	
		SC_OutEnergy(	std::ostream& out);
       	        SC_OutEnergy(	std::ostream& out, size_t outFreq, size_t previousCount);
       	        SC_OutEnergy(	HDF5TrajWriter* writer, size_t outFreq, size_t previousCount);
		
		virtual ~SC_OutEnergy();
	
		  //! This function is used to track all accepted intermediate States.
		  //! @param s the accepted State
		virtual void add(const State& s);

		//! print the last state. 
		virtual void outputLast();
			
	};

}

#endif /*SC_OUTSTREAMENERGY_HH_*/
