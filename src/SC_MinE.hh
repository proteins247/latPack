#ifndef SC_MINE_HH_
#define SC_MINE_HH_


// #include "ell/StateCollector.hh"
#include "SC_CountTarget.hh"

enum OUT_MODE { OUT_NO, OUT_E, OUT_ES };

namespace ell
{

	/*! A StateCollector that stores the minimal energy seen so far.
	 * 
	 * @author Martin Mann
	 */
	class SC_MinE : public SC_CountTarget
	{
	protected:
		  //! the minimal energy seen so far
		double minE;
		
	public:
	
		SC_MinE();
		SC_MinE(size_t previousCount);
		
		virtual ~SC_MinE();
		
		virtual void add(const State& s);
	
		  //! access to the minimal energy seen so far
		virtual double getMinE() const;

		//! print the last state. For SC_MinE, the function is empty
		virtual void outputLast();
		
	};

}

#endif /*SC_MINE_HH_*/
