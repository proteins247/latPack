#ifndef SC_COUNTTARGET_HH_
#define SC_COUNTTARGET_HH_

#include <set>
#include "ell/StateCollector.hh"
#include <biu/LatticeDescriptor.hh>

namespace ell
{

    /*! A StateCollector that counts the number of times that
     *! a target structure is added.
     * 
     * @author Victor Zhao
     */
    class SC_CountTarget : public SC_Counting
    {
    protected:
	//! whether or not this SC does anything
	bool targetDefined;

	int moveStrSize;
	size_t targetCount;
	size_t stepsToReachTarget = 0;
	double targetEnergy;

	//! the set of all absolute move strings symmetric to the given 
	//! structure
	std::set<std::string> absMoveStrings;
		
	// //! the letters of the current move alphabet to derive the move 
	// //! string from the state string representation
	// std::string moveAlphabet;
		
    public:
	
	SC_CountTarget();
	SC_CountTarget(size_t previousCount);
		
	virtual ~SC_CountTarget();
		
	virtual bool hasTarget() const;

	virtual void defineTarget(const std::string& absMoves,
				  const biu::LatticeDescriptor& latticeDescr);

	// This function is used to track all added intermediate States.
	// @param s the added State
	virtual void add(const State& s);
	
	//! access to target count
	virtual size_t getTargetCount() const;

	//! access to steps to reach target
	virtual size_t getStepsToReachTarget() const;

	virtual double getTargetEnergy() const;

	//! access to target count as fraction of count
	virtual double getTargetCountFraction() const;

    };

}

#endif /*SC_COUNTTARGET_HH_*/
