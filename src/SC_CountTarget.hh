#ifndef SC_COUNTTARGET_HH_
#define SC_COUNTTARGET_HH_

#include <set>
#include "ell/StateCollector.hh"
#include <biu/LatticeDescriptor.hh>

namespace ell
{

    /*! A StateCollector that counts the number of times that
     *! a target structure is added.
     *!
     *! The functions in this class help with calculating
     *! biophysically relevant properties.
     * 
     * @author Victor Zhao
     */
    class SC_CountTarget : public SC_Counting
    {
    protected:
	//! whether or not this SC does anything
	bool targetDefined;

	//! How long the target move string is
	int moveStrSize;

	//! How many times the target has been encountered.
	size_t targetCount;

	//! First passage time to target
	size_t stepsToReachTarget = 0;

	//! Energy of target conformation
	double targetEnergy;
	double targetEnergyFound = false;

	// The following variables have been added to track protein degradation
	bool trackSurvival = false;
	double degradationRate = 0;
	double survivalSum = 0;	

	//! the set of all absolute move strings symmetric to the given 
	//! structure
	std::set<std::string> absMoveStrings;
		
	// From WAC_Final, but not needed.
	// //! the letters of the current move alphabet to derive the move 
	// //! string from the state string representation
	// std::string moveAlphabet;

	//! Calculate the probability of surviving till now.
	virtual double calculateSurvival() const;
		
    public:
	
	SC_CountTarget();
	SC_CountTarget(size_t previousCount);
		
	virtual ~SC_CountTarget();
		
	virtual bool hasTarget() const;

	virtual void defineTarget(const std::string& absMoves,
				  const biu::LatticeDescriptor& latticeDescr);

	virtual void setDegradationRate(double degradationScale);

	// This function is used to track all added intermediate States.
	// @param s the added State
	virtual void add(const State& s);
	
	//! access to target count
	virtual size_t getTargetCount() const;

	//! access to steps to reach target
	virtual size_t getStepsToReachTarget() const;

	//! access to the energy of the target conformation
	virtual double getTargetEnergy() const;

	//! access to target count as fraction of stateCount (pnat)
	virtual double getTargetCountFraction() const;

	//! access survivalSum as fraction of total
	virtual double getSurvivalSumFraction() const;

	//! Calculate an extrapolated protein output.
	virtual double getExtrapolatedOutput() const;

    };

}

#endif /*SC_COUNTTARGET_HH_*/
