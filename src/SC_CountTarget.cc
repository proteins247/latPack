#include "SC_CountTarget.hh"

#include "ell/protein/S_LP.hh"

#include <biu/assertbiu.hh>

#include <cmath>

namespace ell
{

    SC_CountTarget::SC_CountTarget()
	: SC_Counting(), targetDefined(false), targetCount(0),
	  targetEnergy(0), absMoveStrings()
    {
    }

    SC_CountTarget::SC_CountTarget(size_t previousCount)
	: SC_Counting(previousCount), targetDefined(false), targetCount(0),
	  targetEnergy(0), absMoveStrings()
    {
    }

    SC_CountTarget::~SC_CountTarget()
    {
    }

    bool
    SC_CountTarget::hasTarget() const {
	return targetDefined;
    }

    // Define the target structure.
    void
    SC_CountTarget::defineTarget(const std::string& absMoves,
				 const biu::LatticeDescriptor& latDescr)
    {
	absMoveStrings.clear();
	// moveAlphabet.clear();

	targetDefined = true;
	targetCount = 0;

	// // gather all single move string representations of the lattice
	// const biu::Alphabet* const lpa = latDescr.getAlphabet();
	// for (size_t i=0; i<lpa->getAlphabetSize(); i++) {
	//     moveAlphabet += lpa->getString(lpa->getElement(i));
	// }
		
	// get all symmetric structures
	typedef std::set<biu::MoveSequence> MSset;
	MSset symmetric = latDescr.getAllSymmetricSequences(latDescr.getSequence(absMoves));
	// convert to move strings
	for (MSset::const_iterator it=symmetric.begin();it!=symmetric.end();it++) {
	    absMoveStrings.insert(latDescr.getString(*it));
	}

	moveStrSize = absMoveStrings.begin()->size();
	
    }

    // Define the target structure.
    void
    SC_CountTarget::setDegradationRate(double degradationScale)
    {
	trackSurvival = true;
	degradationRate = 1 / degradationScale;
    }

    // This function is used to track all added intermediate States.
    // @param s the added State
    void
    SC_CountTarget::add(const State& s) {
	//  call handler of superclass
	SC_Counting::add(s);

	if (targetDefined)
	{
	    assertbiu(dynamic_cast<const S_LP*>(&s)!=NULL,
		      "last state in state collector is no S_LP instance");
	    // get string representation
	    // We substring because s.toString() gives a string
	    //   that's more than just the move string.
	    std::string curMoves = s.toString().substr(0, moveStrSize);

	    // Code from WAC_Final, not needed
	    // // derive absolute move string
	    // curMoves = curMoves.substr( 0
	    // 				, curMoves.find_first_not_of(moveAlphabet));
		
	    // check whether or not this move string is present among the
	    // symmetric ones to check for ...
	    if (absMoveStrings.find(curMoves) != absMoveStrings.end())
	    {
		++targetCount;

		if (!targetEnergyFound)
		{
		    targetEnergy = s.getEnergy();
		}

		// Record first passage time if this is first time at target
		if (targetCount == 1)
		{
		    stepsToReachTarget = totalCount - 1;
		    // We subtract 1 because totalCount - 1 is the step number
		}

		if (trackSurvival)
		{
		    survivalSum += calculateSurvival();
		}
	    }
	
	}
	else
	{
	    // do nothing
	}
    }

    // access to target count
    size_t
    SC_CountTarget::getTargetCount() const {
	return targetCount;
    }

    // access to steps to reach target count
    size_t
    SC_CountTarget::getStepsToReachTarget() const {
	return stepsToReachTarget;
    }

    // access to energy of target conf
    double
    SC_CountTarget::getTargetEnergy() const {
	return targetEnergy;
    }

    // Access to target count as fraction of stateCount
    //   since the target was reached.
    double
    SC_CountTarget::getTargetCountFraction() const  {
	return targetCount / (double)(totalCount - stepsToReachTarget);
	// If target not reached, targetCount = stepsToReachTarget = 0
	// return value is 0.
    }

    // access survivalSum as fraction of total
    double
    SC_CountTarget::getSurvivalSumFraction() const  {
	return survivalSum / (double)stateCount;
    }

    // Calculate an extrapolated protein output.
    double
    SC_CountTarget::getExtrapolatedOutput() const  {
	if (stepsToReachTarget == 0)
	    return 0;
	double output = exp(
	    -degradationRate
	    * (stepsToReachTarget - (totalCount - stateCount)));
	double pnat = getTargetCountFraction();
	// In order to avoid divide by zero, pnat can't be 1
	if (pnat == 1)
	{
	    pnat = 1 - 1e-10;
	    // A double has ~15-16 digits of precision
	}
	output *= pnat / (degradationRate * (1 - pnat));
	return output;
    }

    // Calculate the probability of surviving till now.
    // (This is a protected function)
    double
    SC_CountTarget::calculateSurvival() const {
	double nonTargetCount = stateCount - targetCount;
	return exp(-degradationRate * nonTargetCount);
    }
       

}
