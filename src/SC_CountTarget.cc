#include "SC_CountTarget.hh"

#include "ell/protein/S_LP.hh"

#include <biu/assertbiu.hh>

// #include <algorithm>
// #include <limits.h>

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

    // 
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

    bool
    SC_CountTarget::hasTarget() const {
	return targetDefined;
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
	    std::string curMoves = s.toString().substr(0, moveStrSize);
	    // // derive absolute move string
	    // curMoves = curMoves.substr( 0
	    // 				, curMoves.find_first_not_of(moveAlphabet));
		
	    // check whether or not this move string is present among the
	    // symmetric ones to check for ...
	    if (absMoveStrings.find(curMoves) != absMoveStrings.end())
	    {
		++targetCount;
		targetEnergy = s.getEnergy();
		if (stepsToReachTarget == 0)
		{
		    stepsToReachTarget = totalCount;
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

    // access to target count as fraction
    double
    SC_CountTarget::getTargetCountFraction() const  {
	return targetCount / (double)stateCount;
    }


}
