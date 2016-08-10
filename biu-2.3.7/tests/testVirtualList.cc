/************************************************************
 *
 * Example usage of virtual lists
 *
 * As an example we define "States" of 01 strings, where the (virtual) list of
 * neighbors consists of all 01-strings in Hamming-distance 1.
 *
 * We show the definition of lists in "canonical" order and lists where
 * elements are iterated in random order.
 *
 *
 * @author Sebastian Will <will@@informatik.uni-freiburg.de>
 *
 * Copyright Sebastian Will, 2006
 *
 ************************************************************/

#include "biu/VirtualList.hh"

#include <iostream>
#include <string>
#include <vector>

// ------------------------------------------------------------
// CLASS State
//
//! here, a state is something that has neighbors and provides means to
//iterate through its neighbors either randomly or in a "canonical" order
//
// usually, we dont't want to create all neighbors before traversing them.
// Hence, we use virtual lists of the neighbors to avoid this (and preserve
// nice syntax at the same time).
//

class State {
public:
	typedef biu::VirtualList<State> NeighborList;

	virtual ~State();
	virtual NeighborList *getNeighborList() const = 0;
	virtual NeighborList *getRandomNeighborList() const = 0;

	//! return string representation
	virtual const std::string &strRep() const = 0;

	virtual State *clone() const = 0;
};

State::~State() {}

// ------------------------------------------------------------
// CLASS Hamming
//
//! implements structure and hamming move on structure. This class is used as
//  "merge in" in order to define HammingState
//

class Hamming {
protected:
	std::string structure;
public:
	Hamming(std::string structure_):structure(structure_) {};
	Hamming(const Hamming &h):structure(h.structure) {};
	virtual ~Hamming();

	virtual const Hamming &applyMove(int pos) {
		structure[pos] = (structure[pos]=='1')?'0':'1';
		return *this;
	}

	virtual const std::string& strRep() const {
		return structure;
	}

};
Hamming::~Hamming() {}



// ------------------------------------------------------------
// CLASS HammingState
//
//! defines State from Hamming, where Neighbors in distance
// of one Hamming move can be iterated

class HammingState: public State, protected Hamming {

public:
	HammingState(std::string structure_):Hamming(structure_) {};
	HammingState(const HammingState &state):Hamming(state.structure) {};

	virtual ~HammingState();

	virtual HammingState *clone() const {
		return new HammingState(*this);
	}

	// ------------------------------------------------------------
	// CLASS NeighborList
	//
	//! exemplifies how to derive virtual lists of neighbors
	//
	class NeighborList : public State::NeighborList {
	private:
		std::string origin;
	public:
		class ItState : public State::NeighborList::ItState {
		public:
			int pos;

			ItState():pos(0) {}

			ItState(const ItState &itstate) {pos=itstate.pos;}

			virtual ItState * clone() const {
				return new ItState(*this);
			}
		};

		NeighborList(std::string origin_):origin(origin_) {}

		virtual ~NeighborList();

		virtual State*
		next(State::NeighborList::ItState *itstate,
		     State* elem) const {

			int &pos = ((ItState*)itstate)->pos;

			pos++;

			if (pos>=(int) origin.length()) return 0;

			((HammingState*)elem)->structure=origin;
			((HammingState*)elem)->applyMove(pos);

			return elem;
		}

		virtual State*
		first(State::NeighborList::ItState** itstate) const {

			if (*itstate) delete *itstate;

			*itstate=new ItState();

			HammingState* s = new HammingState(origin);
			s->applyMove(0);

			return s;
		}
	};


	// ------------------------------------------------------------
	// CLASS RandomNeighborList
	//
	//! exemplifies how to derive virtual lists of neighbors
	// that are traversed in a *random* order
	//
	class RandomNeighborList:public State::NeighborList {
		std::string origin;
	public:
		class ItState : public State::NeighborList::ItState {
		public:
			std::vector<int> unseenPos;

			ItState(unsigned int n) {
				for (unsigned int i=0; i<n; i++)
					unseenPos.push_back(i);
			}

			ItState(const ItState &itstate) {
				unseenPos=itstate.unseenPos;
			}

			virtual ItState * clone() const {
				return new ItState(*this);
			}
		};

		RandomNeighborList(std::string origin_):origin(origin_) {}

		virtual ~RandomNeighborList();

		virtual State*
		next(State::NeighborList::ItState *itstate,
		     State* elem) const {

			std::vector<int> &unseenPos =
				((ItState*)itstate)->unseenPos;

			int n=unseenPos.size();

			if (n==0) return 0;

			  // random is system dependent
//			int k=random()%n;
			  // make test system independent
			int k = n/2;
			int pos=unseenPos[k];

			unseenPos[k]=unseenPos[n-1];
			unseenPos.resize(n-1);

			((HammingState*)elem)->structure=origin;
			((HammingState*)elem)->applyMove(pos);
			return elem;
		}

		virtual State*
		first(State::NeighborList::ItState** itstate) const {

			if (*itstate) delete *itstate;

			*itstate=new ItState(origin.length());
			HammingState* s = new HammingState(origin);

			return next(*itstate,s);
		}

	};

	// ------------------------------------------------------------
	// virtual methods of HammingState
	//

	virtual NeighborList *getNeighborList() const {
		return new NeighborList(structure);
	}

	virtual RandomNeighborList *getRandomNeighborList() const {
		return new RandomNeighborList(structure);
	}

	virtual const std::string& strRep() const {
		return Hamming::strRep();
	}

};

HammingState::~HammingState() {}
HammingState::NeighborList::~NeighborList() {}
HammingState::RandomNeighborList::~RandomNeighborList() {}


// nice type name for list of neighbors
typedef biu::VirtualList<State> NeighborList;

// ------------------------------------------------------------
// MAIN
//!

int main() {
	// srandom(time(0L));

	State *state = new HammingState("11111");

	NeighborList *neighbors = state->getRandomNeighborList();

	for (NeighborList::Iterator it = neighbors->begin();
	     neighbors->end() != it; ++it) {
		NeighborList *neighbors_it=it->getRandomNeighborList();
		for (NeighborList::Iterator it2=neighbors_it->begin();
		     neighbors_it->end() != it2; ++it2) {
			std::cout << it->strRep() << " "
					  << it2->strRep() << std::endl;
		}
	}

	delete neighbors;
	delete state;

	return 0;
}
