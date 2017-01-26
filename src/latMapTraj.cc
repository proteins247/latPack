
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include <biu/OptionParser.hh>
#include <biu/LatticeDescriptorSQR.hh>
#include <biu/LatticeDescriptorCUB.hh>
#include <biu/LatticeDescriptorFCC.hh>
#include <biu/LatticeDescriptorCKW.hh>
#include <biu/LatticeModel.hh>
#include <biu/Alphabet.hh>

#include <biu/LatticeProteinUtil.hh>
#include <biu/SuperPos_Kabsch.hh>

#include "version.hh"

/*!
 * Performs structural calculations on a latFold trajectory.
 * Based on latMap program
 * 
 * Calculations on latfold/latfoldvec output text files are made
 * relative to a reference structure. Measurements available are
 * 
 * - cRMSD
 * - GDT_TS
 * - GDT_HA
 * - native contacts
 * - non-native contacts
 * 
 * For the superpositioning the Kabsch algorithm is used.
 * 
 * @author (c) 2009 Martin Mann http://www.bioinf.uni-freiburg.de/~mmann/
 * @author victor zhao
 * 
 */


//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
// DEFINITIONS
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////

 /*! The parameter setup for this tools.
  * 
  * @param allowedArgs the provided parameters to fill
  * @param infoText the additional information text to setup 
  */
void
initAllowedArguments(biu::OptionMap & allowedArgs, std::string &infoText );

//////////////////////////////////////////////////////////////////////////

std::string
parseTrajLine(std::string const &line, const biu::Alphabet* const alphabet) {
	std::istringstream iss(line);
	int steps;
	float energy;
	std::string moveStr;
	iss >> steps >> energy >> moveStr;
	if (!alphabet->isAlphabetString(moveStr))
		moveStr.clear();
	return moveStr;
}

//////////////////////////////////////////////////////////////////////////

void
printPoints( const biu::DPointVec & p, std::ostream& out = std::cout ) {
	for (size_t i=0; i<p.size(); i++) {
		out <<" (" <<p.at(i) <<"),";
	}
	out <<std::endl;
}


//////////////////////////////////////////////////////////////////////////

void
printEvaluation( const biu::DPointVec & pos1
				, const biu::DPointVec & pos2
				, std::ostream& out = std::cout
				, const size_t outPrec = 3)
{
	out
		<<std::fixed <<std::setprecision(outPrec)
		<<"  " <<std::setw(6)<< biu::LatticeProteinUtil::cRMSD( pos1, pos2 )
		<<"  " <<biu::LatticeProteinUtil::GDT_TS( pos1, pos2 )
		<<"  " <<biu::LatticeProteinUtil::GDT_HA( pos1, pos2 );
		;

}

//////////////////////////////////////////////////////////////////////////

void
printContacts( const biu::IPointVec & ref
	       , const biu::IPointVec & pos
	       , biu::LatticeModel * lattice
	       , std::ostream& out = std::cout)
{
	int nativeContacts;
	int nonNativeContacts;
	biu::LatticeProteinUtil::countContacts(ref, pos, nativeContacts, nonNativeContacts, lattice);
	out     <<"  " <<std::setw(2) << nativeContacts
		<<"  " <<std::setw(2) << nonNativeContacts
		;

}



//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
// IMPLEMENTATIONS
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

int
main( int argc, char** argv ) 
{

	//////////////////////////////////////////////////////////////
	// parameter parsing and checking
	//////////////////////////////////////////////////////////////

	biu::OptionMap allowedArgs;
	std::string infoText;
	initAllowedArguments(allowedArgs,infoText);	// init

		// parse programm arguments
	biu::COptionParser opts = biu::COptionParser(	allowedArgs, argc, argv, infoText);
		// check if arguments parseable and all mandatory arguments given
	if (!opts.noErrors()) {
		return -1;
	}
		// help output
	if (opts.getBoolVal("help")) {
		opts.coutUsage();
		return 0;
	}
	  // version information requested
	if (opts.getBoolVal("version")) {
		giveVersion();
		return 0;
	}
	
	//////////////////////////////////////////////////////////////////////
	// setup parameters etc.
	//////////////////////////////////////////////////////////////////////
	
	// latMapTraj doesn't deal with side chains
	const bool sideChain = false;

	const double cAlphaDist = opts.getDoubleVal("cAdist");
	if (cAlphaDist <= 0.0) {
		std::cerr <<"\n   ERROR : C_alpha distance has to be greater than zero!\n";
		return -1;
	}
	// out precision fixed to 3
	const int outPrec = 3;
	// if (outPrec < 0) {
	// 	std::cerr <<"\n   ERROR : output precision has to be at least zero!\n";
	// 	return -1;
	// }
	
	// note, verbosity currently not used by latMapTraj
	enum OUTMODE { OUT_SILENT, OUT_NORMAL, OUT_VERBOSE };
	OUTMODE outMode = OUT_NORMAL;
	if (opts.argExist("s") && opts.argExist("v")) {
		std::cerr <<"\n   ERROR : cannot be silent and verbose at the same time!\n";
		return -1;
	}
	if (opts.argExist("s")) {
		outMode = OUT_SILENT;
	}
	if (opts.argExist("v")){
		outMode = OUT_VERBOSE;	
	}
	
	biu::LatticeDescriptor* latDescr = NULL;
	std::string latStr = opts.getStrVal("lat");
	if (latStr.compare("SQR") == 0)
		latDescr = new biu::LatticeDescriptorSQR();
	else if (latStr.compare("CUB") == 0)
		latDescr = new biu::LatticeDescriptorCUB();
	else if (latStr.compare("FCC") == 0)
		latDescr = new biu::LatticeDescriptorFCC();
	else if (latStr.compare("210") == 0)
		latDescr = new biu::LatticeDescriptorCKW();
	else {
		std::cerr <<"\n   ERROR : Unknown lattice type '"+latStr+"'\n\n";
		return -1;
	}
	// get alphabet
	const biu::Alphabet* const alphabet = latDescr->getAlphabet();
	  // create lattice
	biu::LatticeModel lattice(latDescr);
	
	  // get ref structure and trajetory file
	const std::string refIn = opts.getStrVal("ref");
	const std::string trajFilename = opts.getStrVal("traj");
	  // the according ref move sequences
	biu::MoveSequence refMoveSeq;

	const std::string whiteChars = " \t\n\r";
	size_t ref_begin = refIn.find_first_not_of( whiteChars , 0);
	size_t ref_end = refIn.find_last_not_of( whiteChars );
	if (ref_begin >= ref_end) {
		std::cerr <<"\n   ERROR : no first move string given!\n";
		return -1;
	}
	const std::string ref = refIn.substr(ref_begin, ref_end - ref_begin + 1);
	refMoveSeq = biu::LatticeProteinUtil::toMoveSequence(
							     ref
							     , *latDescr
							     , sideChain
							     );	
	// open I/O
	std::ostream* outstream = &std::cout; // could replace with file one day
	std::ifstream trajstream(trajFilename);
	if (!trajstream.is_open()) {
		std::cerr << "\n   ERROR : can not open trajectory file '" << trajFilename << "' !\n\n";
		return -1;
	}
	if (trajstream.bad()) {
		std::cerr << "\n   Error: cannot read trajectory file from '" << trajFilename << "' !\n";
		return -1;
	}

	//////////////////////////////////////////////////////////////////////
	// process the in file line by line...
	//////////////////////////////////////////////////////////////////////
	std::string line;
	std::string currMoveStr;
	biu::MoveSequence currMoveSeq;
	size_t moveStrLength;
	while (std::getline(trajstream, line)) {
		currMoveStr = parseTrajLine(line, alphabet);
		if (currMoveStr.empty()) {
			// wasn't a trajectory line; output the line unchanged
			*outstream << line << std::endl;
			continue;
		}

		*outstream << line;
		
		// now we do processing on the structure and append data
		moveStrLength = currMoveStr.length();
		if (moveStrLength > ref.length()) {
			*outstream << std::endl;
			continue;
		}

		*outstream << std::string(ref.length() - moveStrLength + 1, ' ');

		biu::MoveSequence::const_iterator first = refMoveSeq.begin();
		biu::MoveSequence::const_iterator last = refMoveSeq.begin() + moveStrLength;
		biu::MoveSequence subRefMoveSeq(first, last);
		currMoveSeq = biu::LatticeProteinUtil::toMoveSequence(
								      currMoveStr
								      , *latDescr
								      , sideChain
								      );
		biu::DPointVec refPos = biu::LatticeProteinUtil::toDblPoints(
									  subRefMoveSeq
									  , *latDescr
									  , sideChain
									  , cAlphaDist
									  );
		biu::DPointVec currPos = biu::LatticeProteinUtil::toDblPoints(
									      currMoveSeq
									      , *latDescr
									      , sideChain
									      , cAlphaDist
									      );
		biu::SuperPos_Kabsch::bestsuperposition( 
							refPos
							, currPos
							, latDescr->getAutomorphisms() );

		printEvaluation( refPos, currPos, *outstream, outPrec );

		biu::IPointVec refIntPos = biu::LatticeProteinUtil::toIntPoints(
										subRefMoveSeq
										, *latDescr
										, sideChain
										);
		biu::IPointVec currIntPos = biu::LatticeProteinUtil::toIntPoints(
										 currMoveSeq
										 , *latDescr
										 , sideChain
										 );
		printContacts( refIntPos, currIntPos, &lattice, *outstream);
		*outstream << std::endl;
	}

	//////////////////////////////////////////////////////////////////////
	// clear data structures
	//////////////////////////////////////////////////////////////////////
	
	delete latDescr;

	return 0;
	// current end
	// latMap code below

	// if (outMode >= OUT_VERBOSE) {
	// 	std::cout <<" abs1 : " 
	// 			<<biu::LatticeProteinUtil::toString( moves1, *latDescr, sideChain )
	// 			<<std::endl;
	// 	std::cout <<" abs2 : " 
	// 			<<biu::LatticeProteinUtil::toString( moves2, *latDescr, sideChain )
	// 			<<std::endl;
	// }


	// if (outMode >= OUT_VERBOSE) {
	// 	std::cout <<" pos1 : ";
	// 	printPoints(pos1, std::cout);
	// 	std::cout <<" pos2 : ";
	// 	printPoints(pos2, std::cout);
	// }
	
	// //////////////////////////////////////////////////////////////////////
	// // run superpositioning
	// //////////////////////////////////////////////////////////////////////
	
	// if (outMode >= OUT_VERBOSE) {
	// 	std::cout <<"\n ==> superpositioning :\n";
	// }
	
	// // TODO: try reflection ! maybe trigger via parameter !
	
	// biu::SuperPos_Kabsch::bestsuperposition( 
	// 								pos1
	// 								, pos2
	// 								, latDescr->getAutomorphisms() );
	
	// if (outMode >= OUT_VERBOSE) {
	// 	std::cout <<" sup1 : ";
	// 	printPoints(pos1, std::cout);
	// 	std::cout <<" sup2 : ";
	// 	printPoints(pos2, std::cout);
	// }
	
	// //////////////////////////////////////////////////////////////////////
	// // calculate distances
	// //////////////////////////////////////////////////////////////////////
	
	// if (outMode >= OUT_VERBOSE) {
	// 	std::cout <<"\n ==> distance :\n";
	// }
	//   // print distance evaluation of whole proteins
	// printEvaluation( pos1, pos2, std::cout, outPrec );
	
	// if (sideChain) {
	// 	std::cout <<"\n ==> backbone data only :\n";
	// 	  // get backbone data only
	// 	biu::DPointVec pos1bb(pos1.size());
	// 	biu::DPointVec pos2bb(pos2.size());
	// 	  // copy backbone data
	// 	for (size_t i=0; i<pos1.size(); i+=2) {
	// 		pos1bb[i/2] = pos1[i];
	// 		pos2bb[i/2] = pos2[i];
	// 	}
	// 	  // superposition backbone data
	// 	biu::SuperPos_Kabsch::bestsuperposition( 
	// 									pos1bb
	// 									, pos2bb
	// 									, latDescr->getAutomorphisms() );
	// 	  // print distance evaluation of protein backbones
	// 	printEvaluation( pos1bb, pos2bb, std::cout, outPrec );
	// }
	
	// if (sideChain) {
	// 	std::cout <<"\n ==> sidechain data only :\n";
	// 	  // get sidechain data only
	// 	biu::DPointVec pos1sc(pos1.size());
	// 	biu::DPointVec pos2sc(pos2.size());
	// 	  // copy sidechain data
	// 	for (size_t i=1; i<pos1.size(); i+=2) {
	// 		pos1sc[i/2] = pos1[i];
	// 		pos2sc[i/2] = pos2[i];
	// 	}
	// 	  // superposition sidechain data
	// 	biu::SuperPos_Kabsch::bestsuperposition( 
	// 									pos1sc
	// 									, pos2sc
	// 									, latDescr->getAutomorphisms() );
	// 	  // print distance evaluation of protein sidechains
	// 	printEvaluation( pos1sc, pos2sc, std::cout, outPrec );
	// }
	
}


void
initAllowedArguments(biu::OptionMap & allowedArgs, std::string &infoText )
{
	allowedArgs.push_back(biu::COption(
							"ref", false, biu::COption::STRING,
							"the first structure in absolute move string representation"));
	allowedArgs.push_back(biu::COption(
							"traj", false, biu::COption::STRING,
							"latFold or latFoldVec trajectory file (-out=S format)"));
	allowedArgs.push_back(biu::COption(
							"lat", false, biu::COption::STRING,
							"lattice of the structures : SQR - 2D square, CUB - 3D cubic, FCC - 3D face-centered-cubic, 210 - 3D chess knights walk)"));
	allowedArgs.push_back(biu::COption(
							"cAdist", true, biu::COption::DOUBLE,
							"the C_alpha atom distance to scale the lattice protein to", "3.8"));
	// allowedArgs.push_back(biu::COption(
	// 						"outPrec", true, biu::COption::INT,
	// 						"output precision, i.e. number of decimal places given", "3"));
	allowedArgs.push_back(biu::COption(
							"v", true, biu::COption::BOOL,
							"do verbose output (does nothing)"));
	allowedArgs.push_back(biu::COption(
							"s", true, biu::COption::BOOL,
							"silent (option does nothing right now)"));
	allowedArgs.push_back(biu::COption(
							"help", true, biu::COption::BOOL,
							"program parameters and help"));
	allowedArgs.push_back(biu::COption(
							"version", true, biu::COption::BOOL,
							"version information of this program"));

	infoText =	std::string(
"Performs structural calculations on a latFold trajectory.\n"
"Based on latMap program"
"\n\n"
"Calculations on latfold/latfoldvec output text files are made\n"
"relative to a reference structure. Measurements available are"
"\n\n"
"- cRMSD\n"
"- GDT_TS\n"
"- GDT_HA\n"
"- native contacts\n"
"- non-native contacts\n"
"\n\n"
"For the superpositioning the Kabsch algorithm is used.\n") ;

} // initArguments


//////////////////////////////////////////////////////////////////////////




