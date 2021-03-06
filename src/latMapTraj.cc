
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <csignal>
#include <cstring>

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
#include "HDF5_support.hh"

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

// additional function added for contact counting capability - VZ
void
printContacts( const biu::IPointVec & ref
	       , const biu::IPointVec & pos
	       , biu::LatticeModel & lattice
	       , std::ostream& out = std::cout)
{
	int nativeContacts;
	int nonNativeContacts;
	float fractionNativeContacts; // the fraction of nativeContacts out of total possible
	// (this is as opposed to fraction of total contacts that are native)
	biu::LatticeProteinUtil::countContacts(ref, pos, nativeContacts, nonNativeContacts, fractionNativeContacts, lattice);
	out     << std::fixed
	        <<"  " <<std::setw(2) << nativeContacts
		<<"  " <<std::setw(2) << nonNativeContacts
		<<"  " <<std::setw(5) << fractionNativeContacts
		;

}

//////////////////////////////////////////////////////////////////////////

// calculate radius of gyration
float
calcRadiusGyration(const biu::DPointVec & DPos)
{
	biu::DblPoint mean(0, 0, 0);
	for (const biu::DblPoint& point : DPos)
	{
		mean = mean + point;
	}
	mean /= (double)DPos.size();

	double radiusGyration = 0;
	for (const biu::DblPoint& point : DPos)
	{
		radiusGyration += pow(point.distance(mean), 2);
	}
	radiusGyration /= (double)DPos.size();

	return sqrt(radiusGyration);
}

//////////////////////////////////////////////////////////////////////////

// later function added to calculate all values in a single function, place into dict
void
evaluateStructure( const biu::DPointVec & refDPos,
		   const biu::DPointVec & DPos,
		   const biu::IPointVec & refIPos,
		   const biu::IPointVec & IPos,
		   biu::LatticeModel & lattice,
		   std::unordered_map<std::string, float> & data)
{
	float rmsd = biu::LatticeProteinUtil::cRMSD( refDPos, DPos );
	float gdt_ts = biu::LatticeProteinUtil::GDT_TS( refDPos, DPos );
	float gdt_ha = biu::LatticeProteinUtil::GDT_HA( refDPos, DPos );

	int nativeContacts;
	int nonNativeContacts;
	float fractionNativeContacts; // the fraction of nativeContacts out of total possible
	biu::LatticeProteinUtil::countContacts(refIPos, IPos, nativeContacts, nonNativeContacts, fractionNativeContacts, lattice);
	data["RMSD"] = rmsd;
	data["GDT_TS"] = gdt_ts;
	data["GDT_HA"] = gdt_ha;
	data["End-end Distance"] = DPos.front().distance(DPos.back());
	data["Radius of gyration"] = calcRadiusGyration(DPos);
	data["Native Contacts"] = (float)nativeContacts;
	data["Non-Native Contacts"] = (float)nonNativeContacts;
	data["Fraction Native"] = fractionNativeContacts;
}		   

//////////////////////////////////////////////////////////////////////////

// csignal for SIGINT, SIGTERM handling
volatile sig_atomic_t stopFlag = 0;

static void
handler(int signum)
{
	stopFlag = signum;
}


//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
// IMPLEMENTATIONS
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

int
main( int argc, char** argv ) 
{

	// setup interrupt catching
	struct sigaction sa;

	memset( &sa, 0, sizeof(sa) );
	sa.sa_handler = handler;
	sigemptyset(&sa.sa_mask);
	sa.sa_flags = SA_RESTART; /* Restart functions if
				     interrupted by handler */
	sigaction(SIGINT, &sa, NULL);
	sigaction(SIGTERM, &sa, NULL);

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
	
	// currently, latMapTraj doesn't deal with side chains
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
	
	  // get ref structure and trajectory file
	const std::string refIn = opts.getStrVal("ref");
	const std::string trajFilename = opts.getStrVal("traj");
	const bool hdf5Traj = opts.getBoolVal("hdf5");
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
	std::string currMoveStr;
	biu::MoveSequence currMoveSeq;
	size_t moveStrLength;

	if (!hdf5Traj) {
		// Text file trajectory
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
		while (std::getline(trajstream, line)) {
			currMoveStr = parseTrajLine(line, alphabet);
			if (currMoveStr.empty()) {
				// wasn't a trajectory line; output the line unchanged
				*outstream << line << std::endl;
				continue;
			}
			
			// print original line
			*outstream << line;
		
			// now we do processing on the structure and append data
			moveStrLength = currMoveStr.length();
			if (moveStrLength > ref.length()) {
				*outstream << std::endl;
				continue;
			}

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

			// print spaces
			*outstream << std::string(ref.length() - moveStrLength + 1, ' ');
			// print analysis measurements
			printEvaluation( refPos, currPos, *outstream, outPrec );
			printContacts( refIntPos, currIntPos, lattice, *outstream);
			*outstream << std::endl;
		}
	} else {
		// handle hdf5
		std::vector<std::string> names = {
			"RMSD",
			"GDT_TS",
			"GDT_HA",
			"End-end Distance",
			"Radius of gyration",
			"Native Contacts",
			"Non-Native Contacts",
			"Fraction Native"};
		std::unordered_map<std::string, float> data;
		for (auto &name : names)
			data.insert({name, 0.0});

		std::unique_ptr<HDF5TrajAnalyzer> analyzer;
		try
		{
			analyzer = std::unique_ptr<HDF5TrajAnalyzer>(
				new HDF5TrajAnalyzer(trajFilename.c_str(), &names));
		}
		catch (File_opening_error) {
			std::cerr << "Failed to open hdf5 file: "
				  << trajFilename << std::endl;
			exit(1);
		}

		size_t groupCount = analyzer->get_group_count();

		// go trajectory by trajectory
		for (size_t trajNum = 1; trajNum <= groupCount; trajNum++) {

			analyzer->open_trajectory_group(trajNum);

			while ( analyzer->read_structure_traj(&currMoveStr) >= 0 && !stopFlag ) {
				moveStrLength = currMoveStr.length();
				if (moveStrLength > ref.length()) {
					std::cerr << "Unexpected moveStrLength > ref length\n";
					raise(SIGINT);
					// we'll exit by SIGINT here to close the HDF5 file
				}

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
				evaluateStructure( refPos, currPos, refIntPos, currIntPos, lattice, data);
				analyzer->write_analysis(data);
				
			}

			if (stopFlag) {
				// explicitly call destructor because we're going to exit by SIGINT
				analyzer.reset(); 
				sa.sa_handler = SIG_DFL;
				sigaction(SIGINT, &sa, NULL);
				raise(SIGINT);
			}
			
			analyzer->close_trajectory_group();
		}
	}
	//////////////////////////////////////////////////////////////////////
	// clear data structures
	//////////////////////////////////////////////////////////////////////
	
	delete latDescr;

	return 0;
	
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
	allowedArgs.push_back(biu::COption(
							"hdf5", true, biu::COption::BOOL,
							"traj file is HDF5 format (as opposed to text)"));
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
"Based on latMap program. Input can either be text output of latFold\n"
"or latFoldVec or HDF5 file generated by those programs.\n"
"\n"
"Calculations on latfold/latfoldvec output text files are made\n"
"relative to a reference structure. Measurements performed include"
"\n\n"
"- RMSD\n"
"- GDT_TS\n"
"- GDT_HA\n"
"- End-end distance\n"
"- Radius of gyration\n"
"- Native contacts\n"
"- Non-native contacts\n"
"\n"
"Note that not all measurements are performed on the text file.\n"
"For the superpositioning the Kabsch algorithm is used.\n") ;

} // initArguments


//////////////////////////////////////////////////////////////////////////
