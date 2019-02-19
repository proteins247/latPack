
#include <assert.h>
#include <ell/protein/S_LP_PivotM.hh>
#include <ell/protein/S_LP_PullM.hh>
#include <ell/protein/WAC_LP.hh>
#include <ell/Walk.hh>
#include <ell/WalkAbortionCriterion.hh>
#include <limits.h>
#include <stddef.h>
#include <time.h>
#include <cmath>
#include <fstream>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iostream>
#include <iterator>
#include <set>
#include <set>
#include <string>
#include <utility>
#include <vector>
#include <csignal>
#include <cstring>

#include "HDF5_support.hh"

#include "../biu-2.3.7/src/biu/assertbiu.hh"
#include "../biu-2.3.7/src/biu/Alphabet.hh"
#include "../biu-2.3.7/src/biu/BioMolecule.hh"
#include "../biu-2.3.7/src/biu/DistanceEnergyFunction.hh"
#include "../biu-2.3.7/src/biu/LatticeDescriptorCUB.hh"
#include "../biu-2.3.7/src/biu/LatticeDescriptorFCC.hh"
#include "../biu-2.3.7/src/biu/LatticeDescriptorSQR.hh"
#include "../biu-2.3.7/src/biu/LatticeModel.hh"
#include "../biu-2.3.7/src/biu/LatticeNeighborhood.hh"
#include "../biu-2.3.7/src/biu/LatticeProtein_Ipnt.hh"
#include "../biu-2.3.7/src/biu/NeighborVector.hh"
#include "../biu-2.3.7/src/biu/OptionParser.hh"
#include "../biu-2.3.7/src/biu/PivotMoveSet.hh"
#include "../biu-2.3.7/src/biu/Point.hh"
#include "../biu-2.3.7/src/biu/PullMoveSet.hh"
#include "../biu-2.3.7/src/biu/RandomNumberFactory.hh"
#include "../biu-2.3.7/src/biu/RandomNumberGenerator.hh"
#include "../biu-2.3.7/src/biu/Timer.hh"
#include "energyFileSupport.hh"
#include "SC_OutAbs.hh"
#include "SC_OutEnergy.hh"
#include "version.hh"



////////////////////////////////////////////////////////////////////////////////

/**
 * Checks whether or not position with index idx of the given structure has a 
 * free neighbored position in the lattice.
 * 
 * @param str the structure in 3D coordinates to check
 * @param idx the index to check
 * @param lattice the lattice in use
 * @return true if the last position has at least one free neighbor, false otherwise
 */
bool
hasFreeNeighbor(	const biu::IPointVec & str,
					const size_t idx,
					const biu::LatticeModel & lattice)
{
	if (str.size() < 3)
		return true;
	
	assertbiu(idx < str.size(), "given index exceeds structure length");
	
	// check if the end position has a free neighbored position in the lattice
	
	biu::IntPoint actCheck = biu::IntPoint(0,0,0);

	// try all possible extensions
	for (	biu::LatticeNeighborhood::const_iterator actNeigh = lattice.getNeighborhood().begin(); 
			actNeigh != lattice.getNeighborhood().end(); actNeigh++ )
	{
		actCheck = str[idx] + *actNeigh;
		  // check if this neighbor is already part of the chain
		bool isNotFree = false;
		for (biu::IPointVec::const_iterator it = str.begin(); !isNotFree && it != str.end(); it++) {
			isNotFree = (*it == actCheck);
		}
		  // this neighbor is not part of the chain --> at least one free neighbor of idx !
		if (!isNotFree)
			return true;
	}
	
	return false;
}

////////////////////////////////////////////////////////////////////////////////


using namespace ell;

////////////////////////////////////////////////////////////////////////////////

//	/*!
//	 * This subclass of StateAcceptor accepts only S_LP states, where the last
//	 * coordinate position has a free neighbor in the lattice, i.e. is 
//	 * extensible during elongation
//	 * 
//	 * @author Martin Mann
//	 */	
//	class SA_FreeEnd : public StateAcceptor
//	{
//	protected:
//		
//		const StateAcceptor & SA_toForward;
//	public:
//		
//		SA_FreeEnd( const StateAcceptor& SA_forward )
//		 :	SA_toForward(SA_forward)
//		{}
//		
//		virtual ~SA_FreeEnd() {} 
//		
//		/*!
//		 * Always returns true.
//		 * @param sc the StateCollector to test against
//		 * @param succ the successor State to test
//		 * @return returns always true
//		 */
//		virtual bool accept(
//				const StateCollector* const sc, 
//				const State& succ) const
//		{
//			const S_LP* surLP = dynamic_cast<const S_LP*>(&succ);
//			assertbiu( surLP != NULL, "no S_LP state present");
//			
//		}
//	};

////////////////////////////////////////////////////////////////////////////////

// error values
static const int PARSE_ERROR = 1;
static const int DATA_ERROR = 2;
static const int IO_ERROR = 3;
// default data

static const bool optional = true;
static const std::string OPTION_CUB = "CUB";
static const std::string OPTION_SQR = "SQR";
static const std::string OPTION_FCC = "FCC";
static const std::string OPTION_PULLM = "PullM";
static const std::string OPTION_PIVOTM = "PivotM";

// default parameters
static const std::string DEFAULT_KT = "0.3";
static const std::string DEFAULT_SEED = "1";
static const double DEFAULT_MINE = (double)INT_MIN;
static const std::string DEFAULT_LATTICE = OPTION_CUB;
static const std::string DEFAULT_MOVES = OPTION_PULLM;
static const std::string DEFAULT_RUNS = "1";
static const std::string DEFAULT_OUTFREQ = "1";
static const std::string DEFAULT_ELEMENTLENGTH = "1";

// constants
static const double DELTA_E_ADD = 0.0001;

 // possible output modes
// enum OUT_MODE { OUT_NO, OUT_E, OUT_ES };
// definition moved to SC_MinE.hh - VZ

// infotexts
static const std::string infotext = 
	"LatFoldVec: co-translational protein folding simulation on different"
	" lattices using the specified move set and a metropolis walk within"
	" chain elongation steps.\n"
	"\n"
	"The minimum protein length that will be simulated is 5."
	"\n"
	"The energy function is contact based and has to be"
	" provided as a text file that contains the alphabet"
	" and the energy contributions in matrix form.\n"
	"\n"
	"An example energy file : \n"
	"------------------------------\n"
	"HPNX\n"
	"-4.0  0.0  0.0  0.0\n"
	" 0.0 +1.0 -1.0  0.0\n"
	" 0.0 -1.0 +1.0  0.0\n"
	" 0.0  0.0  0.0  0.0\n"
	"------------------------------\n"
	"\n" 
	;
	
static const std::string seqInfo =
	"protein sequence (sequence valid for alphabet from energy file). length >= 4 (although length = 4 will not produce MC dynamics)";
static const std::string titleInfo =
	"Give simulations a title/label; most relevant for HDF5 output";
static const std::string ktInfo =
	"kT parameter for metropolis criterion "
	"(double value from [0, infinity))";
static const std::string maxStepsInfo =
	"each intermediate simulation ends after [maxSteps] steps. It's called maxSteps because simulation can be interrupted by minE or final criteria";
static const std::string elongationInfo =
	"a custom elongation schedule, given as a sequence of L - 4 integers, "
	"where L is the length of the protein sequence. example: if the sequence "
	"(-seq) has length of 8, \"-elongationSchedule=50 55 60 50\" is valid. "
	"This option overrides -maxSteps. Note that the i'th number in "
	"-elongationSchedule corresponds to the number of simulations steps when "
	"i+4 residues are present.\nTo wit, the 1st number corresponds to the "
	"number of steps when 5 residues are present, and the last number "
	"corresponds to the number of steps when all residues are present. "
	"Each number therefore corresponds to the translation time for the "
	"(n+5)th residue.";
static const std::string maxStepsIncreaseInfo =
	"length-proportional MC simulation. Either the custom elongation schedule or maxSteps is multiplied by sequence length to obtain the number of MC steps for each protein length";
static const std::string minenergyInfo =
	"walk ends if energy gets below or equal to [minE]";
static const std::string finalInfo =
	"if present, each folding simulation run is aborted if the given "
	"(or a symmetric) structure is visited.";
static const std::string targetInfo =
	"if present: count number of times folding simulation visits the given"
	" (or symmetric) structure.";
static const std::string fullLengthStepsInfo =
	"after full elongation, sim will run additional number of steps (or until final structure is reached if -final is present). By default, the number of steps run on the final structure is given by -maxSteps or last number in -elongationSchedule. If -ribosomeRelease is present, chain will be free.";
static const std::string seedInfo =
	"seed for random number generator. should be <= 2^31 - 1 "
	"[uses biu::RNG_ARS4x32, a counter-based generator from the Random123 library]";
static const std::string runsInfo =
	"number of folding simulations to perform";
static const std::string latticeInfo =
	"which lattice to use: CUB, SQR or FCC";
static const std::string ofileInfo =
	"write output of simulations to filename (HDF5 format). if argument is 'STDOUT' or this argument is not present, output is written to standard output (plain text).";
static const std::string timingInfo =
	"print cpu-time used";
static const std::string verbosityInfo =
	"be verbose";
static const std::string vvInfo =
	"be extra verbose";
static const std::string helpInfo =
	"display program parameters and help";
static const std::string moveSetInfo =
	"which move set to use: PullM or PivotM";
static const std::string ribosomeInfo =
	"simulate with a non-interacting wall that acts as a barrier and anchors the last residue. follows convention that synthesis is from N term to C term";
static const std::string ribosomeReleaseInfo =
	"simulate as unbound from ribosome (no constraints) after last residue added. (requires -fullLengthSteps)";

// csignal for SIGINT, SIGTERM handling
volatile sig_atomic_t stopFlag = 0;

static void handler(int signum)
{
	stopFlag = signum;
}


int main(int argc, char** argv) {

	// If we receive a SIGINT or SIGTERM, stopFlag will change to nonzero
	//  but we will continue execution until a point where we can handle it.
	struct sigaction sa;

	memset( &sa, 0, sizeof(sa) );
	sa.sa_handler = handler;
	sigemptyset(&sa.sa_mask);
	sa.sa_flags = SA_RESTART; /* Restart functions if
				     interrupted by handler */
	sigaction(SIGINT, &sa, NULL);
	sigaction(SIGTERM, &sa, NULL);
	
	/*
	 * parse input
	 */
	biu::OptionMap options;
	
	options.push_back(biu::COption(
			"seq", !optional, biu::COption::STRING, seqInfo));
	
	options.push_back(biu::COption(
			"title", optional, biu::COption::STRING, titleInfo));

	options.push_back(biu::COption(	
			"energyFile", !optional, biu::COption::STRING, 
			"the contact energy function file that contains also the alphabet"));

	options.push_back(biu::COption(	
			"elementLength", optional, biu::COption::INT,
			"the character length of a single alphabet element in the energy file",
			DEFAULT_ELEMENTLENGTH));

	options.push_back(biu::COption(	
			"energyForDist",  true, biu::COption::BOOL, 
			"if present a distance based energy function is used, otherwise a contact energy function is applied"));
	
	options.push_back(biu::COption(	
			"energyCalphaDist",  true, biu::COption::DOUBLE, 
			"if '-energyForDist' is present, this value is used to scale the C_alpha distances of the given energy function to the C_alpha monomer distances in the used lattice model",
			"3.8"));

	options.push_back(biu::COption(
			"kT",  optional, biu::COption::DOUBLE, ktInfo, DEFAULT_KT));
	
	options.push_back(biu::COption(
			"maxSteps", optional, biu::COption::INT, maxStepsInfo, "20"));

	options.push_back(biu::COption(
			"elongationSchedule", optional, biu::COption::STRING, elongationInfo));
	
	options.push_back(biu::COption(
			"maxStepsIncrease", optional, biu::COption::BOOL, maxStepsIncreaseInfo));
	
	options.push_back(biu::COption(
			"fullLengthSteps", optional, biu::COption::INT, fullLengthStepsInfo));
		
	options.push_back(biu::COption(
			"minE", 
			optional,
			biu::COption::DOUBLE,
			minenergyInfo ));
	
	options.push_back(biu::COption(
			"final", optional, biu::COption::STRING, finalInfo));
		
	options.push_back(biu::COption(
			"countTarget", optional, biu::COption::STRING, targetInfo));
		
	options.push_back(biu::COption(
			"seed", optional, biu::COption::INT, seedInfo, 
			DEFAULT_SEED));
	
	options.push_back(biu::COption(
			"runs", optional, biu::COption::INT, runsInfo, DEFAULT_RUNS));
	
	options.push_back(biu::COption(
			"lat", optional, biu::COption::STRING, latticeInfo, DEFAULT_LATTICE));
	
	options.push_back(biu::COption(
			"moveSet", optional, biu::COption::STRING, moveSetInfo, DEFAULT_MOVES));
	
	options.push_back(biu::COption(
			"ribosome", optional, biu::COption::BOOL, ribosomeInfo));
	
	options.push_back(biu::COption(
			"ribosomeRelease", optional, biu::COption::BOOL, ribosomeReleaseInfo));
	
	options.push_back(biu::COption(
			"out", optional, biu::COption::CHAR,
			"output mode along the folding simulation: (N)o, (E)nergy, (S)tructure+Energy",
			"N"));
	
	options.push_back(biu::COption(
			"outFreq", optional, biu::COption::INT,
			"output frequency. print information every INT steps. has effect only if -out=E or -out=S",
			DEFAULT_OUTFREQ));

	options.push_back(biu::COption(
			"outFile", optional, biu::COption::STRING, ofileInfo, "STDOUT"));
	
	options.push_back(biu::COption(
			"outTiming", optional, biu::COption::BOOL, timingInfo));
	
	options.push_back(biu::COption(
			"v", optional, biu::COption::BOOL, verbosityInfo));
	
	options.push_back(biu::COption(
			"vv", optional, biu::COption::BOOL, vvInfo));
	
	options.push_back(biu::COption(
			"help", optional, biu::COption::BOOL, helpInfo));
	options.push_back(biu::COption(
			"version", true, biu::COption::BOOL,
			"version information of this program"));


	biu::COptionParser parser(options, argc, argv, infotext);
	
	// values of those depend on command line arguments
	biu::Alphabet * alph = NULL;
	biu::EnergyMatrix * energyMatrix = NULL;
	biu::DistanceEnergyFunction * energy = NULL;
	std::string title;
	std::string seqStr;
	std::string absMoveStr;
	std::string absMoveStrFinal;
	std::string targetMoveStr;
	std::string moves;
	double kT;
	std::vector<unsigned int> elongationSchedule;
	bool sequenceDependentSimLength = false;
	double minEnergy;
	biu::LatticeDescriptor* latticeDescriptor = NULL;
	biu::LatticeModel * lattice = NULL;
	size_t seed;
	int runs;
	unsigned int outFreq;
	unsigned int alphElementLength;
	std::ostream* outstream = &std::cout;
	ell::SC_MinE* sc;
	double minE;
	double targetCountFraction;
	double targetConfEnergy;
	size_t stepsToReachTarget;
	bool ribosome;
	bool ribosomeRelease;
	bool timing;
	size_t verbosity;
	OUT_MODE simOutMode = OUT_NO;
	bool outHDF = false;
	std::unique_ptr<HDF5TrajWriter> hdf5writer;	
	// std::ostream* simOut = &std::cout; // will remove this
	std::string alphString = ""; // temporary alphabet string representation
	double cAlphaDist = 3.8;	//!< the C_alpha distance used to scale the distances of the energy function
	
	if (parser.noErrors()) {
		
		if (parser.argExist("help")) {
			parser.coutUsage();
			return 0;
		}
		if (parser.argExist("version")) {
			giveVersion();
			return 0;
		}

		if (parser.argExist("title"))
			title = parser.getStrVal("title");
		
		if (parser.argExist("lat"))
		{
			std::string lattice = parser.getStrVal("lat");
			if (lattice == OPTION_CUB)
			{
				latticeDescriptor = new biu::LatticeDescriptorCUB();
			}
			else if (lattice == OPTION_SQR)
			{
				latticeDescriptor = new biu::LatticeDescriptorSQR();
			}
			else if (lattice == OPTION_FCC)
			{
				latticeDescriptor = new biu::LatticeDescriptorFCC();
			}
			else
			{
				std::cerr 	<< "Error: lattice must be one of the following: "
							<< OPTION_CUB << ", "
							<< OPTION_SQR << ", "
							<< OPTION_FCC << ". Is: "
							<< lattice << std::endl;
				return PARSE_ERROR;
			}
		}
		
		/*
		 * Building Lattice related objects.
		 */
		lattice = new biu::LatticeModel(latticeDescriptor);
		
		
			// init energy function
		std::string energyFile = parser.getStrVal("energyFile");
		alphElementLength = parser.getIntVal("elementLength");
		if ( energyFile.size() == 0 ) {
			std::cerr <<"\n   Error: no energy file given ('-energyFile=XXX') !\n";
			return -1;
		}
		{
			  // temporary data structures
			std::ifstream *inFile = NULL;
			std::istream* in = &std::cin;
			  // open stream if file given
			if (energyFile.compare("STDIN") != 0) {
				inFile = new std::ifstream( energyFile.c_str() );
				if (!inFile->is_open()) {
					std::cerr <<"\n   ERROR : can not open energy file '"+energyFile+"' !\n\n";
					return -1;
				}
				in = inFile;
			}
			if (in->bad()) {
				std::cerr <<"\n   Error: cannot read energy function from '" <<energyFile <<"' !\n";
				return -1;
			}
			if (parser.argExist("energyForDist")) {
			/////// DISTANCE BASED ENERGY FUNCTION ///////////////////////////////////
		
				double cAlphaDistScale = (lattice->getNeighborhood().getElementByIndex(0).distance(biu::IntPoint(0,0,0)))
											/ cAlphaDist;
				
				// do parsing
				if (initIntervalEnergyFunction( alph, energy, cAlphaDistScale, *in, alphElementLength) != 0) {
					std::cerr <<"\n   Error: the given energy file '"<<energyFile <<"' is not valid !\n";
					return -1;
				}
				
			} else {
			/////// CONTACT BASED ENERGY FUNCTION ////////////////////////////////////
				  // do parsing
				if (initContactEnergyFunction( alph, energyMatrix, *in, alphElementLength) != 0) {
					std::cerr <<"\n   Error: the given energy file '"<<energyFile <<"' is not valid !\n";
					return -1;
				}
				  // create energy function
				energy = new biu::ContactEnergyFunction(alph, energyMatrix, lattice);
			}
			  // close stream if necessary
			if (energyFile.compare("STDIN") != 0) {
				in = &std::cin;
				inFile->close();
				delete inFile;
			}
		}
		  // init string representation of the used alphabet for parameter output
		{
			biu::Sequence tmp;
			for (size_t i = 0; i < alph->getAlphabetSize(); i++) {
				tmp.push_back(alph->getElement(i));
			}
			alphString = alph->getString(tmp);
		}

		  // parse sequence
		if (parser.argExist("seq")) 
		{
		    seqStr = parser.getStrVal("seq");
		    // check about alphabet
		    if ( !alph->isAlphabetString(seqStr)) {
		    	std::cerr	<<"Error: given sequence '" 
		    				<<seqStr 
		    				<<"' is no valid sequence for alphabet '" 
		    				<<alphString 
		    				<<"' !\n";
		    	return PARSE_ERROR;
		    }
		    // check size
		    if (seqStr.size()/alphElementLength < 4) {
		    	std::cerr	<< "Error: seq ("
					<< seqStr
					<< ") must have at least length 4."
		    				<< std::endl;
		    	return PARSE_ERROR;
		    }
		}
		
		
		if (parser.argExist("kT")) 
		{
			kT = parser.getDoubleVal("kT");
			if (kT <= 0) {
				std::cerr << "Error: kT must be > 0." << std::endl;
				return PARSE_ERROR;
			}
		}
		
		if (parser.argExist("maxSteps")) // Note, an argument always exists if it has a default value
		{
			// By default, we run the same number of steps for each elongation.
			// So we fill our elongationSchedule vector with the same value
			int maxSteps = parser.getIntVal("maxSteps"); 
			if (maxSteps <= 0) {
				std::cerr << "Error: maxSteps must be > 0" << std::endl;
				return PARSE_ERROR;
			}
			// makes elongationSchedule size seqStr.size()
			elongationSchedule = std::vector<unsigned int>(seqStr.size(), (unsigned int)maxSteps);
		}

		if (parser.argExist("elongationSchedule"))
		{
			// But if the user gave a custom schedule,
			// parse the elongationSchedule string into a vector of uints
			std::istringstream scheduleReader(parser.getStrVal("elongationSchedule"));
			std::string buffer;
			std::vector<unsigned int> schedule;
			unsigned int value;
			while (scheduleReader >> buffer) {
				try {
					value = (unsigned int) std::stoi(buffer);
				} catch (std::invalid_argument) {
					std::cerr << "Error: In reading elongation schedule, '"
						  << buffer
						  << "' could not be converted to int."
						  << std::endl;
					return PARSE_ERROR;
				}
				schedule.push_back(value);
			}
			if (schedule.size() != seqStr.size() - 4) {
				std::cerr << "Error: -elongationSchedule requires sequence length - 4 (="
					  << seqStr.size() - 4
					  << ") ints"
					  << std::endl;
				return PARSE_ERROR;
			}
			std::copy(schedule.begin(), schedule.end(), elongationSchedule.begin() + 4);
			// i.e. copy `schedule` to elongationSchedule[4:]
			// Note the offset, since later, we access using: elongationSchedule[currentLength]
		}
		
		sequenceDependentSimLength = parser.argExist("maxStepsIncrease");
		
		if (parser.argExist("minE"))
		{
			minEnergy = parser.getDoubleVal("minE");
		} else {
			minEnergy = (double)INT_MIN;
		}
		
		if (parser.argExist("final"))
		{
			absMoveStrFinal = parser.getStrVal("final");
			if (absMoveStrFinal.size() != seqStr.size() - 1) {
				std::cerr << "Warning: given '-final' move string, "
					  << absMoveStrFinal
					  << ", doesn't correspond to given protein: "
					  << seqStr
					  << std::endl;
			}
			if (!latticeDescriptor->getAlphabet()->isAlphabetString(absMoveStrFinal)) {
				std::string moveAlphStr = latticeDescriptor->getAlphabet()->getString(latticeDescriptor->getAlphabet()->getElement(0));
				for (size_t i=1; i<latticeDescriptor->getAlphabet()->getAlphabetSize();i++) {
					moveAlphStr += ",";
					moveAlphStr += latticeDescriptor->getAlphabet()->getString(latticeDescriptor->getAlphabet()->getElement(i));
				}
				std::cerr	<<"Error: given '-final' move string '"
						<<absMoveStrFinal
						<<"' is not valid for the absolute move alphabet {"
						<<moveAlphStr
						<<"} !\n";
				return PARSE_ERROR;
			}
		}

		if (parser.argExist("countTarget")) {
			targetMoveStr = parser.getStrVal("countTarget");

			if (targetMoveStr.size() != seqStr.size() - 1) {
				std::cerr << "Warning: given '-countTarget' move string, "
					  << targetMoveStr
					  << ", doesn't correspond to given protein: "
					  << seqStr
					  << std::endl;
			}
			if (!latticeDescriptor->getAlphabet()->isAlphabetString(targetMoveStr)) {
				std::string moveAlphStr = latticeDescriptor->getAlphabet()->getString(latticeDescriptor->getAlphabet()->getElement(0));
				for (size_t i=1; i<latticeDescriptor->getAlphabet()->getAlphabetSize();i++) {
					moveAlphStr += ",";
					moveAlphStr += latticeDescriptor->getAlphabet()->getString(latticeDescriptor->getAlphabet()->getElement(i));
				}
				std::cerr	<<"Error: given '-countTarget' move string '"
						<<targetMoveStr
						<<"' is not valid for the absolute move alphabet {"
						<<moveAlphStr
						<<"} !\n";
				return PARSE_ERROR;
			}
		}

		if (parser.argExist("fullLengthSteps")) {
			unsigned int fullLengthSteps = parser.getIntVal("fullLengthSteps");
			elongationSchedule.push_back(fullLengthSteps);
			// lengthens elongationSchedule by 1
		}

		if (parser.argExist("seed")) 
		{
			seed = parser.getIntVal("seed");
		}
		
		if (parser.argExist("runs"))
		{
			runs = parser.getIntVal("runs");
			if (runs < 1) {
				std::cerr << "Error: runs must be at least 1." << std::endl;
				return PARSE_ERROR;
			}
			
		}
		
		if (parser.argExist("moveSet"))
		{
			moves = parser.getStrVal("moveSet");
			if (moves != OPTION_PIVOTM && moves != OPTION_PULLM) {
				std::cerr 	<< "Error: moveSet must be one of the following: "
							<< OPTION_PULLM << ", "
							<< OPTION_PIVOTM << ". Is: "
							<< moves << std::endl;
				return PARSE_ERROR;
			}
		}
		
		// ribosome options
		ribosome = parser.argExist("ribosome");
		ribosomeRelease = parser.argExist("ribosomeRelease");

		  // check for simulation output mode
		switch (parser.getCharVal("out")) {
		case 'N' :
			simOutMode = OUT_NO;
			break;
		case 'E' : 
			simOutMode = OUT_E;
			break;
		case 'S' :
			simOutMode = OUT_ES;
			break;
		default  :
			std::cerr	<<"Error: given output mode '"
						<<parser.getCharVal("out")
						<<"' is not supported !\n";
			return PARSE_ERROR;
		}
		  // setup HDF5 output if outFile argument is not STDOUT
		if (	simOutMode != OUT_NO 
				&& parser.argExist("outFile") 
				&& parser.getStrVal("outFile").compare("STDOUT") != 0)
		{
			outHDF = true;
			std::string filename = parser.getStrVal("outFile");
			try {
				hdf5writer = std::unique_ptr<HDF5TrajWriter>(new HDF5TrajWriter(filename.c_str(), simOutMode,
												seqStr.size()-1));
				// seqStr.size()-1 is the size of the move string for the structure
			}
			catch (File_exists_error) {
				std::cerr << "Could not create output file at "
					  << filename << std::endl;
				return IO_ERROR;
			}
		} 
		  // check if at least some output hast to be produced
		if ( simOutMode == OUT_NO && !parser.argExist("minE") && !parser.argExist("final")) {
			std::cerr <<"\nError: neither simulation output \n\t NOR minimal energy \n\t NOR final structure given !\n\t What to compute ?\n";
			return PARSE_ERROR;
		}
		
		 // set output frequency
		if (	simOutMode != OUT_NO
			        && parser.argExist("outFreq")   ) {
	        	outFreq = parser.getIntVal("outFreq");
			if (outFreq <= 0)
				std::cerr << "Error: given output frequencey '"
					  << outFreq
					  << "' is not supported. Require outFreq > 0\n";
		}

		verbosity = 0;
		if (parser.argExist("v")) verbosity = 1;
		if (parser.argExist("vv")) verbosity = 2;
		
		timing = parser.argExist("outTiming");
		
	}
	else
	{
		return PARSE_ERROR;
	}
	
	biu::PullMoveSet::PullMoveDecoder pmd(lattice);

	/*
	 * Preparing Randomizer
	 */
	biu::RandomNumberGenerator* rng = new biu::RNG_ARS4x32();
	biu::RNF::setRNG( rng );
	delete rng;
	// To avoid using the same RN sequence for all simulations with the same user seed
	//   (but different sequences and/or temperature), a modified seed is generated
	//   based on the user seed, as well as temperature and sequence.
	// See Sindhikara et al. JCTC 2009, 5
	// "Bad Seeds Sprout Perilous Dynamics: Stochastic Thermostat Induced Trajectory Synchronization in Biomolecules"
	// Robustness of multiplication and addition after hash uncertain.
	unsigned int modified_seed = 0;
	for (std::string::iterator it=seqStr.begin(); it!=seqStr.end(); ++it) {
		// (sdbm hash)
		modified_seed = *it + (modified_seed << 6) + (modified_seed << 16) - modified_seed;
	}
	modified_seed *= (unsigned int)(kT * 100);
	biu::RNF::getRNG().setSeed(modified_seed + seed);
	
	  // output parameter setting
	if (verbosity > 0) {
		if (parser.argExist("title"))
			*outstream << "Title : " << title;
		*outstream	<< "\n Parameter setup :"
					<< "\n ================="
					<< "\n  - Version     : " << BIN_PACKAGE_VERSION
					<< "\n  - Lattice     : " << lattice->getDescriptor()->getName()
					<< "\n  - Energy file : " << parser.getStrVal("energyFile")
					<< "\n  - Alphabet    : " << alphString
					<< "\n  - Elt. length : " << alphElementLength
					<< "\n  - Sequence    : " << seqStr
					<< "\n  - Move set    : " << moves
         				<< "\n  - Ribosome?   : " << (ribosome ? "tethered" : "untethered")
         				<< "\n  - Ribo release: " << (ribosomeRelease ? "yes" : "no")
					<< "\n  - Simulations : " << runs
             				<< "\n  - Seed (rand) : " << seed << " (" << modified_seed + seed << ")"
					<< "\n  - kT (MC)     : " << kT;
		if (parser.argExist("elongationSchedule")) {
			*outstream << "\n  - Elongation  : ";
			std::ostringstream schedule;
			for (unsigned int i = 4; i < seqStr.size(); ++i) {
				schedule << elongationSchedule[i] << " ";
			}
			*outstream << schedule.str();
		}
		else {
			*outstream << "\n  - Max. steps  : " << elongationSchedule[0];
		}
		if (sequenceDependentSimLength)
			*outstream << "\n  + steps incr. : (steps * seqLength)";
		if (parser.argExist("fullLengthSteps"))
			*outstream << "\n  + Full steps  : " << parser.getIntVal("fullLengthSteps");
		if (minEnergy != DEFAULT_MINE)
			*outstream << "\n  - Min. Energy : " << minEnergy;
		if (simOutMode != OUT_NO)
		        *outstream << "\n  - Out freq.   : " << outFreq;
		if (parser.argExist("countTarget"))
			*outstream << "\n  - Target str. : " << targetMoveStr;
		if (parser.argExist("final"))
			*outstream << "\n  - Final str.  : " << absMoveStrFinal;
		*outstream <<std::endl;
	}
	
	if (outHDF) {
		if (parser.argExist("title"))
			hdf5writer->write_attribute("Title", title);
		hdf5writer->write_attribute("Version", BIN_PACKAGE_VERSION);
		hdf5writer->write_attribute("Lattice", lattice->getDescriptor()->getName());
		hdf5writer->write_attribute("Energy file", parser.getStrVal("energyFile"));
		hdf5writer->write_attribute("Alphabet", alphString);
		hdf5writer->write_attribute("Elt. length", alphElementLength);
		hdf5writer->write_attribute("Sequence", seqStr);
		hdf5writer->write_attribute("Move set", moves);
		hdf5writer->write_attribute("Ribosome", (ribosome ? "tethered" : "untethered"));
		hdf5writer->write_attribute("Ribosome release", (ribosomeRelease ? "yes" : "no"));
		hdf5writer->write_attribute("Simulations", (unsigned int)runs);
		hdf5writer->write_attribute("Seed", (unsigned int)seed);
		hdf5writer->write_attribute("Modified seed", (unsigned int)(modified_seed + seed));
		hdf5writer->write_attribute("kT", (float)kT);
		hdf5writer->write_attribute("Max. steps", elongationSchedule[0]);

		{
			std::ostringstream schedule;
			if (parser.argExist("elongationSchedule")) {
				for (unsigned int i = 4; i < seqStr.size(); ++i)
					schedule << elongationSchedule[i] << " ";
			}
			hdf5writer->write_attribute("Custom elongation schedule",
						    parser.argExist("elongationSchedule") ?
						    schedule.str() : "no");
		}

		if (parser.argExist("fullLengthSteps"))
			hdf5writer->write_attribute("Full length steps",
						    (unsigned int)parser.getIntVal("fullLengthSteps"));

		hdf5writer->write_attribute("Proportional sim",
					    sequenceDependentSimLength ? "true" : "false" );
		if (minEnergy != DEFAULT_MINE)
			hdf5writer->write_attribute("Min. energy", (float)minEnergy);
		if (simOutMode != OUT_NO)
			hdf5writer->write_attribute("Out freq.", outFreq);
		if (parser.argExist("final"))
			hdf5writer->write_attribute("Final str", absMoveStrFinal);
		if (parser.argExist("countTarget"))
			hdf5writer->write_attribute("Count target str", targetMoveStr);
	}

	if (verbosity > 0) {
		*outstream	<< "\n Folding simulations :"
					<< "\n ====================="
					<<std::endl;
	}

	int foundMinE = 0;
	int foundFinalStructure = 0;
	WAC_LP_final *wac_final = NULL;
	if (parser.argExist("final")) {
		wac_final = new WAC_LP_final(absMoveStrFinal,*latticeDescriptor);
	}

	  // set up timer
	biu::Timer globalTime;
	globalTime.start();
	  // extra time counter for aborted runs (no elongatable chains)
	double timeAbortedRuns = 0.0;
	
	  // start simulation runs
	for (int doneRuns=0; doneRuns<runs; doneRuns++)
	{
		if (outHDF)
			hdf5writer->create_trajectory_group();

		///////////////////  BEGIN CHAIN GROWTH ///////////////////////////
			
		  // holds the energy of the current structure to extend 
		  // (update after folding simulation needed)
		double curEnergy = 0.0;
		  // will be true if a chain was not elongatable
		bool runWasAborted = false;
		  // flag to store if the minimal energy was reached within this run
		bool successfulRunMinE = false;
		bool successfulRunFinal = false;
		
		  // start local time measurement
		biu::Timer localTime;
		localTime.start();
		
		if (verbosity > 0)
			*outstream 	<< " performing co-translational folding simulation " << (doneRuns+1)
						<< "\n";
	
		  // hold the current chain length
	//	size_t curLength = 2;
		size_t curLength = 3;

		// track total number of steps in simulation
		size_t totalSteps = 0;		
		
		
		  // initialize with a ONE move
		absMoveStr = latticeDescriptor->getAlphabet()->getString(
						biu::MoveSequence( 1,
								latticeDescriptor->getNeighborhood().getElementByIndex(0).getMove())
						);
		  // create initial structure based on CT only
		for (size_t cL=2; cL <= curLength; cL++)
		{
			
			  // get subsequence of current intermediate chain
			biu::Sequence seq = alph->getSequence( seqStr.substr(0,(cL+1)*alphElementLength) );
			
			  // the set of all valid elongations including Boltzmann weight
			std::set< std::pair<std::string,double> > elongations;
			typedef std::set< std::pair<std::string,double> > ESET;
			  // partition function of all entries in elongations
			double curElongZ = 0.0;
			
			//////////////////// CALCULATE ALL ELONGATIONS //////////////////
	
			  // convert current data
			biu::MoveSequence actMoves = lattice->parseMoveString(absMoveStr);
			  // get 3D points of current relative move string from store
			biu::IPointVec actStruct = lattice->absMovesToPoints( actMoves );
			
			  // extend structure with dummy monomer
			actStruct.push_back(biu::IntPoint(0,0,0));
			actMoves.push_back(latticeDescriptor->getAlphabet()->getElement(0));
			  // try all possible extensions
			for (	biu::LatticeNeighborhood::const_iterator actNeigh = lattice->getNeighborhood().begin(); 
					actNeigh != lattice->getNeighborhood().end(); actNeigh++ )
			{
				bool isSelfavoiding = true;
				double energyGain = 0.0;
				
				if (ribosome && (actNeigh->getX() < 0)) {
					// Ribosome is a plane parallel to yz axes
					// Disallow extension in -x.
					continue;
				}
				
				  // get correct last position according to current neighbor vector
				actStruct[cL] = actStruct[cL-1] + *actNeigh;
				
				// do selfavoiding check and calculate energy gain
				for (size_t i = 0; i < cL-1; i++) {
					  // selfavoidingness check
					if ( actStruct[i] == actStruct[cL] ) {
						isSelfavoiding = false;
						break;
					}
					  // get distance dependent energy contribution
					energyGain += energy->getEnergy(	seq[i], 
														seq[cL],
														actStruct[i],
														actStruct[cL]);
				}
				  // store if necessary
				if (isSelfavoiding 
					&& (actStruct.size() == seqStr.size() || hasFreeNeighbor(actStruct,cL,*lattice))) 
				{
					  // set according move
					*(actMoves.rbegin()) = actNeigh->getMove();
					  // get move string
					const std::string curElongStr = lattice->getString(actMoves);
					  // get Boltzmann weight of energy change of current elongation structure
					const double BW = std::exp( - energyGain / kT );
					  // store
					elongations.insert( ESET::value_type( curElongStr, BW ));
					  // update partition function if current structure was not present
					curElongZ += BW;
				}
			}
			
			//////////////////// CHECK IF ONE ELONGATION POSSIBLE //////////////////
			
			if ( elongations.size() == 0 )
			{
				
				if (verbosity > 0)
				{
					*outstream 	<< "  = ABORTED due to non-elongatable chain !!! --> retrying ..."
								<< "\n";
				}
				  // get "lost" time of unsuccessful run
				timeAbortedRuns += localTime.stop();
				  // set flag to abort elongation iteration
				runWasAborted = true;
				  // reduce the number of done runs since this one was not successful
				doneRuns--;
				  // break current co-translational simulation / run
				break;  /////////////    BREAK ELONGATION  /////////////////
			}
			
			//////////////////// SELECT ONE ELONGATION //////////////////
			
			  // get random number in interval [0,Z]
			const double randNumber = curElongZ * (double(biu::RNF::getRN()) / double(biu::RNF::getMaxRN()));
			  // select neighbor according to Boltzmann weight sum
			double curBWsum = 0.0;
			for (ESET::const_iterator it=elongations.begin(); it!=elongations.end(); it++) {
				  // update Boltzmann weight sum
				curBWsum += it->second;
				  // check if selection found
				if (curBWsum >= randNumber) {
					  // update move string
					absMoveStr = it->first;
					  // update energy
					bool seqShared = true;
					bool isAbsMove = true;
					  // create lattice protein instance
					biu::LatticeProtein_Ipnt 
						latProt(lattice,energy,&seq,seqShared,absMoveStr,isAbsMove,ribosome);
					curEnergy = latProt.getEnergy();
	
					break;
				}
			}
		}
	
		
		assertbiu( absMoveStr.size() == (latticeDescriptor->getAlphabet()->getElementLength()*curLength), "absolute moves initialization was not successful");
		
		// Stopping criterion. Increase by 1 in the case of free chain simulation following translation.
		unsigned int stopLength = seqStr.size() / alphElementLength;
		if (parser.argExist("fullLengthSteps"))
			++stopLength;

		// curLength = 3 --> curLength = 4
		// "current length" is a misnomer. For each iteration, a residue is added.
		// Thus, the chain length simulated is actually curLength + 1.
		for (curLength++; (!runWasAborted) && curLength < stopLength; curLength++) {
			
			// TODO : walk accepts only elongatable structures (SAC specialization)
			// TODO : mem-leak ?!? (no global folding mode [maxSteps == 0])
			
			
			  // get subsequence of current intermediate chain
			biu::Sequence seq = alph->getSequence( seqStr.substr(0,(curLength+1)*alphElementLength) );
	
			///////////////////  ELONGATE CHAIN  ///////////////////////////

			if (curLength < seqStr.size() / alphElementLength)
			{
				  // the set of all valid elongations including Boltzmann weight
				std::set< std::pair<std::string,double> > elongations;
				typedef std::set< std::pair<std::string,double> > ESET;
				  // partition function of all entries in elongations
				double curElongZ = 0.0;
				
				//////////////////// CALCULATE ALL ELONGATIONS //////////////////
	
				  // convert current data
				biu::MoveSequence actMoves = lattice->parseMoveString(absMoveStr);
				  // get 3D points of current relative move string from store
				biu::IPointVec actStruct = lattice->absMovesToPoints( actMoves );
				
				  // extend structure with dummy monomer
				actStruct.push_back(biu::IntPoint(0,0,0));
				actMoves.push_back(latticeDescriptor->getAlphabet()->getElement(0));
				  // try all possible extensions
				for (	biu::LatticeNeighborhood::const_iterator actNeigh = lattice->getNeighborhood().begin(); 
						actNeigh != lattice->getNeighborhood().end(); actNeigh++ )
				{
					bool isSelfavoiding = true;
					double energyGain = 0.0;

					if (ribosome && (actNeigh->getX() < 0)) {
						// Ribosome is a plane parallel to yz axes
						// Disallow extension in -x.
						continue;
					}
				
					  // get correct last position according to current neighbor vector
					actStruct[curLength] = actStruct[curLength-1] + *actNeigh;
					
					// do selfavoiding check and calculate energy gain
					for (size_t i = 0; i < curLength-1; i++) {
						  // selfavoidingness check
						if ( actStruct[i] == actStruct[curLength] ) {
							isSelfavoiding = false;
							break;
						}
						  // get distance dependent energy contribution
						energyGain += energy->getEnergy(	seq[i], 
															seq[curLength],
															actStruct[i],
															actStruct[curLength]);
					}
					  // store if necessary
					if (isSelfavoiding 
						&& (actStruct.size() == seqStr.size() || hasFreeNeighbor(actStruct,curLength,*lattice))) 
					{
						  // set according move
						*(actMoves.rbegin()) = actNeigh->getMove();
						  // get move string
						const std::string curElongStr = lattice->getString(actMoves);
						  // get Boltzmann weight of energy change of current elongation structure
						const double BW = std::exp( - energyGain / kT );
						  // store
						elongations.insert( ESET::value_type( curElongStr, BW ));
						  // update partition function if current structure was not present
						curElongZ += BW;
					}
				}
				
				//////////////////// CHECK IF ONE ELONGATION POSSIBLE //////////////////
				
				if ( elongations.size() == 0 )
				{
					
					if (verbosity > 0)
					{
						*outstream 	<< "  = ABORTED due to non-elongatable chain !!! --> retrying ..."
									<< "\n";
					}
					  // get "lost" time of unsuccessful run
					timeAbortedRuns += localTime.stop();
					  // set flag to abort elongation iteration
					runWasAborted = true;
					
					  // break current co-translational simulation / run
					break;  /////////////    BREAK ELONGATION  /////////////////
				}
				
				//////////////////// SELECT ONE ELONGATION //////////////////
				
				  // get random number in interval [0,Z]
				const double randNumber = curElongZ * (double(biu::RNF::getRN()) / double(biu::RNF::getMaxRN()));
				  // select neighbor according to Boltzmann weight sum
				double curBWsum = 0.0;
				for (ESET::const_iterator it=elongations.begin(); it!=elongations.end(); it++) {
					  // update Boltzmann weight sum
					curBWsum += it->second;
					  // check if selection found
					if (curBWsum >= randNumber) {
						  // update move string
						absMoveStr = it->first;
						break;
					}
				}
			}
			assertbiu( absMoveStr.size() == (latticeDescriptor->getAlphabet()->getElementLength()*curLength), "resulting move string is of wrong length");
		
			/*
			 * Building Protein related objects.
			 */
			bool seqShared = true;
			bool isAbsMove = true;
			  // create lattice protein instance
			biu::LatticeProtein_Ipnt latProt(lattice,energy,&seq,seqShared,absMoveStr,isAbsMove,ribosome);

			if (ribosome && !latProt.isRibosomeValid())
			{
				std::cerr      << "Error: cannot anchor move sequence \'"
					       << absMoveStr << "\'"
					       << " to a wall. Not ribosome valid"
					       << std::endl;
				return DATA_ERROR;
			}

			if (curLength == seqStr.size() && ribosomeRelease && ribosome) {
				// Release from ribosome if chain is complete.
				latProt = biu::LatticeProtein_Ipnt(lattice,energy,&seq,seqShared,absMoveStr,isAbsMove,!ribosome);
			}

			  // update energy
			curEnergy = latProt.getEnergy();
	
			///////////////////  BEGIN FOLDING SIMULATION  ///////////////////////////
			
			/*
			 * Building State related objects
			 */
			// shared PullMoveDecoder
			biu::LatticeMoveSet* moveSet;
			S_LP* s;
			if (moves == OPTION_PULLM) {
				moveSet = new biu::PullMoveSet(&pmd, true);
				s = new S_LP_PullM(&latProt, moveSet);
			}
			else if (moves == OPTION_PIVOTM) {
				moveSet = new biu::PivotMoveSet(lattice);
				s = new S_LP_PivotM(&latProt, moveSet);
			}
			
			// set simulation time
			size_t simTime = elongationSchedule[curLength];
			// set simulation length depending on sequence length
			if (sequenceDependentSimLength) {
				simTime *= seq.size();
			}

			/*
			 * Building Walk related objects
			 */
			WAC_Signal wac_s(&stopFlag);    // stop upon change to stopFlag
			WAC_MinEnergy wac_e(minEnergy); // stop if minEnergy (default minEnergy is (double)INT_MIN)
			WAC_OR wac_or(wac_s, wac_e);	// contains wac_s and wac_e
			WAC_MaxLength wac_l(simTime);	// stop if length reached
			WAC_OR wac(wac_or, wac_l);	// contains wac_e, wac_l, and wac_s

				
			/* 
			 * Executing walk
			 */
		
			// build StateCollector
			switch(simOutMode)
				{
				case OUT_ES:
					sc = outHDF ? new SC_OutAbs(hdf5writer.get(), absMoveStr.length(), outFreq, totalSteps)
						: new SC_OutAbs(*outstream, absMoveStr.length(), outFreq, totalSteps);
					// sc = new SC_OutAbs(*simOut, absMoveStr.length(), outFreq, totalSteps);
					break;
				case OUT_E:
					sc = outHDF ? new SC_OutEnergy(hdf5writer.get(), outFreq, totalSteps)
						: new SC_OutEnergy(*outstream, outFreq, totalSteps);
					// sc = new SC_OutEnergy(*simOut, outFreq, totalSteps);
					break;
				case OUT_NO:
					sc = new SC_MinE();
					break;
				}
				
			if (((curLength+1) >= seqStr.size()) && parser.argExist("countTarget"))
			{
				sc->defineTarget(targetMoveStr, *latticeDescriptor);
			}

				
			// perform simulation...
			// if protein is full length, and last stage of simulation has a final
			//   target structure
			if ( ((curLength+1) >= seqStr.size()) && parser.argExist("final") ) {

				// curWac: wac_signal or wac_energy or normal length or final structure
				WAC_OR curWac(wac, *wac_final);
				WalkMC::walkMC(s, sc, &curWac, kT);
				successfulRunFinal = wac_final->abort(sc);
			} else {
				// "normal run"
				WalkMC::walkMC(s, sc, &wac, kT);
			}
				
			// if SIGINT or SIGTERM was caught
			if (stopFlag) {
				if (outHDF)
					hdf5writer.reset(); // unique_ptr reset deletes object
				sa.sa_handler = SIG_DFL;
				sigaction(SIGINT, &sa, NULL);
				raise(SIGINT);
			}

			// normal printing of final structure (if excluded by outFreq)
			if ( ((sc->size()-1) % outFreq) )
				sc->outputLast();

			// update 
			curEnergy = sc->getLastAdded()->getEnergy();
			absMoveStr = sc->getLastAdded()->toString();
			totalSteps += sc->size() - 1;
			if (absMoveStr.find("(") != std::string::npos) {
				absMoveStr = absMoveStr.substr(0,absMoveStr.find("("));
			}
			assertbiu( absMoveStr.size() == (latticeDescriptor->getAlphabet()->getElementLength()*curLength), "resulting move string is of wrong length");
				
			// check if lower energy bound was reached
			if (parser.argExist("minE") && curEnergy < minEnergy+DELTA_E_ADD) {
				successfulRunMinE = true;
			}
			minE = sc->getMinE();
			targetCountFraction = sc->getTargetCountFraction();
			targetConfEnergy = sc->getTargetEnergy();
			stepsToReachTarget = sc->getStepsToReachTarget();

			delete sc;
			delete s;  // who is deleting this? --> get segfault for explicit deletion
			delete moveSet;
				
			///////////////////  END FOLDING SIMULATION  ///////////////////
			
		} // for all elongations
	
		if (outHDF)
			hdf5writer->close_trajectory_group(
				successfulRunMinE, successfulRunFinal,
				targetCountFraction, targetConfEnergy,
				stepsToReachTarget);

		  // check if last run was successful or aborted
		if (runWasAborted) {
			  // reduce number of successful runs
			doneRuns--;
			  // delete the trajectory group
			if (outHDF)
				hdf5writer->delete_last_group();
			  // skip remainder of current run
			continue;
		}

		
		  // check if lower energy bound was reached
		if (!successfulRunMinE && parser.argExist("minE") && curEnergy < minEnergy+DELTA_E_ADD) {
			successfulRunMinE = true;
		}
		
		  // check if final structure was reached
		if (!successfulRunFinal && parser.argExist("final")) {
			  // get all symmetries of desired final structure
			typedef std::set<biu::MoveSequence> MSset; 
			MSset sym = latticeDescriptor->getAllSymmetricSequences(
							latticeDescriptor->getSequence(absMoveStrFinal));
			  // check if the final structure is among them
			successfulRunFinal = sym.find(latticeDescriptor->getSequence(absMoveStr))
									!= sym.end();
		}
		

		  // check if minimal energy was reached --> set counter
		if (successfulRunMinE) {
			foundMinE++;
		}
		  // check if final structure was reached --> set counter
		if (successfulRunFinal) {
			foundFinalStructure++;
		}
	
		  // get and print normalized final move string 
		absMoveStr = latticeDescriptor->getString(
						latticeDescriptor->normalizeSequence( 
							latticeDescriptor->getSequence(absMoveStr) ) );
		*outstream	<< "  - final structure / E : " <<absMoveStr <<" " << curEnergy
					<< "\n";

		if (parser.argExist("countTarget"))
			*outstream      << "  - target count fraction  : " << targetCountFraction <<"\n";
		
		  // print out local time used for current simulation
		if (verbosity > 0 && timing)
			*outstream	<< "  - computation time  : " 
						<< localTime.stop() <<" ms"
						<< "\n";
		if (verbosity > 1)
		{
			*outstream
						<< "  - simulation steps  : " << totalSteps <<"\n"
						<< "  - minimal E reached : " << minE
						<< std::endl;
		}

		///////////////////  END ELONGATION  ///////////////////
	
	} // for all runs
	
	///////////////////  END REPEATED RUNS  ///////////////////
	
	double globalTimeElapsed = globalTime.stop() - timeAbortedRuns;
	
	if (verbosity > 0) {
		*outstream	<< "\n Results :"
					<< "\n ========="
					<< "\n";
	}

	// print out overall time used
	if (timing)
		*outstream	<< "\n  overall computation time : " 
					<< globalTimeElapsed <<" ms"
					<< "\n";
	
	// print out hit percentage
	if (parser.argExist("minE")) {
		*outstream 	<< "\n  simulations reaching the minimal energy : "
					<< foundMinE << " / " << runs << " = "
					<< ((float) foundMinE / (float) runs)*100 << "%" 
					<< "\n";
	}
	if (parser.argExist("final")) {
		*outstream 	<< "\n  simulations reaching the final structure : "
					<< foundFinalStructure << " / " << runs << " = "
					<< ((float) foundFinalStructure / (float) runs)*100 << "%" << std::endl;
	}
	
	// note hdf5writer will be deleted at the end (smart pointer)

	  // final outstream clearing
	*outstream <<std::endl;

	  // clear memory
	if (energy != NULL) delete energy;
	if (energyMatrix != NULL) delete energyMatrix;
	if (alph != NULL) delete alph;
	if (lattice != NULL) delete lattice;
	if (latticeDescriptor != NULL) delete latticeDescriptor;
	if (wac_final != NULL) delete wac_final;
	
	return 0;
}
