
#include <iostream>
#include <fstream>
#include <time.h>
#include <csignal>
#include <cstring>

#include <biu/OptionParser.hh>
#include <biu/RandomNumberGenerator.hh>
#include <biu/RandomNumberFactory.hh>
#include <biu/LatticeDescriptorCUB.hh>
#include <biu/LatticeDescriptorFCC.hh>
#include <biu/LatticeDescriptorSQR.hh>
#include <biu/LatticeModel.hh>
#include <biu/LatticeProtein_Ipnt.hh>
#include <biu/PullMoveSet.hh>
#include <biu/PivotMoveSet.hh>
#include <biu/Timer.hh>

#include <ell/Walk.hh>
#include <ell/WalkAbortionCriterion.hh>
#include <ell/protein/S_LP_PullM.hh>
#include <ell/protein/S_LP_PivotM.hh>
#include <ell/protein/WAC_LP.hh>

#include "HDF5_support.hh"

#include "SC_OutEnergy.hh"
#include "SC_OutAbs.hh"
#include "energyFileSupport.hh"
#include "version.hh"

using namespace ell;

// error values
static const int PARSE_ERROR = 1;
static const int DATA_ERROR = 2;
static const int IO_ERROR = 3;

// default parameters
static const std::string DEFAULT_MAXSTEPS = "50";
static const std::string DEFAULT_KT = "0.3";
static const std::string DEFAULT_SEED = "1";
static const double DEFAULT_MINE = (double)INT_MIN;
static const std::string DEFAULT_LATTICE = "CUB";
static const std::string DEFAULT_MOVES = "PullM";
static const std::string DEFAULT_RUNS = "1";
static const std::string DEFAULT_OUTFREQ = "1";
static const std::string DEFAULT_ELEMENTLENGTH = "1";

// default data

static const bool optional = true;
static const std::string OPTION_CUB = "CUB";
static const std::string OPTION_SQR = "SQR";
static const std::string OPTION_FCC = "FCC";
static const std::string OPTION_PULLM = "PullM";
static const std::string OPTION_PIVOTM = "PivotM";

// constants
static const double DELTA_E_ADD = 0.0001;

 // possible output modes
// enum OUT_MODE { OUT_NO, OUT_E, OUT_ES };
// definition moved to SC_MinE.hh - VZ

// infotexts
static const std::string infotext = 
	"LatFold: lattice protein folding simulation on different"
	" lattices using the specified move set and a metropolis walk.\n"
	"\n"
	"The simulation stops if a given minimal energy is reached or"
	" a given number of simulation steps have been performed.\n"
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
	"protein sequence (sequence valid for alphabet from energy file)";
static const std::string absInfo =
	"absolute move sequence "
	"(the moves are encoded using: "
	"F/B:+-x, "
	"L/R:+-y, "
	"U/D:+-z) "
	"[defaults to a series of equal moves e.g. 'F' for SQR-lattice]";
static const std::string titleInfo =
	"Give trajectory a title; most relevant for HDF5 output";
static const std::string ktInfo =
	"kT parameter for metropolis criterion "
	"(double value from [0, infinity))";
static const std::string maxstepsInfo =
	"walk ends after accepting [maxSteps] states";
static const std::string minenergyInfo =
	"walk ends if energy gets below or equal to [minE]";
static const std::string seedInfo =
	"seed for random number generator "
	"[uses biu::RNG_ARS4x32, a counter-based generator from the Random123 library]";
static const std::string runsInfo =
	"number of folding simulations to perform";
static const std::string latticeInfo =
	"which lattice to use: CUB, SQR or FCC";
static const std::string ofileInfo =
	"write output of simulations to filename (HDF5). if equal to 'STDOUT' it is written to standard output (text).";
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
	"simulate with a non-interacting wall that acts as a barrier and anchors the last residue. 3D-only";

// csignal for SIGINT, SIGTERM handling
volatile sig_atomic_t stopFlag = 0;

static void handler(int signum)
{
	stopFlag = signum;
}


int main(int argc, char** argv) {
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
			"seq", !optional, 
			biu::COption::STRING, 
			seqInfo));
	
	options.push_back(biu::COption(
			"abs", 
			optional, 
			biu::COption::STRING, 
			absInfo));

	options.push_back(biu::COption(
			"title", 
			optional, 
			biu::COption::STRING, 
			titleInfo));
	
	options.push_back(biu::COption(	
			"energyFile", 
			!optional, 
			biu::COption::STRING, 
			"the contact energy function file that contains also the alphabet"));

	options.push_back(biu::COption(	
			"elementLength", 
			optional, 
			biu::COption::INT,
			"the character length of a single alphabet element in the energy file (not a tested feature)",
			DEFAULT_ELEMENTLENGTH));

	options.push_back(biu::COption(	
			"energyForDist", 
			true, 
			biu::COption::BOOL, 
			"if present a distance based energy function is used, otherwise a contact energy function is applied"));
	
	options.push_back(biu::COption(	
			"energyCalphaDist", 
			true, 
			biu::COption::DOUBLE, 
			"if '-energyForDist' is present, this value is used to scale the C_alpha distances of the given energy function to the C_alpha monomer distances in the used lattice model",
			"3.8"));

	options.push_back(biu::COption(
			"kT", 
			optional, 
			biu::COption::DOUBLE, 
			ktInfo, 
			DEFAULT_KT));
	
	options.push_back(biu::COption(
			"maxSteps", 
			optional, 
			biu::COption::DOUBLE,
			maxstepsInfo, 
			DEFAULT_MAXSTEPS));
	
	{
	std::ostringstream stringStream;
	stringStream << DEFAULT_MINE;
	options.push_back(biu::COption(
			"minE", 
			optional,
			biu::COption::DOUBLE,
			minenergyInfo ));
	}
	
	options.push_back(biu::COption(
			"final", optional, biu::COption::STRING,
			"if present: each folding simulation run is aborted if the given"
			" (or a symmetric) structure is visited."
			));

	options.push_back(biu::COption(
			"countTarget", optional, biu::COption::STRING,
			"if present: count number of times folding simulation visits the given"
			" (or symmetric) structure."
			));
		
	options.push_back(biu::COption(
			"targetDegradationScale", optional, biu::COption::DOUBLE,
			"if present: sets scale for number associated with cumulative survival of "
			"structure. (probability of survival falls with the number of steps "
			"that the structure is not native.)"
			));
		
	options.push_back(biu::COption(
			"seed", 
			optional, 
			biu::COption::INT,
			seedInfo, 
			DEFAULT_SEED));
	
	options.push_back(biu::COption(
			"runs",
			optional,
			biu::COption::INT,
			runsInfo,
			DEFAULT_RUNS));
	
	options.push_back(biu::COption(
			"lat",
			optional,
			biu::COption::STRING,
			latticeInfo,
			DEFAULT_LATTICE));
	
	options.push_back(biu::COption(
			"moveSet",
			optional,
			biu::COption::STRING,
			moveSetInfo,
			DEFAULT_MOVES));

	options.push_back(biu::COption(
			"ribosome",
			optional,
			biu::COption::BOOL,
			ribosomeInfo));
	
	options.push_back(biu::COption(
			"out",
			optional,
			biu::COption::CHAR,
			"output mode along the folding simulation: (N)o, (E)nergy, (S)tructure+Energy",
			"N"));
	
	options.push_back(biu::COption(
			"outFreq",
			optional,
			biu::COption::INT,
			"output frequency. print information every INT steps. in effect only if -out=E or -out=S",
			DEFAULT_OUTFREQ));
	
	options.push_back(biu::COption(
			"outFile",
			optional,
			biu::COption::STRING,
			ofileInfo,
			"STDOUT"));
	
	options.push_back(biu::COption(
			"outTiming",
			optional,
			biu::COption::BOOL,
			timingInfo));
	
	options.push_back(biu::COption(
			"s",
			optional,
			biu::COption::BOOL,
			"silent mode: only final mfe hit statistics printed"));
	
	options.push_back(biu::COption(
			"v",
			optional,
			biu::COption::BOOL,
			verbosityInfo));
	
	options.push_back(biu::COption(
			"vv",
			optional,
			biu::COption::BOOL,
			vvInfo));
	
	options.push_back(biu::COption(
			"help",
			optional,
			biu::COption::BOOL,
			helpInfo));
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
	double degradationScale;
	size_t maxLength;
	double minEnergy;
	biu::LatticeDescriptor* latticeDescriptor = NULL;
	biu::LatticeModel * lattice = NULL;
	size_t seed;
	int runs;
	unsigned int outFreq;
	unsigned int alphElementLength;
	std::ostream* outstream = &std::cout;
	ell::SC_MinE* sc;
	bool ribosome;
	bool timing;
	int verbosity;
	OUT_MODE simOutMode = OUT_NO;
	bool outHDF = false;
	std::unique_ptr<HDF5TrajWriter> hdf5writer;
	// std::ostream* simOut = &std::cout; will remove this
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
		    if ( seqStr.size()/alphElementLength <3) {
		    	std::cerr	<< "Error: seq must have at least length 3."
		    				<< std::endl;
		    	return PARSE_ERROR;
		    }
		}
		
		
		if (parser.argExist("abs")) 
		{
			absMoveStr = parser.getStrVal("abs");
			// check alphabet
			if (!latticeDescriptor->getAlphabet()->isAlphabetString(absMoveStr)) {
				std::string moveAlphStr = latticeDescriptor->getAlphabet()->getString(
											latticeDescriptor->getAlphabet()->getElement(0));
				for (size_t i=1; i<latticeDescriptor->getAlphabet()->getAlphabetSize();i++) {
					moveAlphStr += ",";
					moveAlphStr += latticeDescriptor->getAlphabet()->getString(
									latticeDescriptor->getAlphabet()->getElement(i));
				}
				std::cerr	<<"Error: given absolute move string '"
							<<absMoveStr
							<<"' is not valid for the absolute move alphabet {"
							<<moveAlphStr
							<<"} !\n";
				return PARSE_ERROR;
			}
			// check size
			if (latticeDescriptor->getAlphabet()->getSequence(absMoveStr).size()+1 != seqStr.size()/alphElementLength)
			{
				std::cerr	<< "Error: sequence and structure differ in size."
							<< std::endl;
				return PARSE_ERROR;
			}
		}
		else {
			absMoveStr = "";
			  // add a sequence of the first move of the move alphabet
			for (size_t i=1; i<seqStr.size()/alphElementLength; i++)
			{
				absMoveStr.append(
							latticeDescriptor->getAlphabet()->getString(
									latticeDescriptor->getAlphabet()->getElement(0))
								);
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
		
		if (parser.argExist("maxSteps")) 
		{
			long long int tmp = (long long int)parser.getDoubleVal("maxSteps");
			if (tmp < 0) {
				std::cerr << "Error: maxSteps must be >= 0." << std::endl;
				return PARSE_ERROR;
			}
			maxLength = tmp; 
		}
		
		if (parser.argExist("minE"))
		{
			minEnergy = parser.getDoubleVal("minE");
		} else {
			minEnergy = (double)INT_MIN;
		}
		
		if (parser.argExist("final"))
		{
			absMoveStrFinal = parser.getStrVal("final");
		}
		
		if (parser.argExist("countTarget"))
		{
			targetMoveStr = parser.getStrVal("countTarget");
		}
		
		if (parser.argExist("targetDegradationScale"))
		{
		        degradationScale = parser.getDoubleVal("targetDegradationScale");
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
												seqStr.size()/alphElementLength-1));
				// seqStr.size()-1 is the size of the move string for the structure
			}
			catch (File_opening_error) {
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
		if (parser.argExist("s")) verbosity = -1;
		
		timing = parser.argExist("outTiming");
		
	}
	else
	{
		return PARSE_ERROR;
	}
	
	/*
	 * Preparing Randomizer
	 */
	biu::RandomNumberGenerator* rng = new biu::RNG_ARS4x32();
	biu::RNF::setRNG( rng );
	delete rng;
	// To avoid using the same RN sequence for all simulations with the same user seed
	//   (but different sequences and/or temperature), a modified seed is generated
	//   based on the user seed, as well as temperature and sequence.
	// At the moment, different starting conformations will have the same seed however.
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
	
	/*
	 * Building Protein related objects.
	 */
	bool seqShared = true;
	bool isAbsMove = true;
	biu::Sequence seq = alph->getSequence(seqStr);
	biu::LatticeProtein_I* latProt = new biu::LatticeProtein_Ipnt
										(lattice,energy,&seq,seqShared,absMoveStr,
										 isAbsMove,ribosome);

	if (!latProt->isSelfAvoiding())
	{
		std::cerr	<< "Error: move sequence \'"
					<< absMoveStr << "\'"
					<< " is not selfavoiding."
					<< std::endl;
		return DATA_ERROR;
	}

	if (ribosome && !latProt->isRibosomeValid())
	{
		std::cerr      << "Error: cannot anchor move sequence \'"
					<< absMoveStr << "\'"
					<< " to a wall. Not ribosome valid"
					<< std::endl;
		return DATA_ERROR;
	}

	/*
	 * Building State related objects
	 */
	// shared PullMoveDecoder
	biu::LatticeMoveSet* moveSet;
	biu::PullMoveSet::PullMoveDecoder* pmd;
	S_LP* s;
	if (moves == OPTION_PULLM) {
		pmd = new biu::PullMoveSet::PullMoveDecoder(lattice);
		moveSet = new biu::PullMoveSet(pmd, true);
		s = new S_LP_PullM(latProt, moveSet);
	}
	else if (moves == OPTION_PIVOTM) {
		moveSet = new biu::PivotMoveSet(lattice);
		s = new S_LP_PivotM(latProt, moveSet);
	}

	
	/*
	 * Building Walk related objects
	 */
	WAC_MinEnergy wac_e(minEnergy);
	WAC_MaxLength wac_l(maxLength);
	WAC_OR wac_or(wac_e, wac_l);
	WAC_Signal wac_s(&stopFlag);
	WAC_OR wac(wac_or, wac_s);
	
	/* 
	 * Executing walk
	 */
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
					<< "\n  - Abs. moves  : " << absMoveStr
					<< "\n  - Move set    : " << moves
         				<< "\n  - Ribosome?   : " << (ribosome ? "tethered" : "untethered")
					<< "\n  - Simulations : " << runs
					<< "\n  - Seed (rand) : " << seed << " (" << modified_seed + seed << ")"
					<< "\n  - kT (MC)     : " << kT
					<< "\n  - Max. steps  : " << maxLength
					;
		if (minEnergy != DEFAULT_MINE)
			*outstream << "\n  - Min. Energy : " << minEnergy;
		if (simOutMode != OUT_NO)
		        *outstream << "\n  - Out freq.   : " << outFreq;
		if (parser.argExist("countTarget"))
		{
			*outstream << "\n  - Target str. : " << targetMoveStr;
			if (parser.argExist("targetDegradationScale"))
			        *outstream << "\n  - Target deg. : " << degradationScale;
		}
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
		hdf5writer->write_attribute("Abs. moves", absMoveStr);
		hdf5writer->write_attribute("Move set", moves);
		hdf5writer->write_attribute("Ribosome", (ribosome ? "tethered" : "untethered"));
		hdf5writer->write_attribute("Simulations", (unsigned int)runs);
		hdf5writer->write_attribute("Seed", (unsigned int)seed);
		hdf5writer->write_attribute("Modified seed", (unsigned int)(modified_seed + seed));
		hdf5writer->write_attribute("kT", (float)kT);
		hdf5writer->write_attribute("Max. steps", (unsigned int)maxLength);
		if (minEnergy != DEFAULT_MINE)
			hdf5writer->write_attribute("Min. energy", (float)minEnergy);
		if (simOutMode != OUT_NO)
			hdf5writer->write_attribute("Out freq.", outFreq);
		if (parser.argExist("final"))
			hdf5writer->write_attribute("Final str", absMoveStrFinal);
		if (parser.argExist("countTarget"))
		{
			hdf5writer->write_attribute("Count target str", targetMoveStr);
			if (parser.argExist("targetDegradationScale"))
			        hdf5writer->write_attribute("Degradation scale", (float)degradationScale);
		}
	}
	
	if (verbosity > 0) {
		*outstream	<< "\n Folding simulations :"
					<< "\n ====================="
					<<std::endl;
	}

	int foundMinE = 0;
	int foundFinalStructure = 0;
	  // set up timer
	biu::Timer globalTime;
	globalTime.start();
	for (int i=1; i<=runs; i++)
	{
		// // restart RNG
		// biu::RNF::getRNG().setSeed(seed + i);

		if (outHDF)
			hdf5writer->create_trajectory_group();

		// build StateCollector
		switch(simOutMode)
		{
		case OUT_ES:
			sc = outHDF ? new SC_OutAbs(hdf5writer.get(), absMoveStr.length(), outFreq, 0)
				: new SC_OutAbs(*outstream, absMoveStr.length(), outFreq, 0);
			break;
		case OUT_E:
			sc = outHDF ? new SC_OutEnergy(hdf5writer.get(), outFreq, 0)
				: new SC_OutEnergy(*outstream, outFreq, 0);
			break;
		case OUT_NO:
			sc = new SC_MinE();
			break;
		default:
			sc = NULL;
			*outstream 	<< "\n RUNTIME ERROR : state collector not set.. aborting here!" <<std::endl;
			exit(-1);
		}

		if (parser.argExist("countTarget"))
		{
			sc->defineTarget(targetMoveStr, *latticeDescriptor);
			if (parser.argExist("targetDegradationScale"))
			{
			        sc->setDegradationRate(degradationScale);
			}
		}
		
		if (verbosity > 0)
		{
			*outstream 	<< "\n  performing folding simulation " << i
						<< std::endl;
		}
		  // start local time measurement
		biu::Timer localTime;
		localTime.start();
		
		bool successfulRunMinE = false;
		bool successfulRunFinal = false;
		
		  // perform simulation
		if (parser.argExist("final")) {
			  // run that checks for final structure
			WAC_LP_final wac_final(absMoveStrFinal,*latticeDescriptor);
			WAC_OR curWac(wac, wac_final);
			WalkMC::walkMC(s, sc, &curWac, kT);
			successfulRunFinal = wac_final.abort(sc);
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
			
		  // get final structure
		double finalEnergy = sc->getLastAdded()->getEnergy();
		std::string finalAbsMoveStr = sc->getLastAdded()->toString();
		if (finalAbsMoveStr.find("(") != std::string::npos) {
			finalAbsMoveStr = finalAbsMoveStr.substr(0,finalAbsMoveStr.find("("));
		}
		  // output of final structure
		if (verbosity >= 0) {
			  // get and print normalized final move string 
			finalAbsMoveStr = latticeDescriptor->getString(
							latticeDescriptor->normalizeSequence( 
								latticeDescriptor->getSequence(finalAbsMoveStr) ) );
			*outstream	<< "  - final structure / E : " 
						<<finalAbsMoveStr <<" " <<finalEnergy
						<< "\n";
			if (parser.argExist("countTarget"))
			{
				*outstream      << "  - target count fraction  : " << sc->getTargetCountFraction() <<"\n";
				if (parser.argExist("targetDegradationScale"))
				    *outstream  << "  - survival sum fraction  : " << sc->getSurvivalSumFraction() <<"\n";
			}
		}

		  // print out local time used for current simulation
		if (verbosity > 0 && timing)
			*outstream	<< "  - computation time  : " 
						<< localTime.stop() <<" ms"
						<< std::endl;
		
		if (verbosity > 1)
		{
			*outstream
						<< "  - simulation steps  : " << (sc->size()-1) <<"\n"
						<< "  - minimal E reached : " << sc->getMinE()
						<< std::endl;
		}
		  // check if lower energy bound was reached
		if (parser.argExist("minE") && sc->getLastAdded()->getEnergy() < minEnergy+DELTA_E_ADD) {
			successfulRunMinE = true;
		}
		
		if (successfulRunMinE) {
			foundMinE++;
		}
		
		if (successfulRunFinal) {
			foundFinalStructure++;
		}
		
		if (outHDF)
			hdf5writer->close_trajectory_group(
				successfulRunMinE, successfulRunFinal,
				sc->getTargetCountFraction(),
				sc->getTargetEnergy(),
				sc->getStepsToReachTarget(),
				sc->getSurvivalSumFraction(),
				sc->getExtrapolatedOutput());
		
		delete sc;
	}
	
	double globalTimeElapsed = globalTime.stop();
	
	if (verbosity > 0) {
		*outstream	<< "\n Results :"
					<< "\n ========="
					<<std::endl;
	}

	// print out overall time used
	if (timing)
		*outstream	<< "\n  overall computation time : " 
					<< globalTimeElapsed <<" ms"
					<< std::endl;
	
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
	
	return 0;
}
