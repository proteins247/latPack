// $Id: OptionParser.cc,v 1.2 2016/08/08 12:41:57 mmann Exp $

// OptionParser.cpp: Implementierung der Klasse COptionParser.
/*

	<author>	Martin Mann					</author>
	<created>	9.10.2005					</created>

	<info>									</info>
*/
//
//////////////////////////////////////////////////////////////////////

#include "biu/OptionParser.hh"
#include <iostream>

namespace biu {
	
//////////////////////////////////////////////////////////////////////
// statics aus COption
//////////////////////////////////////////////////////////////////////

std::string COption::DEF_INIT = "§$%&/()=";
std::vector<std::string> COption::TYPE_NAME = COption::initTypeNames();
int COptionParser::OUTPUT_LINE_LENGTH = 80;

//////////////////////////////////////////////////////////////////////
// Konstruktion
//////////////////////////////////////////////////////////////////////


// Konstruktor
COptionParser::COptionParser(OptionMap _options, int argc, char** argv, std::string infoText) 
: 	opt(_options), 
	programName(std::string(argv[0])), 
	infoText(infoText), 
	errorOccured(false), 
	maxOName(0) {
	if (programName.find_last_of("/\\")!=std::string::npos)
		programName = programName.substr(programName.find_last_of("/\\")+1);
	parseOpt(argc, argv);


	for (OptionMap::size_type i=0; i<opt.size(); i++) {
		mopt[opt[i].option] = i;				// merke mit zuordnung fuer "direkten" zugriff
		if (opt[i].option.size() > maxOName)	// merke laenge des laengesten optionsnamens
			maxOName = opt[i].option.size();
	}
} // construction()

// zeigt ob Probleme beim Parsen vorliegen
bool COptionParser::noErrors() {
	return !errorOccured;
}


void COptionParser::coutLineBreaking(std::string text, std::ostream &os,const int emptyHeadSize,const int lineLength) const {

	std::string line = "", emptyHead = std::string(emptyHeadSize, ' ');

	size_t nextCutPos = std::max(lineLength-emptyHeadSize,0);
	std::ostringstream oss;

	bool firstOut = true, cutted = false;

	while ((int)text.size() >  (lineLength-emptyHeadSize)) {
		cutted = false;
			// suche zeilenumbruch zum trennen
		        // (search line break to separate)
		nextCutPos = text.find_last_of("\n",lineLength-emptyHeadSize);
		if (nextCutPos == text.npos) { 	// nixs -> suche leerzeichen zum trennen (search blank to separate)
			nextCutPos = text.find_last_of(" ",lineLength-emptyHeadSize);
		}
		if (nextCutPos == text.npos) {	// notwendige worttrennung (necessary word separation)
			nextCutPos = lineLength-emptyHeadSize;
			cutted = true;
		}
		line += text.substr(0, nextCutPos);
		text = text.substr(nextCutPos+(cutted?0:1));
		// currently there is a bug
		os << line << std::endl;
		// std::cerr << line << std::endl;
		firstOut = false;
		line = emptyHead;
	}
	if (text.size() > 0)
		os <<(firstOut?"":emptyHead) <<text<<std::endl;
	os.flush();

}

// zeigt Parameter an
void COptionParser::coutUsage() const {
	std::cout <<"\n => usage:\t" <<programName ;
	OptionMap::size_type i;

	// ausgabe von nicht optionalen parametern
	for (i=0; i<opt.size(); i++) {
		if (opt[i].optional)
			continue;
		std::cout <<" -"<<opt[i].option <<(opt[i].retType==COption::BOOL?"":"=...");
	}
	std::cout <<" [...optional parameters]\n\n => complete parameterlist:\n";
	std::cout <<std::endl;

	int maxName = maxOName + 7, j=0, preStrLength = maxName +2;
	std::string pInfo = "", emptyPreStr = std::string(preStrLength,' ');

	// komplette parameterliste
	for (i=0; i<opt.size(); i++) {
		pInfo = std::string(opt[i].optional?"  [-":"   -") + opt[i].option + std::string(opt[i].retType==COption::BOOL?"":"=..");
		for (j = pInfo.size(); j<maxName; j++)
			pInfo += " ";
		std::cout <<pInfo <<(opt[i].optional?"] ":"  ");
		coutLineBreaking(opt[i].description, std::cout, preStrLength, OUTPUT_LINE_LENGTH);
		std::cout <<emptyPreStr<<"type=" <<COption::TYPE_NAME[opt[i].retType];
		if (opt[i].defValue.compare(COption::DEF_INIT) != 0) { // dann defaultwert vorhanden
			std::cout <<"  default=\"" <<opt[i].defValue <<"\"";
		}
		std::cout <<std::endl;
	}

	// zusaetzliche informationen
	std::cout <<"\n\n => informations:\n\n";
	coutLineBreaking(infoText, std::cout, 0, OUTPUT_LINE_LENGTH);
	std::cout <<"\n"<<std::endl;

	// nochmal usage
	std::cout <<" => usage:\t" <<programName ;

	// ausgabe von nicht optionalen parametern
	for (i=0; i<opt.size(); i++) {
		if (opt[i].optional)
			continue;
		std::cout <<" -"<<opt[i].option <<(opt[i].retType==COption::BOOL?"":"=...");
	}
	std::cout <<" [...optional parameters]"<<std::endl;
} // coutUsage()

// zeigt ob der entsprechende Parameter angegeben wurde
bool COptionParser::argExist(std::string option) {
	if (mopt.find(option) != mopt.end()) {// dann existiert die option
		return (opt[mopt[option]].exist || !(opt[mopt[option]].defValue.compare(COption::DEF_INIT) == 0)); // dann auch angegeben
	}
	return false;
}

// Parse-Funktionen

std::string COptionParser::getStrVal(std::string arg) {
	if (argExist(arg)) {
		return opt[mopt[arg]].strValue;
	} else {
		return "";
	}
}
char COptionParser::getCharVal(std::string arg) {
	if (argExist(arg)) {
		return opt[mopt[arg]].strValue[0];
	} else {
		return ' ';
	}
}
int COptionParser::getIntVal(std::string arg) {
	if (argExist(arg)) {
		int ret;
		std::istringstream is;
		is.str(opt[mopt[arg]].strValue);
		is >> ret;
		return ret;
	} else {
		return 0;
	}
}
float COptionParser::getFloatVal(std::string arg) {
	if (argExist(arg)) {
		float ret;
		std::istringstream is;
		is.str(opt[mopt[arg]].strValue);
		is >> ret;
		return ret;
	} else {
		return 0.0;
	}
}
double COptionParser::getDoubleVal(std::string arg) {
	if (argExist(arg)) {
		double ret;
		std::istringstream is;
		is.str(opt[mopt[arg]].strValue);
		is >> ret;
		return ret;
	} else {
		return 0.0;
	}
}
bool COptionParser::getBoolVal(std::string arg) {
		return argExist(arg);
}






// liefert formatierte Fehlerausgabe
void COptionParser::coutError(int error, std::string optionName, std::string errormsg) {
	errorOccured = true;
	std::cout <<"\n\tERROR : ";
	switch (error) {
		case ERR_NO_OPT : std::cout <<"unknown option : "; break;
		case ERR_WR_USE : std::cout <<"wrong usage : "; break;
		case ERR_WR_VAL : std::cout <<"wrong argument type : "; break;
		case ERR_NO_ARG : std::cout <<"needed argument not given : "; break;
		default : std::cout <<" errorcode = " <<error <<" : "; break;
	}
	std::cout <<"'" <<optionName <<"' " <<errormsg <<std::endl<<std::endl;
} // coutError()

// prueft, ob Parameter parsebar
bool COptionParser::isCastable(std::string val, int type) const {
	std::istringstream ss;
	ss.str(val);
	std::string dump;

	switch (type) {
		case COption::STRING :
			return true;
		case COption::CHAR	:
			char ac;
			if (!(ss >> ac)) {
				return false;
			} else {
				if ((ss >> dump) && dump.size() > 0 && dump.find_first_not_of(" \t")!=std::string::npos)
					return false;
				else
					return true;
			}
		case COption::INT :
			int ai;
			if (!(ss >> ai)) {
				return false;
			} else {
				if ((ss >> dump) && dump.size() > 0 && dump.find_first_not_of(" \t")!=std::string::npos)
					return false;
				else
					return true;
			}
		case COption::DOUBLE :
			double ad;
			if (!(ss >> ad)) {
				return false;
			} else {
				if ((ss >> dump) && dump.size() > 0 && dump.find_first_not_of(" \t")!=std::string::npos)
					return false;
				else
					return true;
			}
		case COption::FLOAT :
			float af;
			if (!(ss >> af)) {
				return false;
			} else {
				if ((ss >> dump) && dump.size() > 0 && dump.find_first_not_of(" \t")!=std::string::npos)
					return false;
				else
					return true;
			}
		default :
			return false;
	}
} // isCastable()

// parst Parameter
void COptionParser::parseOpt(int argc, char** argv) {
	std::vector<int> optStart; // Stores indices where argv holds an option flag
	OptionMap::size_type i,j,last;
	std::string actOpt, actOptName, actVal;
	for (int k=1; k < argc; k++) { // merke mir alle potentiellen optionsstarts
		if (argv[k][0] == '-')
		optStart.push_back(k);
	}
	std::string help1 = "-help", help2 = "--help";
	for (i=0; i<optStart.size(); i++) { // check for 'help' parameter
		if (	 help1.compare(argv[optStart[i]]) == 0 
				|| help2.compare(argv[optStart[i]]) == 0) 
		{
			coutUsage();
			errorOccured = true; // for program abort if wanted
			return;
		}
	}
	for (i=0; i<optStart.size(); i++) { // gehe optionen durch
		actOpt = std::string(argv[optStart[i]]); // current option flag
		if (i+1<optStart.size())  // entweder argumente bis zur naechsten option
			                  // (either argue arguments to the next option
			last = optStart[i+1];
		else  // oder alle nachfolgenden argumente bis zum ende zusammenfuegen
		      // or all subsequent arguments to the end)
			last = argc;

		for (j=optStart[i]+1; j<last; j++)  // erzeuge komplette option (create complete option)
			                            // appends all argv up to the next flag to actOpt
			actOpt.append(" "+std::string(argv[j]));

		actOptName = actOpt.substr(1, actOpt.find("=",0)-1); // optionsname zum pruefen (option name to check)
		for (j=0; j<opt.size(); j++) { // suche eintrag in opt (search entry in opt)
			                       // (finds int j such that opt[j] is our option)
			if (opt[j].option == actOptName)
				break;
		}
		if (j != opt.size()) {	// dann option tatsaechlich verfuegbar (Then option is actually available)
			if (opt[j].retType == COption::BOOL) { // extrabehandlund der booleans (handling booleans)
				if (actOpt.size() > actOptName.size()+1) { // dann falsche benutzung (wrong use)
					coutError(ERR_WR_USE, actOptName, "is a boolean and no input argument = '"+actOpt+"'");
				} else { // merken
					opt[j].strValue = "true";
					opt[j].exist = true;
				}
			} else {
				if (actOpt.find("=",0) == std::string::npos) { // dann keine zuweisung obwohl kein bool
					                                       // (then no assignment although no bool)
					coutError(ERR_WR_USE, actOptName, "is an input argument. use '='");
				} else {
					actVal = actOpt.substr(actOptName.size()+2);
					if (actVal.size() == 0) {
						coutError(ERR_WR_VAL, actOptName, "is given with empty input");
					} else if (isCastable(actVal, opt[j].retType)) {  // wenn typecast moeglich
						                                          // (typecast possible)
						opt[j].strValue = actVal; // dann speichere (then save)
						opt[j].exist = true;
					} else {  // fehler da typecast nicht moeglich
						  // (error, typecast not possible)
						coutError(ERR_WR_VAL, actOptName, "in argument '"+actOpt+"'");
					}
				}
			}
		} else { // fehlermeldung (error message)
			coutError(ERR_NO_OPT, actOptName, "in argument '" + actOpt + "'");
			break;
		}
	}
	for (i=0; i<opt.size(); i++) {
		if(!(opt[i].optional) && !(opt[i].exist)) {
			coutError(ERR_NO_ARG,opt[i].option,"check usage");
		}
	}
} // parseoption()

} //namespace biu
