

#include "biu/OptionParser.hh"
#include <iostream>

int main(int argc, char **argv) {

	biu::OptionMap allowedArgs;
	allowedArgs.push_back(biu::COption("text", false, biu::COption::STRING, "ein string parameter"));
	allowedArgs.push_back(biu::COption("max", false, biu::COption::INT, "ein integer parameter"));
	
	std::string infotext = "Das ist der Infotext zu diesem tollen programm!";

	// konstruktor macht im fehlerfall usageoutput
	biu::COptionParser opts = biu::COptionParser(allowedArgs, argc, argv, infotext);

	if (opts.noErrors()) {
		if (opts.argExist("text")) {
		    std::string text = opts.getStrVal("text");
		    std::cout <<" parameter 'text'='" <<text <<"'" <<std::endl;
		}
	}
	
	return 0;
}
