
#include <iostream>

#include "biu/util/biu-string.h"
void test_biu_string() {

	std::string tstr = "lalalahierllhierala";
	std::string query = "hier";
	
	biu::util::findQueryString(tstr, query);
	
}


#include "biu/util/Util_String.h"
void test_util_string() {
	std::string test = "lululuuuuu";
	std::cout <<"string='" <<test <<"' toUCase='" <<(biu::util::Util_String::str2upperCase(test)) <<"'" <<std::endl;
}


int main(int argc, char **argv) {
	
	test_biu_string();
	test_util_string();
	
	return 0;
}
