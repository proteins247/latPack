
#include "biu/Matrix.hh"
#include <iostream>


int main(int argc, char ** argv) {
	
	typedef biu::Matrix<int> mathInt;
	
	mathInt m1 = mathInt(4,4), m2 = mathInt(4,2);
	
	std::cout << "matrix m1 = \n" <<m1 <<"matrix m2 = \n" <<m2 <<std::endl;
	std::cout << "m1*m2 = \n" <<m1*m2 <<std::endl;
	std::cout << "m1 is null ? " <<((m1==0)?"true":"false") <<std::endl;
	
	return 0;
}
