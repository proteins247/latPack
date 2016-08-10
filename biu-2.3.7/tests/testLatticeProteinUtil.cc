

#include <iostream>

#include "LatticeProteinUtil.hh"

void
testMoveStringOutput( const std::string moveString );

int main( int argc, char ** argv ) {
	
	testMoveStringOutput( "FLBBRRB" );
	testMoveStringOutput( "RFRBBRFRF" );
	
	return 0;
}

void
testMoveStringOutput( const std::string moveString )
{
	biu::LatticeDescriptorSQR ld;
	
	std::cout <<" derive move sequence of   : '" <<moveString <<"'" <<std::endl;
	biu::MoveSequence moves = ld.getSequence( moveString );
	
	std::cout <<" coordinate representation : \n" <<std::endl;
	
	std::cout <<toASCII( moves ) <<std::endl;
}
