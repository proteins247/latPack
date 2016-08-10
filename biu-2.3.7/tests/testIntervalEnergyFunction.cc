#include <iostream>
#include <limits.h>
#include <biu/DistanceEnergyFunction.hh>


int main(int argc, char **argv) {

	using namespace biu;

	std::cout	<<"========== BEGIN OF TEST =================\n"
				<<" test class : IntervalEnergyFunction\n"
				<<"==========================================" <<std::endl;
				
		// erzeugen	
	std::cout <<"=> creating biu::Alphabet(ABCD,1) = alph" <<std::endl;
	Alphabet alph("ABCD",1);
	std::cout <<"=> creating biu::IntervalEnergyFunction( alph)" <<std::endl;
	IntervalEnergyFunction f(&alph);
	
		// getAlphabet()
	std::cout	<<"->getAlphabet()->getAlphabetSize() = " 
				<<f.getAlphabet()->getAlphabetSize() 
				<<std::endl;

	std::cout <<"=> getIntervalNum() = " <<f.getIntervalNum() <<std::endl;
	std::cout <<"=> getInterval( 4.0 ) == UINT_MAX ?  " 
				<<(f.getInterval( 4.0 )==UINT_MAX?true:false) <<std::endl;

	
	std::cout <<"=> creating biu::EnergyMatrix(4,4,1.0) = E1" <<std::endl;
	EnergyMatrix E1 = EnergyMatrix(	alph.getAlphabetSize(),
									alph.getAlphabetSize(),
									1.0 );
	std::cout <<E1 <<std::endl;
	
	std::cout <<"=> addInterval( E1, 2.0 ) = " <<f.addInterval(E1,2.0) <<std::endl;
	
	std::cout <<"=> creating biu::EnergyMatrix(4,4,2.0) = E1" <<std::endl;
	EnergyMatrix E2 = EnergyMatrix(	alph.getAlphabetSize(),
									alph.getAlphabetSize(),
									2.0 );
	std::cout <<E2 <<std::endl;

	std::cout <<"=> addInterval( E2, 3.0 ) = " <<f.addInterval(E2,3.0) <<std::endl;

	std::cout <<"=> creating biu::EnergyMatrix(4,4,3.0) = E1" <<std::endl;
	EnergyMatrix E3 = EnergyMatrix(	alph.getAlphabetSize(),
									alph.getAlphabetSize(),
									3.0 );
	std::cout <<E3 <<std::endl;
	
	std::cout <<"=> addInterval( E3, 14.0 ) = " <<f.addInterval(E3,14.0) <<std::endl;

	std::cout <<"=> getIntervalNum() = " <<f.getIntervalNum() <<std::endl;
	for (size_t i=0; i<f.getIntervalNum(); i++) {
		std::cout <<"=> getIntervalMax( "<<i <<") = " << f.getIntervalMax( i ) <<std::endl;
	}
	std::cout <<"=> getInterval( 0.0 ) = " << f.getInterval( 0.0 ) <<std::endl;
	std::cout <<"=> getInterval( 2.0 ) = " << f.getInterval( 2.0 ) <<std::endl;
	std::cout <<"=> getInterval( 2.4 ) = " << f.getInterval( 2.4 ) <<std::endl;
	std::cout <<"=> getInterval( 13.0 ) = " << f.getInterval( 13.0 ) <<std::endl;

		// getEnergy( distance )
	{
		double dist = 0.1;
		while (dist < 5.0) {
			std::cout <<"->getEnergy( A, C, "<<dist<<" )" 
						<<f.getEnergy( alph.getElement("A"), alph.getElement("C"), dist )<<std::endl;
			dist += 1.0;
		}
	}
		// getEnergy( DblPoint )
	{
		biu::DblPoint p1(1.2,3.0,1.2);
		biu::DblPoint p2(4.3,2.0,1.8);
		std::cout <<"=> new DblPoint p1 = " << p1 <<std::endl;
		std::cout <<"=> new DblPoint p2 = " << p2 <<std::endl;
		std::cout <<"=> distance of p1 and p2 = " << p1.distance(p2) <<std::endl;
		std::cout <<"->getEnergy( A, C, p1, p2 ) = " 
					<<f.getEnergy( alph.getElement("A"), alph.getElement("C"), p1, p2 )<<std::endl;
	}

		// getEnergy( IntPoint )
	{
		biu::IntPoint p1(1,3,1);
		biu::IntPoint p2(4,2,1);
		std::cout <<"=> new IntPoint p1 = " << p1 <<std::endl;
		std::cout <<"=> new IntPoint p2 = " << p2 <<std::endl;
		std::cout <<"=> distance of p1 and p2 = " << p1.distance(p2) <<std::endl;
		std::cout <<"->getEnergy( A, C, p1, p2 ) = " 
					<<f.getEnergy( alph.getElement("A"), alph.getElement("C"), p1, p2 )<<std::endl;
	}

		// ==
	std::cout	<<"-> (f == f) = " 
				<<(f == f?"true":"false")
				<<std::endl;
		// !=
	std::cout	<<"-> (f != f) = " 
				<<(f != f?"true":"false")
				<<std::endl;

	  // make a copy
	std::cout <<"=> copy f : biu::IntervalEnergyFunction f2 = f" <<std::endl;
	biu::IntervalEnergyFunction f2 = f;
	
		// ==
	std::cout	<<"-> (f == f2) = " 
				<<(f == f2?"true":"false")
				<<std::endl;
		// !=
	std::cout	<<"-> (f != f2) = " 
				<<(f != f2?"true":"false")
				<<std::endl;

	std::cout <<"=> addInterval( E3, 15.0 ) = " <<f.addInterval(E3,15.0) <<std::endl;
		// ==
	std::cout	<<"-> (f == f2) = " 
				<<(f == f2?"true":"false")
				<<std::endl;
		// !=
	std::cout	<<"-> (f != f2) = " 
				<<(f != f2?"true":"false")
				<<std::endl;


	std::cout	<<"============= END OF TEST ================\n" <<std::endl;
	
	return 0;
}
