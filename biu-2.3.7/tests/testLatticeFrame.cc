

#include <iostream>
#include <biu/LatticeFrame.hh>
#include <biu/LatticeDescriptorCUB.hh>


int main(int argc, char **argv) {

	using namespace biu;

	std::cout	<<"========== BEGIN OF TEST =================\n"
				<<" test class : LatticeFrame\n"
				<<"==========================================" <<std::endl;
				
		// erzeugen	
	std::cout <<"=> creating biu::LatticeDescriptorCUB() = ld" <<std::endl;
	LatticeDescriptorCUB ld;
	std::cout <<"=> creating biu::LatticeFrame( ld, 4)" <<std::endl;
	LatticeFrame lf(&ld, 4);
	
		// getFrameSize()
	std::cout	<<"->getFrameSize()  = " << lf.getFrameSize()
				<<std::endl;

		// isInFrame() == false
		IntPoint p(3,4,5);
	std::cout	<<"-> IntPoint = " <<p <<" = p" <<std::endl;
	std::cout	<<"->isInFrame(p)  = " << (lf.isInFrame(p)?"true":"false")
				<<std::endl;

		// setFrameSize()
	std::cout	<<"->setFrameSize(7)" <<std::endl;
		lf.setFrameSize(7);							// setFrameSize()

		// getFrameSize()
	std::cout	<<"->getFrameSize()  = " << lf.getFrameSize()
				<<std::endl;

		// isInFrame() == true
	std::cout	<<"->isInFrame(p)  = " << (lf.isInFrame(p)?"true":"false")
				<<std::endl;

		// getPoint( getIndex( p )) == p
	std::cout	<<"->(getPoint(getIndex(p)) == p)  = "
				<< (lf.getPoint(lf.getIndex(p)) == p?"true":"false")
				<<std::endl;

		// getCenter()
	std::cout	<<"->getCenter()  = " << lf.getCenter()
				<<std::endl;

		// getIndexedNeighborhood()
	std::cout	<<"->getIndexedNeighborhood().size()  = " 
				<< lf.getIndexedNeighborhood().size()
				<<std::endl;
	
	std::cout	<<"->getIndexedNeighborhood() : ";
	for (biu::LatticeFrame::index_set::const_iterator it=lf.getIndexedNeighborhood().begin();
		it != lf.getIndexedNeighborhood().end(); it++ )
	{
		std::cout <<*it <<",";
	}
	std::cout <<std::endl;
		
		// == & !=
	std::cout	<<"-> (*this == *this)   = " << (lf==lf?"true":"false")
				<<std::endl;
	std::cout	<<"-> (*this != *this)   = " << (lf!=lf?"true":"false")
				<<std::endl;


	std::cout	<<"============= END OF TEST ================\n" <<std::endl;
	
	return 0;
}
