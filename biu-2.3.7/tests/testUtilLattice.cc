
#include <iostream>
#include <biu/LatticeDescriptorCUB.hh>
#include <biu/LatticeModel.hh>
#include <biu/util/point.hh>


int
main() {
	using namespace biu;


	std::cout	<<"========== BEGIN OF TEST =================\n"
				<<" testing : lattice utilities\n"
				<<"==========================================" <<std::endl;
	
	biu::IPointSet s;
	biu::IPointSet p;
	
	std::cout <<"  s = " <<s <<std::endl;
	std::cout <<"  p = " <<p <<std::endl;
	std::cout <<"  isSubSet(s,p) = " <<(isSubset<int>(s,p)?"true":"false") <<std::endl;
	
	s.insert(biu::IntPoint(0,0,0));
	s.insert(biu::IntPoint(0,1,0));
	s.insert(biu::IntPoint(0,1,1));
	s.insert(biu::IntPoint(0,2,0));
	s.insert(biu::IntPoint(0,0,1));
	
	p.insert(biu::IntPoint(2,2,2));
	p.insert(biu::IntPoint(2,1,2));

	std::cout <<"  s = " <<s <<std::endl;
	std::cout <<"  p = " <<p <<std::endl;
	std::cout <<"  isSubSet(s,p) = " <<(isSubset<int>(s,p)?"true":"false") <<std::endl;
	std::cout <<"  isSubSet(p,s) = " <<(isSubset<int>(p,s)?"true":"false") <<std::endl;

	p.insert(biu::IntPoint(2,1,4));
	std::cout <<"  s = " <<s <<std::endl;
	std::cout <<"  p = " <<p <<std::endl;
	std::cout <<"  isSubSet(s,p) = " <<(isSubset<int>(s,p)?"true":"false") <<std::endl;

	std::cout	<<"============= END OF TEST ================\n" <<std::endl;
	
	return 0;
}
