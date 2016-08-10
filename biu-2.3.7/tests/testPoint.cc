

#include <iostream>
#include <biu/Point.hh>


int main(int argc, char **argv) {

	using namespace biu;

	std::cout	<<"========== BEGIN OF TEST =================\n"
				<<" test template class : Point3D\n"
				<<"==========================================" <<std::endl;
				
		// erzeugen	
		IntPoint p1(1,2,3);
	std::cout <<"=> creating biu::IntPoint("<<p1<<") = p1" <<std::endl;
		IntPoint p2(1,0,-1);
	std::cout <<"=> creating biu::IntPoint("<<p2<<") = p2" <<std::endl;

	std::cout <<std::endl;
		// == true
	std::cout	<<"->(p1 == p1)  = "
				<<(p1 == p1?"true":"false")
				<<std::endl;

		// == false
	std::cout	<<"->(p1 == p2)  = "
				<<(p1 == p2?"true":"false")
				<<std::endl;

		// != false
	std::cout	<<"->(p1 != p1)  = "
				<<(p1 != p1?"true":"false")
				<<std::endl;

		// != true
	std::cout	<<"->(p1 != p2)  = "
				<<(p1 != p2?"true":"false")
				<<std::endl;

		// +
	std::cout	<<"->(p1 +  p2)  = "
				<<(p1 + p2)
				<<std::endl;

		// -
	std::cout	<<"->(p1 -  p2)  = "
				<<(p2 - p1)
				<<std::endl;

		// +=
	std::cout	<<"->(p1 += 1)   = "
				<<(p1 += 1)
				<<std::endl;

		// -=
	std::cout	<<"->(p1 -= 1)   = "
				<<(p1 -= 1)
				<<std::endl;

		// *=
	std::cout	<<"->(p1 *= 2)   = "
				<<(p1 *= 2)
				<<std::endl;

		// /=
	std::cout	<<"->(p1 /= 2)   = "
				<<(p1 /= 2)
				<<std::endl;

		// =
		IntPoint p3 = p1;
	std::cout	<<"->(p3 = p1)   = "
				<<p3
				<<std::endl;

		// < == false
	std::cout	<<"->(p1 < p2)   = "
				<<(p1<p2?"true":"false")
				<<std::endl;

		// > == true
	std::cout	<<"->(p1 > p2)   = "
				<<(p1>p2?"true":"false")
				<<std::endl;

		// < == true
	std::cout	<<"->(p2 < p1)   = "
				<<(p2<p1?"true":"false")
				<<std::endl;

		// > == false
	std::cout	<<"->(p2 > p1)   = "
				<<(p2>p1?"true":"false")
				<<std::endl;

		// += point
	std::cout	<<"->(p2 += IntPoint(0,1,1))  = "
				<<(p2 += IntPoint(0,1,1))
				<<std::endl;

		// -= point
	std::cout	<<"->(p2 -= IntPoint(0,1,0))  = "
				<<(p2 -= IntPoint(0,1,0))
				<<std::endl;

		// vectorLength()
	std::cout	<<"->p2.vectorLength()  = "
				<<p2.vectorLength()
				<<std::endl;

		// distance()
	std::cout	<<"->p2.distance(IntPoint(2,0,0))  = "
				<<p2.distance(IntPoint(2,0,0))
				<<std::endl;

		// isEven = false
	std::cout	<<"->p2.isEven()  = "
				<<(p2.isEven()?"true":"false")
				<<std::endl;

		// isEven = true
	std::cout	<<"->p1.isEven()  = "
				<<(p1.isEven()?"true":"false")
				<<std::endl;

		// unary -
	std::cout	<<"->(-p1)  = "
				<<-p1
				<<std::endl;


	std::cout	<<"\n============= END OF TEST ================\n" <<std::endl;
	
	return 0;
}
