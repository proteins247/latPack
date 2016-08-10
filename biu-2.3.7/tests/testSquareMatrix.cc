
#include "biu/SquareMatrix.hh"
#include <iostream>

int main(int argc, char** argv) {
	typedef biu::SquareMatrix<int,3> SqMat;

	SqMat matrix;


	std::cout <<"matrix dim = " <<matrix.numRows() <<"x" <<matrix.numColumns() <<"\n" <<matrix <<std::endl;

	typedef biu::SquareMatrix<int, 3> SM3;

	int m1[3][3] = {{11,12,13}, {21,22,23},{31,32,33}};
	SM3 mat1 = SM3::createSquareMatrix(m1);
	int m2[3][3] = {{1,2,3}, {1,2,3},{1,2,3}};
	SM3 mat2 = SM3::createSquareMatrix(m2);
	SM3 mat3 = mat1-mat2;

	biu::Matrix<int> mat4(3,3,2);


	std::cout <<"squarematrix \n" <<mat1 <<"\n" <<mat2 <<"\n" <<mat3 <<std::endl;
	std::cout <<"prod \n" <<(mat2*mat4) <<std::endl;
	SqMat m5 = mat2;
	std::cout <<" m5 \n" <<m5 <<"\n" <<mat4;
	m5 *= mat4;
	std::cout <<" m5 \n" <<m5 <<std::endl;

	return 0;
}
