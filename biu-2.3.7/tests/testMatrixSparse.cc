
#include "biu/MatrixSparse.hh"
#include <iostream>


	
int testMatrixSparseR_int() {
	

	std::cout	<<"========== BEGIN OF TEST =================\n"
				<<" test template class : MatrixSparseR<int>\n"
				<<"==========================================" <<std::endl;

	typedef biu::MatrixSparseR<int> mathInt;
	
	mathInt m1 = mathInt(4,4,0), m2 = mathInt(4,2,0);
	
	std::cout << "matrix m1 = \n" <<m1 <<"matrix m2 = \n" <<m2 <<std::endl;
	
	std::cout <<"matrix m3(m2) = "<<std::endl;
	mathInt m3(m2);
	std::cout <<m3 <<std::endl;
	std::cout <<"matrix m3 = m1 : "<<std::endl;
	m3 = m1;
	std::cout <<m3 <<std::endl;
	
	std::cout <<"m2.resize(3,3,1)" <<std::endl;
	m2.resize(3,3,1);
	std::cout <<"matrix m2 = \n" <<m2 <<std::endl;

	std::cout <<"m1 == m2 = " <<(m1 == m2 ? "true" : "false") <<std::endl;
	
	std::cout <<"m2.numRows() = " <<m2.numRows() <<std::endl;
	std::cout <<"m2.numColumns() = " <<m2.numColumns() <<std::endl;
	std::cout <<"m2.resize(4,4)" <<std::endl;
	m2.resize(4,4);
	std::cout <<"m2.numRows() = " <<m2.numRows() <<std::endl;
	std::cout <<"m2.numColumns() = " <<m2.numColumns() <<std::endl;
	std::cout <<"m1 == m2 = " <<(m1 == m2 ? "true" : "false") <<std::endl;
	
	std::cout <<"m1.getDefaultValue() = " <<m1.getDefaultValue() <<std::endl;
	std::cout <<"m1.setDefaultValue(1)" <<std::endl;
	m1.setDefaultValue(1);
	std::cout <<"m1.getDefaultValue() = " <<m1.getDefaultValue() <<std::endl;
	std::cout <<"m1 == m2 = " <<(m1 == m2 ? "true" : "false") <<std::endl;
	std::cout <<"m2 == m2 = " <<(m2 == m2 ? "true" : "false") <<std::endl;
	
	std::cout <<"m1.exists(2,1) = " <<(m1.exists(2,1) ? "true" : "false") <<std::endl;
	std::cout <<"m1.atConst(2,1) = " <<m1.atConst(2,1) <<std::endl;
	std::cout <<"m1.exists(2,1) = " <<(m1.exists(2,1) ? "true" : "false") <<std::endl;
	std::cout <<"(const m1).at(2,1) = " <<((const mathInt)m1).at(2,1) <<std::endl;
	std::cout <<"m1.exists(2,1) = " <<(m1.exists(2,1) ? "true" : "false") <<std::endl;
	std::cout <<"m1.at(2,1) += 3" <<std::endl;
	m1.at(2,1) += 3;
	std::cout <<"m1.exists(2,1) = " <<(m1.exists(2,1) ? "true" : "false") <<std::endl;
	std::cout <<"m1.at(2,1) = " <<m1.at(2,1) <<std::endl;
	
	std::cout <<"m1.rowVec(..)"<<std::endl;
	for (size_t r=0; r<m1.numRows(); r++) {
		const std::vector<int> row = m1.rowVec(r);
		for (size_t i=0; i<row.size(); i++) {
			std::cout <<row.at(i) <<' ';
		}
		std::cout <<std::endl;
	}
	std::cout <<"m1.rowValues(..)"<<std::endl;
	for (size_t r=0; r<m1.numRows(); r++) {
		const mathInt::EntryMap row = m1.rowValues(r);
		std::cout <<r <<" : ";
		for (mathInt::EntryMap::const_iterator it=row.begin(); it!=row.end(); it++) {
			std::cout <<it->first <<"=" <<it->second <<", ";
		}
		std::cout <<std::endl;
	}
	std::cout <<"m1.columnVec(..)"<<std::endl;
	std::vector< std::vector<int> > cols(m1.numColumns());
	for (size_t c=0; c<m1.numColumns(); c++) {
		cols[c] = m1.columnVec(c);
	}
	for (size_t r=0; r<m1.numRows(); r++) {
		for (size_t c=0; c<m1.numColumns(); c++) {
			std::cout <<cols.at(c).at(r) <<' ';
		}
		std::cout <<std::endl;
	}
	std::cout <<"m1.columnValues(..)"<<std::endl;
	for (size_t c=0; c<m1.numColumns(); c++) {
		const mathInt::EntryMap col = m1.columnValues(c);
		std::cout <<c <<" : ";
		for (mathInt::EntryMap::const_iterator it=col.begin(); it!=col.end(); it++) {
			std::cout <<it->first <<"=" <<it->second <<", ";
		}
		std::cout <<std::endl;
	}

	std::cout	<<"\n========== END OF TEST =================\n"
				<<std::endl;
	
	return 0;
}

	
int testMatrixSparseC_int() {

	

	std::cout	<<"========== BEGIN OF TEST =================\n"
				<<" test template class : MatrixSparseC<int>\n"
				<<"==========================================" <<std::endl;

	typedef biu::MatrixSparseC<int> mathInt;
	
	mathInt m1 = mathInt(4,4,0), m2 = mathInt(4,2,0);
	
	std::cout << "matrix m1 = \n" <<m1 <<"matrix m2 = \n" <<m2 <<std::endl;
	
	std::cout <<"matrix m3(m2) = "<<std::endl;
	mathInt m3(m2);
	std::cout <<m3 <<std::endl;
	std::cout <<"matrix m3 = m1 : "<<std::endl;
	m3 = m1;
	std::cout <<m3 <<std::endl;
	
	std::cout <<"m2.resize(3,3,1)" <<std::endl;
	m2.resize(3,3,1);
	std::cout <<"matrix m2 = \n" <<m2 <<std::endl;

	std::cout <<"m1 == m2 = " <<(m1 == m2 ? "true" : "false") <<std::endl;
	
	std::cout <<"m2.numRows() = " <<m2.numRows() <<std::endl;
	std::cout <<"m2.numColumns() = " <<m2.numColumns() <<std::endl;
	std::cout <<"m2.resize(4,4)" <<std::endl;
	m2.resize(4,4);
	std::cout <<"m2.numRows() = " <<m2.numRows() <<std::endl;
	std::cout <<"m2.numColumns() = " <<m2.numColumns() <<std::endl;
	std::cout <<"m1 == m2 = " <<(m1 == m2 ? "true" : "false") <<std::endl;
	
	std::cout <<"m1.getDefaultValue() = " <<m1.getDefaultValue() <<std::endl;
	std::cout <<"m1.setDefaultValue(1)" <<std::endl;
	m1.setDefaultValue(1);
	std::cout <<"m1.getDefaultValue() = " <<m1.getDefaultValue() <<std::endl;
	std::cout <<"m1 == m2 = " <<(m1 == m2 ? "true" : "false") <<std::endl;
	std::cout <<"m2 == m2 = " <<(m2 == m2 ? "true" : "false") <<std::endl;
	
	std::cout <<"m1.exists(2,1) = " <<(m1.exists(2,1) ? "true" : "false") <<std::endl;
	std::cout <<"m1.atConst(2,1) = " <<m1.atConst(2,1) <<std::endl;
	std::cout <<"m1.exists(2,1) = " <<(m1.exists(2,1) ? "true" : "false") <<std::endl;
	std::cout <<"(const m1).at(2,1) = " <<((const mathInt)m1).at(2,1) <<std::endl;
	std::cout <<"m1.exists(2,1) = " <<(m1.exists(2,1) ? "true" : "false") <<std::endl;
	std::cout <<"m1.at(2,1) += 3" <<std::endl;
	m1.at(2,1) += 3;
	std::cout <<"m1.exists(2,1) = " <<(m1.exists(2,1) ? "true" : "false") <<std::endl;
	std::cout <<"m1.at(2,1) = " <<m1.at(2,1) <<std::endl;
	
	std::cout <<"m1.rowVec(..)"<<std::endl;
	for (size_t r=0; r<m1.numRows(); r++) {
		const std::vector<int> row = m1.rowVec(r);
		for (size_t i=0; i<row.size(); i++) {
			std::cout <<row.at(i) <<' ';
		}
		std::cout <<std::endl;
	}
	std::cout <<"m1.rowValues(..)"<<std::endl;
	for (size_t r=0; r<m1.numRows(); r++) {
		const mathInt::EntryMap row = m1.rowValues(r);
		std::cout <<r <<" : ";
		for (mathInt::EntryMap::const_iterator it=row.begin(); it!=row.end(); it++) {
			std::cout <<it->first <<"=" <<it->second <<", ";
		}
		std::cout <<std::endl;
	}
	std::cout <<"m1.columnVec(..)"<<std::endl;
	std::vector< std::vector<int> > cols(m1.numColumns());
	for (size_t c=0; c<m1.numColumns(); c++) {
		cols[c] = m1.columnVec(c);
	}
	for (size_t r=0; r<m1.numRows(); r++) {
		for (size_t c=0; c<m1.numColumns(); c++) {
			std::cout <<cols.at(c).at(r) <<' ';
		}
		std::cout <<std::endl;
	}
	std::cout <<"m1.columnValues(..)"<<std::endl;
	for (size_t c=0; c<m1.numColumns(); c++) {
		const mathInt::EntryMap col = m1.columnValues(c);
		std::cout <<c <<" : ";
		for (mathInt::EntryMap::const_iterator it=col.begin(); it!=col.end(); it++) {
			std::cout <<it->first <<"=" <<it->second <<", ";
		}
		std::cout <<std::endl;
	}
	
	std::cout	<<"\n========== END OF TEST =================\n"
				<<std::endl;
	
	return 0;
}

int main(int argc, char ** argv) {

	testMatrixSparseR_int();
	testMatrixSparseC_int();
}
