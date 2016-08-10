
#include <biu/LimitedHash.hh>

#include <string>
#include <iostream>
#include <set>
#include <iterator>
#include <algorithm>

int main(int argc, char** argv) {
	

	std::cout	<<"========== BEGIN OF TEST =================\n"
				<<" test template class : LimitedHash\n"
				<<"==========================================" <<std::endl;
				
	const size_t size = 6; 
	std::cout <<"\n create biu::LimitedHash< size_t, size_t > hash("<<size<<")"
			<<std::endl;
	biu::LimitedHash< size_t, size_t > hash(size, 2);
	typedef biu::LimitedHash< size_t, size_t >::iterator it_type;
	
	std::cout	<<"\n hash.max_size() = "
				<<hash.max_size()
				<<std::endl;
	
	std::cout	<<std::endl;
	
	std::set<size_t> output;
	
	//////////////////////////////////////////////////////////////////////////
	std::cout <<" hash [begin,end] = ";
	for (it_type it=hash.begin(); it!=hash.end(); it++) {
		output.insert(it->first);
	}
	std::copy(output.begin(), output.end(), std::ostream_iterator<int>(std::cout, ", "));
	std::cout <<"\n" <<std::endl;
	output.clear();

	std::cout <<" hash.insert( 1, 2 ) " <<std::endl;
	hash.insert( std::make_pair(1,2));
	std::cout <<" hash [begin,end] = ";
	for (it_type it=hash.begin(); it!=hash.end(); it++) {
		output.insert(it->first);
	}
	std::copy(output.begin(), output.end(), std::ostream_iterator<int>(std::cout, ", "));
	std::cout <<"\n" <<std::endl;
	output.clear();
	
	//////////////////////////////////////////////////////////////////////////
	size_t i = 0;
	for (i=1; i<=2; i++) {
		std::cout <<" hash[" <<i <<"] = " <<i <<std::endl;
		hash[i] = i;
	}
	std::cout <<" hash.size() = " <<hash.size() <<std::endl;
	std::cout <<" hash [begin,end] = ";
	for (it_type it=hash.begin(); it!=hash.end(); it++) {
		output.insert(it->first);
	}
	std::copy(output.begin(), output.end(), std::ostream_iterator<int>(std::cout, ", "));
	std::cout <<"\n" <<std::endl;
	output.clear();
	
	
	//////////////////////////////////////////////////////////////////////////
	for (i=1; i<=6; i++) {
		std::cout <<" hash[" <<i <<"] = " <<i <<std::endl;
		hash[i] = i;
	}
	std::cout <<" hash.size() = " <<hash.size() <<std::endl;
	std::cout <<" hash [begin,end] = ";
	for (it_type it=hash.begin(); it!=hash.end(); ++it) {
		output.insert((*it).first);
	}
	std::copy(output.begin(), output.end(), std::ostream_iterator<int>(std::cout, ", "));
	std::cout <<"\n" <<std::endl;
	output.clear();

	
	//////////////////////////////////////////////////////////////////////////
	i=7;
	std::cout <<" hash[" <<i <<"] = " <<i <<std::endl;
	hash[i] = i;
	std::cout <<" hash.size() = " <<hash.size() <<std::endl;
	std::cout <<" hash [begin,end] = ";
	for (it_type it=hash.begin(); it!=hash.end(); ++it) {
		output.insert((*it).first);
	}
	std::copy(output.begin(), output.end(), std::ostream_iterator<int>(std::cout, ", "));
	std::cout <<"\n" <<std::endl;
	output.clear();

	return 0;
}

